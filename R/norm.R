#' An Accurate RNA-seq Transcript Quantification Method
#' Preserving Biological Asymmetric Differentiation
#'
#'
#' MUltiple REference Normalizer (MUREN) of RNA-seq counts.
#'
#' @usage
#' muren_norm(reads, refs = 'saturated', single_param = TRUE,
#'   pairwise_method = 'lts', maxiter = 70, workers = 2, ...)
#'
#' @param reads data.frame or matrix. A tabular counts of gene/transcript

#' @param refs name characters or integers of sample orders of reference samples.
#'   Or, 'saturated' selects all samples as references.
#'
#' @param pairwise_method string. \code{lts} or \code{mode}
#' @param single_param logical. single parameter form (scaling) or double paramter
#' form (non-linear)
#' @param res_return type of values returned. \code{counts} (normalized counts),
#'   \code{log_counts} (log2 (normalized counts + 1)), or \code{library_size} (scaling
#'   coefficients, valid when \code{single_param = TRUE})
#' @param trim numeric. trim genes with counts less than \code{trim} in
#'   all samples
#' @param maxiter ingeter. the maximum number of iterations
#'   in the median polish.
#' @param workers integer. the number of nodes in the parallel cluster.
#' @param ... additional parameters passing on to \code{pairwise_method}
#'
#' @details
#' The \code{reads} should be a tabular format that contains each sample
#' in the columns. Columns of non-numeric annotations are kept in the
#' returned values. Note than all the numeric columns all treated as samples.
#'
#' You can specify reference samples with their names (column names of \code{reads}),
#' or the order indices, e.g., \code{1:3} means the first 3 samples. Or, take all
#' samples as references (\code{refs = 'saturated'}).
#'
#' Single parameter form (\code{single_param = TRUE}) (default and recommended)
#' is a scaling normalization, in this case, \code{lts} (Least Trimmed Squares regression)
#' or \code{mode} are applicable. Double parameter form
#' (\code{single_param = FALSE}) (slower) is a non-linear (power function)
#' normalization.
#'
#' The returned library_size (\code{res_return = 'library_size'}) should divide
#' the raw read counts to get normalized counts. Never re-adjust it by total sample
#' counts or other qualities.
#'
#' Parallel computation is implemented with \code{doSNOW} and its dependencies.
#'
#' @return depends \code{res_return}
#'
#' @examples
#' ## A rat toxicogenomics RNA-seq data
#'
#' data(rat_tox_thi)
#'
#' # normalzied counts
#' thi_norm = muren_norm(rat_tox_thi)
#'
#' # library sizes (scaling coefficients)
#' (thi_coeff = muren_norm(rat_tox_thi, res_return = 'library_size'))
#'
#'
#' @include utils.R
#' @import doSNOW
#' @import MASS
#'
#' @import magrittr
#' @importFrom assertthat assert_that
#' @rdname muren_norm
#' @name muren_norm




#' @export
muren_norm <- function(reads,
                       # design = NULL,
                       refs = 'saturated',
                       pairwise_method = "lts",
                       single_param = TRUE,
                       res_return = 'counts',
                       trim = 10,
                       maxiter = 70,
                       workers = 2,
                       ...
                       ) {




  ##  check params

  assert_that(is.logical(single_param))

  assert_that(is.character(pairwise_method))
  assert_that(any(pairwise_method %in% c('lts', 'mode')))

  assert_that(maxiter > 0 & is.wholenumber(maxiter))
  assert_that(workers > 0 & is.wholenumber(workers))

  # if(!is.null(design)) assert_that(is.data.frame(design))
  assert_that(is.data.frame(reads) | is.matrix(reads))
  assert_that(trim >= 0 )


  # too few samples
  assert_that(nrow(reads) > 1, msg = "Bad input data !")

  # filter the un-numeric cols
  raw_reads = reads

  log_raw_reads_mx =
    raw_reads %>%
    sapply(is.numeric) %>%
    extract(raw_reads, .) %>%
    as.matrix %>%
    lg

  log_raw_reads_df = as.data.frame(log_raw_reads_mx)


  reads %<>%
    sapply(is.numeric) %>%
    extract(reads, .) %>%
    as.matrix %>%
    # filter the low reads counts genes in param estimation
    extract(., rowMaxs(.) >= trim,)


  # sample ids in reads and design don't coincide
  # if(!is.null(design))
  #   assert_that (! anyDuplicated(colnames(reads)) &
  #                ! anyDuplicated(design[[1]]) &
  #                  ncol(reads) == nrow(design) &
  #                  all(colnames(reads) %in% design[[1]]),
  #                    msg = "Bad sample IDs in `reads` and/or `design`! ")


  REFSERR = "Bad specification of references !"

  # specify reference strategy

  if (is.character(refs)) {
    if (length(refs) == 1) {
      # if (!is.null(design)){
      #
      #   # first sample in each level conbination
      #   # sepcified by experimental design
      #   # override refs
      #   get_refs <- function(design, k) {
      #     design %>%
      #       extract(-1) %>%
      #       sapply(as.character) %>%
      #       apply(1, paste0, collapse = '') %>%
      #       duplicated %>%
      #       not %>%
      #       which %>%
      #       setdiff(k) %>%
      #       return
      #   }
      #
      # }

      if (refs == "saturated") {
        get_refs <- function(k) {
          reads %>%
            ncol %>%
            seq.int %>%
            is_in(k) %>%
            not %>%
            which %>%
            setdiff(k) %>%
            return
          }
      }

      # single ref specified by sample name
      else if (refs %in% colnames(reads)) {
        i = which(colnames(reads) %in% refs)
        get_refs <- function(k) {
          return(i)
        }
      }

      else {
        stop(REFSERR)
      }
    }

    # length(refs) != 1
    else{
      # refs is a string vector
      if (all(refs %in% colnames(reads))){
        i = which(colnames(reads) %in% refs)
        get_refs <- function(k) {
          return(setdiff(i, k))
        }
      }
      else{
        stop(REFSERR)
      }
    }

  }

  # is.character(refs) == FALSE

  else if (all(is.wholenumber(refs))) {

    # integer refs
    if (max(refs) <= ncol(reads)){

      # multiple integers
      if (length(refs) > 1){
        get_refs <- function(k) {
          return(setdiff(refs, k))
        }
      }
      else {

        # one integer
        get_refs <- function(k) {
          return(refs)
        }
      }

    }

    else {
      stop(REFSERR)
    }
  }

  # is.integer(refs) == FALSE
  else{
    stop(REFSERR)
  }



  # specify the regression method
  if (single_param){
    if(pairwise_method == 'mode')
      reg_wapper = 'mode_sp'
    else
      reg_wapper = 'reg_sp'
  }
  else
    reg_wapper = 'reg_dp'



  # log2(x + 1)
  # keep the data.frame and matrix forms

  reads = as.data.frame(reads)
  log_reads_mx = lg(as.matrix(reads))
  log_reads_df = as.data.frame(log_reads_mx)



  # number of experiments
  n_exp = ncol(reads)

  # number of genes
  n_gene = nrow(reads)
  N_gene = nrow(raw_reads)

  # refs for each sample
  refs = lapply(1:n_exp,
                function(k){
                  return(get_refs(k))
                }
  )

  # the used (as refs) samples
  used_refs = sort(unique(unlist(refs)))
  unused_refs = (1:n_exp)[!(1:n_exp %in% used_refs)]



  ## Step 1: regressions with each ref

  # task quene for robust regression (ltsreg)
  task_quene = unlist(lapply(1:n_exp, function(k)
    return(
      get_tasks(k, reg_wapper, refs[[k]])
    )))

  # regressions in parallel
  cl <- makeCluster(workers, type = "SOCK")
  registerDoSNOW(cl)

  coeffs <-  foreach(
    task = iter(task_quene),
    .packages = c("MASS"),
    .combine = cbind,
    .export = c("log_reads_df", "reg_sp", "reg_dp", "mode_sp")
  ) %dopar%
    eval(task)



  ## Step 2: Median polish for each gene

  # For each gene, the fitted reads of different refs and samples
  # lie in the same row of matrix fitted_reads.
  # A matrix of two factors 'ref X sample' (rs_mx) is needed for median polish.
  # The levels of refs are set to 1,2,...,n_exp for brief, and
  # then remove the unused refs.


  # Skip median polish if there is one ref
  if (refs %>% sapply(length) %>% max == 1){

    # transform to raw scale
    fitted_reads  = ep(fitted_reads)
    # 0 -> 0

    # fitted_reads[abs(log_raw_reads_mx) < TOL] = 0

    rownames(fitted_reads) = NULL
    colnames(fitted_reads) = colnames(reads)

    t = raw_reads %>%
      sapply(is.numeric) %>%
      not %>%
      extract(raw_reads, .) %>%
      cbind(fitted_reads)

    return(t)
  }

  # locations in rs_mx
  locations = rep(0, length(task_quene))
  k = 1
  for (i in 1:n_exp) {
    # i-sample
    for (j in refs[[i]]) {
      # j-ref
      locations[k] = j + (i - 1) * n_exp
      k = k + 1
    }
  }


  # coeffs = as.matrix(coeffs)

  # genewise polish when dp

  if (!single_param) {
    alpha_dp = matrix(coeffs[1,], byrow = T, nrow = length(used_refs))
    beta_dp = matrix(coeffs[2,], byrow = T, nrow = length(used_refs))

    polished_mx = foreach(
      n = 1:N_gene,
      .combine = rbind,
      .export = c("polish_one_gene")
    ) %dopar%
      polish_one_gene(log_raw_reads_mx[n,],

                      # pre-normalized log counts
                      as.vector(t(log_raw_reads_mx[n,]*beta_dp+alpha_dp)),
                      n_exp,
                      locations,
                      unused_refs,
                      maxiter)

    # 0 -> 0 when dp
    polished_mx[log_raw_reads_mx < TOL | polished_mx < 0] = 0

    # return
    if(res_return == 'counts'){
      polished_mx = ep(polished_mx)
    }


  }
  stopCluster(cl)


  if(single_param){
    coef_sp = polish_one_gene(1,
                    coeffs,
                    n_exp,
                    locations,
                    unused_refs,
                    maxiter)

    # return library size that be divided by raw counts
    coef_sp = 2^(as.vector(coef_sp))
    names(coef_sp) = colnames(reads)

    if(res_return == 'library_size'){

      return(1/coef_sp)
    }

    polished_mx =  ep(log_raw_reads_mx) %*% diag(coef_sp)

    if(res_return == 'log_counts'){
      polished_mx = lg(polished_mx)
    }

  }



  rownames(polished_mx) = NULL
  colnames(polished_mx) = colnames(reads)

  #

  ## add anno if raw_reads is data.frame

  if(is.data.frame(raw_reads)){
    i_count = sapply(raw_reads, is.numeric)
    if(!all(i_count)){
      res = raw_reads
      res[, i_count] = polished_mx

      return(res)
    }

  }

  return(as.data.frame(polished_mx))


}
