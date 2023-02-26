#' MUREN: a Robust and Multi-reference Approach of RNA-seq
#' Transcript Normalization
#'
#'
#' MUltiple REference Normalizer (MUREN) of RNA-seq counts.
#'
#' @usage
#' muren_norm(reads, refs = 'saturated', single_param = TRUE,
#'   pairwise_method = 'lts', maxiter = 70, workers = 2, ...)
#'
#' @param reads data.frame or matrix. A tabular counts of gene/transcript X sample

#' @param refs reference samples by name characters or integers of sample orders.
#'   Or, 'saturated' selects all samples as references.
#'
#' @param pairwise_method string. \code{lts} (Least Trimmed Squares regression)
#'   or \code{mode}
#' @param single_param logical. single parameter form (scaling) or double paramter
#' form (non-linear)
#' @param res_return type of values returned. \code{counts} (normalized counts),
#'   \code{log_counts} (log2 (normalized counts + 1)), or \code{scaling_coeff} (scaling
#'   coefficients, valid when \code{single_param = TRUE})
#' @param filter_gene logical. whether filter genes with rare raw counts, the filter
#'   strategy simply filters genes with the maximum count less than \code{trim}
#' @param trim numeric. filter genes with counts less than \code{trim} in
#'   all samples
#' @param maxiter integer. the maximum number of iterations
#'   in the median polish.
#' @param workers integer. the number of nodes in the parallel cluster (\code{snow}).
#' @param ... additional parameters passing on to \code{pairwise_method}
#'
#' @details
#' The \code{reads} should be a tabular type that contains each sample
#' in the columns. Columns of non-numeric annotations are kept in the
#' returned values. Notice than all the numeric columns all treated as samples.
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
#' The returned scaling coeff (\code{res_return = 'scaling_coeff'}) should divide
#' the raw read counts to get normalized counts. Never re-adjust it by total sample
#' counts or other quantities.
#'
#' Parallel computation is implemented with \code{doSNOW} and its dependencies.
#'
#' @return depends \code{res_return}
#'
#' @examples
#' # A rat toxicogenomics RNA-seq data
#' # from Munro et al (2014)
#'
#' data(rat_tox_thi)
#'
#' # normalized counts (raw scale)
#' thi_norm = muren_norm(rat_tox_thi)
#'
#' # scaling coefficients
#' (thi_coeff = muren_norm(rat_tox_thi, res_return = 'scaling_coeff'))
#'
#' @references
#' Munro SA, Lund SP, Pine PS, et al. Assessing technical
#' performance in differential gene expression experiments with
#' external spike-in RNA control ratio mixtures. Nat Commun.
#' 2014;5:5125. doi:10.1038/ncomms6125
#'
#' @include utils.R
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom MASS ltsreg
#' @importFrom matrixStats rowMaxs
#' @import magrittr
#' @importFrom assertthat assert_that
#' @importFrom iterators iter
#'
#' @rdname muren_norm
#' @name muren_norm



#' @export
muren_norm <- function(reads,
                       refs = 'saturated',
                       pairwise_method = "lts",
                       single_param = TRUE,
                       res_return = 'counts',
                       filter_gene = TRUE,
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

  # numeric cols
  i_sample = rep(TRUE, ncol(raw_reads))
  if(is.data.frame(raw_reads)){
    i_sample = sapply(raw_reads, is.numeric)
  }


  reads = raw_reads %>%
    extract(, i_sample) %>%
    as.matrix

  # filter the low reads counts genes in param estimation
  i_gene = filter_gene_l(reads, trim)

  if(filter_gene){
    reads %<>% extract(i_gene, )
  }

  log_raw_reads_mx = reads %>% lg

  # log_raw_reads_df = as.data.frame(log_raw_reads_mx)

  # #samples
  n_exp = ncol(reads)


  REFSERR = "Bad specification of references !"

  # specify reference strategy

  if (is.character(refs)) {
    if (length(refs) == 1) {

      if (refs == "saturated") {
        get_refs <- function(k) {
          1:n_exp
          # reads %>%
          #   ncol %>%
          #   seq.int %>%
          #   is_in(k) %>%
          #   not %>%
          #   which %>%
          #   # setdiff(k) %>%
          #   return
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
          # return(setdiff(i, k))
          return(i)
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
    if (max(refs) <= n_exp){

      # multiple integers
      if (length(refs) > 1){
        get_refs <- function(k) {
          # return(setdiff(refs, k))
          return(refs)
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


  # number of genes
  n_gene = nrow(reads)
  N_gene = nrow(raw_reads)

  # refs for each sample
  refs = lapply(1:n_exp, get_refs)

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

  res_pairwise <-  foreach(
    task = iter(task_quene),
    .packages = c("MASS"),
    .combine = cbind,
    .export = c("log_raw_reads_mx", "reg_sp", "reg_dp", "mode_sp")
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

  # gene-wise polish in dp and general (spline/smooth) cases

  if (!single_param) {
    # alpha_dp = matrix(coeffs[1,], byrow = T, nrow = length(used_refs))
    # beta_dp = matrix(coeffs[2,], byrow = T, nrow = length(used_refs))

    polished_mx = foreach(
      n = 1:n_gene,
      .combine = rbind,
      .export = c("polish_one_gene")
    ) %dopar%
      polish_one_gene(
        res_pairwise[n,], # pre-normalized log counts
        n_exp,
        locations,
        unused_refs,
        maxiter)

    # 0 -> 0 when dp
    polished_mx[log_raw_reads_mx < TOL | polished_mx < 0] = 0

    # return value
    if(res_return == 'counts'){
      polished_mx = ep(polished_mx)
    }


  }
  stopCluster(cl)


  if(single_param){
    coef_sp = polish_coeff(res_pairwise, # scaling coeffs
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

    polished_mx =  as.matrix(raw_reads[, i_sample]) %*% diag(coef_sp)

    if(res_return == 'log_counts'){
      polished_mx = lg(polished_mx)
    }

  }



  rownames(polished_mx) = NULL
  if(!is.null(rownames(raw_reads))){
    if(single_param) rownames(polished_mx) = rownames(raw_reads)
    else rownames(polished_mx) = rownames(raw_reads)[i_gene]
  }
  colnames(polished_mx) = colnames(reads)

  #

  ## add anno if raw_reads is data.frame

  if(is.data.frame(raw_reads)){
    i_count = sapply(raw_reads, is.numeric)
    if(!all(i_count)){
      res = raw_reads
      if(single_param){
        res[, i_count] = polished_mx
      }
      else{
        res[i_gene, i_count] = polished_mx
      }


      return(res)
    }

    return(as.data.frame(polished_mx))
  }

  return(polished_mx)


}
