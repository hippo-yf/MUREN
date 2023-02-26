

#' @importFrom magrittr %>%
#' @importFrom magrittr not
#' @import matrixStats
#'
#'


# constants
TOL = .Machine$double.eps ^ 0.5


is.wholenumber <-
  function(x, tol = TOL)
    abs(x - round(x)) < tol

lg <- function(x) log2(1 + x)
ep <- function(x) 2**x - 1

filter_gene_l <- function(reads, trim){

  rowMaxs(reads) >= trim
}


# single parameter regression

reg_sp <- function(
  # raw_s_k,      # un-filtered sample k
  s_k,          # objective sample
  s_r,          # ref
  # pairwise_method,   # robust regression
  # t,            # lg(trim)
  ...
) {


  ltsreg(s_r - s_k ~ 1, ...)$coefficients[1]
}


# mode
mode_sp <- function(s_k, s_r, ...){
  d = density(s_r - s_k, ...)
  d$x[which.max(d$y)]
}

# pairwise regressions return fitted values
# it is expansive
# except the sp case


# double parameter regression

reg_dp <- function(s_k, s_r, ...) {

  ltsreg(s_r ~ s_k, ...)$fitted.values

}

# tasks for sample k
get_tasks <- function(k, reg_wapper, refs) {
  # refs = get_refs(design, k)
  parse(
    text = paste(
      reg_wapper,
      "(log_raw_reads_mx[, ",
      k ,
      "],log_raw_reads_mx[, ",
      refs,
      '], ...)',
      collapse = '\n',
      sep = ''
    )
  )
}


polish_coeff <-
  function(
           fitted_n,       # scaling coeff
           n_exp,          # number of experiments
           locations,      # line to matrix
           unused_refs,    # samples don't act as refs
           maxiter         # max iteration
  ) {
    # for nth-gene

    # rearrangement
    rs_mx = matrix(NA, nrow = n_exp, ncol = n_exp)
    rs_mx[locations] = fitted_n

    # pad the diagonals
    # rs_mx[seq.int(1, n_exp ^ 2, n_exp + 1)] = raw_n


    # remove unused refs (rows)
    if (length(unused_refs) > 0)
      rs_mx = rs_mx[-unused_refs,]

    # medpolish
    m = medpolish(rs_mx,
                  na.rm = TRUE,
                  trace.iter = FALSE,
                  maxiter = maxiter)

    # return polished sample (column) effects
    m$overall + m$col

  }

polish_one_gene <-
  # function(n) { # n-gene
  function(
    # raw_n,         # raw (log) expr of nth-gene
    fitted_n,       # fitted expr in the regression
    n_exp,          # number of experiments
    locations,      # line to matrix
    unused_refs,    # samples don't act as refs
    maxiter         # max iteration
  ) {

  # rearrangement
  rs_mx = matrix(NA, nrow = n_exp, ncol = n_exp)
  rs_mx[locations] = fitted_n

  # pad the diagonals
  # rs_mx[seq.int(1, n_exp ^ 2, n_exp + 1)] = log_reads_mx[n,]

  # remove unused refs (rows)
  if (length(unused_refs) > 0)
    rs_mx = rs_mx[-unused_refs,]

  # medpolish
  m = medpolish(rs_mx, na.rm = TRUE, trace.iter = FALSE, maxiter = maxiter)

  # return polished sample (column) effects
  m$overall + m$col

}
