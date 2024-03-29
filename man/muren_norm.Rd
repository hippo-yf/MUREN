% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm.R
\name{muren_norm}
\alias{muren_norm}
\title{MUREN: a Robust and Multi-reference Approach of RNA-seq
Transcript Normalization}
\usage{
muren_norm(
  reads,
  refs = "saturated",
  pairwise_method = "lts",
  single_param = TRUE,
  res_return = "counts",
  filter_gene = TRUE,
  trim = 10,
  maxiter = 70,
  workers = 2,
  ...
)
}
\arguments{
\item{reads}{data.frame or matrix. A tabular counts of gene/transcript X sample}

\item{refs}{reference samples by name characters or integers of sample orders.
Or, 'saturated' selects all samples as references.}

\item{pairwise_method}{string. \code{lts} (Least Trimmed Squares regression)
or \code{mode}}

\item{single_param}{logical. single parameter form (scaling) or double paramter
form (non-linear)}

\item{res_return}{type of values returned. \code{counts} (normalized counts),
\code{log_counts} (log2 (normalized counts + 1)), or \code{scaling_coeff} (scaling
coefficients (counterpart of library size), valid when \code{single_param = TRUE})}

\item{filter_gene}{logical. whether filter genes with rare raw counts, the filter
strategy simply filters genes with the maximum count less than \code{trim}}

\item{trim}{numeric. filter genes with counts less than \code{trim} in
all samples}

\item{maxiter}{integer. the maximum number of iterations
in the median polish.}

\item{workers}{integer. the number of nodes in the parallel cluster (\code{snow}).}

\item{...}{additional parameters passing on to \code{pairwise_method}}
}
\value{
depends \code{res_return}
}
\description{
MUltiple REference Normalizer (MUREN) of RNA-seq counts.
}
\details{
The \code{reads} should be a tabular type that contains each sample
in the columns. Columns of non-numeric annotations are kept in the
returned values. Notice than all the numeric columns all treated as samples.

You can specify reference samples with their names (column names of \code{reads}),
or the order indices, e.g., \code{1:3} means the first 3 samples. Or, take all
samples as references (\code{refs = 'saturated'}).

Single parameter form (\code{single_param = TRUE}) (default and recommended)
is a scaling normalization, in this case, \code{lts} (Least Trimmed Squares regression)
or \code{mode} are applicable. Double parameter form
(\code{single_param = FALSE}) (slower) is a non-linear (power function)
normalization.

The returned scaling coeff (\code{res_return = 'scaling_coeff'}), which is
counterpart of library size, should divide
the raw read counts to get normalized counts. Never re-adjust it by total sample
counts or other numbers

Parallel computation is implemented with \code{doSNOW} and its dependencies.
}
\examples{
# A rat toxicogenomics RNA-seq data
# from Munro et al (2014)

data(rat_tox_thi)

# normalized counts (raw scale)
thi_norm = muren_norm(rat_tox_thi)

# scaling coefficient is a counterpart of library size
(thi_coeff = muren_norm(rat_tox_thi, res_return = 'scaling_coeff'))

# normalize manually
thi_norm2 = cbind(rat_tox_thi[1],
              as.matrix(rat_tox_thi[-1]) \%*\% diag(1/thi_coeff)
            )

}
\references{
Munro SA, Lund SP, Pine PS, et al. Assessing technical
performance in differential gene expression experiments with
external spike-in RNA control ratio mixtures. Nat Commun.
2014;5:5125. doi:10.1038/ncomms6125
}
