# MUlti-REference Normalizer (MUREN)

MUREN: a Robust and Multi-reference Approach of RNA-seq Transcript Normalization

## Background

Normalization of RNA-seq expression data aims at identifying biological expression differentiation between samples by removing unwanted confounding factors. Explicitly or implicitly, the justification of normalization requires a set of housekeeping genes, as provided by biological a priori knowledge. However, the existence of housekeeping genes common for a very large collection of samples, especially under a wide range of conditions, is questionable.

## Results

We propose to carry out pairwise normalization with respect to multiple references, selected from representative samples. Then the pairwise intermediates are integrated based on a linear model that adjusts the reference effects. Motivated by the notion of housekeeping genes and their statistical counterparts, we adopt the robust least trimmed squares (LTS) regression in pairwise normalization. The effectiveness of the proposed method (MUREN) is compared with other existing tools using some standard data sets. The goodness of normalization emphasizes on preserving possible asymmetric differentiation, whose biological significance is exemplified by a single cell cycle data.

## Conclusion

MUREN performs the RNA-seq normalization using a two-step statistical regression induced from a general principle. We propose that the densities of pairwise differentiations are used to evaluate the goodness of normalization. MUREN adjusts the mode of differentiation toward zero while preserves the skewness due to biological asymmetric differentiation. Moreover, by robustly integrating pre-normalized counts with respect to multiple reference, MUREN is immune to individual outlier samples.

## Installation

`remotes::install_github("hippo-yf/MUREN")`

## Dependencies

- R (>= 4.0.2)
- assertthat (>= 0.2.1)
- doSNOW (>=1.0.18)
- iterators (>=1.0.12)
- magrittr (>=1.5)
- MASS (>=7.3-51.6)
- matrixStats (>=0.56.0)

## Examples

A rat toxicogenomics RNA-seq data
from Munro et al (2014)

```R
library(MUREN)
data(rat_tox_thi)
```

Normalized counts

```R
thi_norm = muren_norm(rat_tox_thi)
```

or, scaling coefficients

```R
(thi_coeff = muren_norm(rat_tox_thi, res_return = 'scaling_coeff'))
```

The following code chunk shows how to combine MUREN and edgeR for differential gene (DE) calling.

The data `rat_tox_thi` has 3 treatment and control samples respectively

```R
library(edgeR)
group = rep(c('treatment', 'control'), each = 3)
thi_de = DGEList(counts = rat_tox_thi[-1], genes = rat_tox_thi[1], group = group)
```

`norm.factors` by TMM

```R
thi_de = calcNormFactors(thi_de)
```

`norm.factos` by MUREN (coerce the median same with that of TMM, **the re-adjustment here is only used for comparison**)

```R
factors_muren = thi_coeff/thi_de$samples$lib.size
factors_muren_adj = factors_muren/median(factors_muren)*median(thi_de$samples$norm.factors)
```

The plot shows the comparison, where the line is identical $y=x$. They may not have much difference when the transcriptional differentiation is symmetric.

<img src=normfactors.png width="80%">

Substitute the `norm.factors`

```R
thi_de$samples$norm.factors = factors_muren
```

Design matrix 

```R
design <- model.matrix(~group)
```

Dispersion

```R
thi_de = estimateDisp(thi_de, design, robust=TRUE)
```

Differential expression

```R
fit <- glmFit(thi_de, design)
lrt <- glmLRT(fit)
topTags(lrt)
```

## Citation

Feng, Y., Li, L.M. MUREN: a robust and multi-reference approach of RNA-seq transcript normalization. BMC Bioinformatics 22, 386 (2021). <https://doi.org/10.1186/s12859-021-04288-0>

## References

Munro SA, Lund SP, Pine PS, et al. Assessing technical performance in differential gene expression experiments with external spike-in RNA control ratio mixtures. Nat Commun. 2014;5:5125. doi:10.1038/ncomms6125

Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010;26(1):139-140.

McCarthy DJ, Chen Y, Smyth GK. Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Res. 2012;40(10):4288-4297.

## Licence

GPL-3


## Change Logs

### version 0.2.2 (2023-02-27)

- fix bug of double parameter model
- fix bug of input data type
- optimize memory usage

### version 0.1.0 (2021)

- first release

