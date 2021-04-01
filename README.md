# MUltiple REference Normalizer (MUREN)

An Accurate RNA-seq Transcript/Gene Quantification Method **Preserving Biological Asymmetric Differentiation**

## Motivation

Normalization of RNA-seq expression data aims at identifying biological expression differentiation between samples by removing unwanted confounding factors. Explicitly or implicitly, the justification of normalization requires a set of housekeeping genes, as provided by biological a priori knowledge. However, the existence of housekeeping genes common for a very large collection of samples, especially under a wide range of conditions, is questionable.

## Results

We propose to carry out pairwise normalization with respect to multiple references, selected from representative samples. Then the pairwise intermediates are integrated based on a linear model that adjusts the reference effects. Motivated by the notion of housekeeping genes and their statistical counterparts, we adopt the robust least trimmed squares (LTS) regression in pairwise normalization. The effectiveness of the proposed method (MUREN) is compared with other existing tools using some standard data sets. The goodness of normalization emphasizes on preserving possible asymmetric differentiation, whose biological significance is exemplified by a single cell cycle data.

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

`data(rat_tox_thi)`

normalized counts

`thi_norm = muren_norm(rat_tox_thi)`

library sizes (scaling coefficients)

`(thi_coeff = muren_norm(rat_tox_thi, res_return = 'library_size'))`

## Citation

Publishing

## Licence

GPL v3
