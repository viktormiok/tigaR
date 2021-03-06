# Scope

The R-package `tigaR` performs temporal integrative differential expression analysis of omics data, both sequencing counts RNAseq and continuous microarray data.

The temporal variation in gene expression is described by a generalised linear mixed model employing low-rank thin-plate splines. Model parameters are estimated with an empirical Bayes procedure, which exploits integrated nested Laplace approximation for fast computation. Iteratively, posteriors of hyperparameters and model parameters are estimated. The empirical Bayes procedure shrinks multiple dispersion-related parameters. Shrinkage leads to more stable estimates of the model parameters, better control of false positives and improvement of reproducibility. In addition, to make estimates of the DNA copy number more stable, model parameters are also estimated in a multivariate way using triplets of features, imposing a spatial prior for the copy number effect.

With the proposed R-package for analysis of time-course multilevel molecular continuous (microarray) and count (RNAseq) data, more profound insight may be gained through:
 - identification of temporal differentially gene expression, where method yields improvements in sensitivity, specificity and reproducibility compared to existing methods.
 - identification of temporal differential expression induced by DNA copy number abnormalities or/and miRNA expression levels.
 - identification of miRNA targets to the mRNA gene expression.

Note: if you have a choice to use either Windows or Unix/Linux, opt for the latter. `tigaR` runs more efficiently under Unix/Linux than under Windows. NOTE:  when running `tigaR` you may see *** WARNINGS ***  from `INLA` (e.g. on eigenvalues, or on convergence, or even something like 18500 Aborted...). They can currently not be surpressed, because they are produced by C-code. Please ignore them. 

# Installation

The package `tigaR` depends on [ShrinkBayes](https://github.com/markvdwiel/ShrinkBayes) and on [R >= 3.0.0](https://cran.r-project.org/) and is available from GitHub. This requires the package [devtools](https://cran.r-project.org/web/packages/devtools/index.html):

``` r
devtools::install_github("viktormiok/tigaR", build_vignettes = TRUE)
```

Please restart R before loading the package and its documentation:

``` r
library(tigaR)
utils::help(tigaR)
utils::vignette("tigaR")
```

# References

Publications related to `tigaR` include:

- Miok, V., Wilting, S.M., Van de Wiel, M.A., Jaspers, A., van Noort, P.I., Brakenhoff, R.H., Snijders, P.J.F., Steenbergen, R.D.M., Van Wieringen, W.N. (2014), "tigaR: integrative significance analysis of temporal differential gene expression induced by genomic abnormalities", *BMC Bioinformatics*, 15, 327. ([doi.org/10.1186/1471-2105-15-327](https://doi.org/10.1186/1471-2105-15-327)).

- Babion, I., Miok, V., Jaspers, A., Huseinovic, A., Steenbergen, R.D., van Wieringen, W.N., Wilting, S.M. (2020), "Identification of Deregulated Pathways, Key Regulators, and Novel miRNA-mRNA Interactions in HPV-Mediated Transformation", *Cancers*, 12(3), 700. ([doi.org/10.3390/cancers12030700](https://doi.org/10.3390/cancers12030700)).

Please cite the relevant publications if you use `tigaR`.

