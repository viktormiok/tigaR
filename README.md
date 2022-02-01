![](https://img.shields.io/badge/language-R-orange.svg) ![version](https://img.shields.io/badge/GiHub_version-1.1.0-519dd9) ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/viktormiok/tigaR) ![GitHub issues](https://img.shields.io/github/issues/viktormiok/tigaR)

![dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-orange)  	![commit](https://img.shields.io/github/last-commit/viktormiok/tigaR) ![GitHub](https://img.shields.io/github/license/viktormiok/tigaR)

[![Edit with Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/viktormiok/tigaR) 



# tigaR

The R-package __`tigaR`__ performs temporal integrative differential expression analysis of omics data, both sequencing counts RNAseq and continuous microarray data.

The temporal variation in gene expression is described by a generalised linear mixed model employing low-rank thin-plate splines. Model parameters are estimated with an empirical Bayes procedure, which exploits integrated nested Laplace approximation for fast computation. Iteratively, posteriors of hyperparameters and model parameters are estimated. The empirical Bayes procedure shrinks multiple dispersion-related parameters. Shrinkage leads to more stable estimates of the model parameters, better control of false positives and improvement of reproducibility. In addition, to make estimates of the DNA copy number more stable, model parameters are also estimated in a multivariate way using triplets of features, imposing a spatial prior for the copy number effect.

With the proposed R-package __`tigaR`__ for analysis of time-course multilevel molecular continuous (microarray) and count (RNAseq) data, more profound insight may be gained through:
 - identification of temporal differentially gene expression, where method yields improvements in sensitivity, specificity and reproducibility compared to existing methods.
 - employing same spline in modeling up/down regulated genes may be identifed over the cell lines, while using different splines more flexiblity may be allowed in capturing the temporal variation over time.
 - identification of temporal differential expression induced by DNA copy number abnormalities or/and miRNA expression levels.
 - identification of miRNA targets to the mRNA gene expression.
 
 
image.png![image](https://user-images.githubusercontent.com/22052679/148564343-38e60761-cb5e-4e1d-966a-77e541a7d1e1.png)


Note: if you have a choice to use either Windows or Unix/Linux, opt for the latter. __`tigaR`__ runs more efficiently under Unix/Linux than under Windows. NOTE:  when running __`tigaR`__ you may see *** WARNINGS ***  from [__`INLA`__](https://www.r-inla.org/) (e.g. on eigenvalues, or on convergence, or even something like 18500 Aborted...). They can currently not be surpressed, because they are produced by C-code. Please ignore them. 

## Installation

The package __`tigaR`__ depends on [__`ShrinkBayes`__](https://github.com/markvdwiel/ShrinkBayes) and on [R >= 3.0.0](https://cran.r-project.org/) and is available from GitHub. This requires the package [__`devtools`__](https://cran.r-project.org/web/packages/devtools/index.html):

``` r
devtools::install_github("viktormiok/tigaR", build_vignettes=TRUE)
```

Please restart R before loading the package and its documentation:

``` r
library(tigaR)
utils::help(tigaR)
utils::vignette("tigaR")
```
## License

__`tigaR`__ is distributed under the MIT license. Please read the license before using __`tigaR`__, which it is distributed in the `LICENSE` file.

## References

Publications related to __`tigaR`__ include:

- Miok, V., Wilting, S.M., Van de Wiel, M.A., Jaspers, A., van Noort, P.I., Brakenhoff, R.H., Snijders, P.J.F., Steenbergen, R.D.M., Van Wieringen, W.N. (2014), "[tigaR: integrative significance analysis of temporal differential gene expression induced by genomic abnormalities](https://doi.org/10.1186/1471-2105-15-327)", *BMC Bioinformatics*, 15, 327.

- Babion, I., Miok, V., Jaspers, A., Huseinovic, A., Steenbergen, R.D., van Wieringen, W.N., Wilting, S.M. (2020), "[Identification of Deregulated Pathways, Key Regulators, and Novel miRNA-mRNA Interactions in HPV-Mediated Transformation](https://doi.org/10.3390/cancers12030700)", *Cancers*, 12(3), 700.

- Babion, I., Miok, V., Jaspers, A., van Wieringen, W.N., Steenbergen, R.D., Wilting, S.M. (2018), "[Comprehensive molecular profiling of HPV-induced transformation over time](https://cancerres.aacrjournals.org/content/78/13_Supplement/5059.short)", *Cancer Research*, 22(3), 200. 


Please cite the relevant publications if you use __`tigaR`__.

