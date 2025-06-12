<img src="https://github.com/viktormiok/viktormiok.me/blob/main/software/tigar.png" align="right" height="200" width="200">

![](https://img.shields.io/badge/language-R-orange.svg) ![version](https://img.shields.io/badge/GiHub_version-1.1.0-519dd9) ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/viktormiok/tigaR) ![GitHub issues](https://img.shields.io/github/issues/viktormiok/tigaR)

![dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-orange)  	![commit](https://img.shields.io/github/last-commit/viktormiok/tigaR) ![GitHub](https://img.shields.io/github/license/viktormiok/tigaR)

[![Edit with Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/viktormiok/tigaR) 



# tigaR

- [Overview](#overview)
  * [Application](#application)
- [Installation](#installation)
- [Docker](#docker)
- [Data](#data)
- [Tutorials](#tutorials)
- [License](#license)
- [References](#references)

## Overview

The R-package __`tigaR`__ performs temporal integrative differential expression analysis of omics data, sequencing count RNAseq and continuous microarray data.

A generalized linear mixed model employing low-rank thin-plate splines describes the temporal variation in gene expression. Model parameters are estimated with an empirical Bayes procedure, which exploits integrated nested Laplace approximation for fast computation. Iteratively, the posteriors of hyperparameters and model parameters are estimated. The empirical Bayes procedure shrinks multiple dispersion-related parameters. Shrinkage leads to more stable estimates of the model parameters, better control of false positives, and improvement of reproducibility. In addition, to make the estimates of the DNA copy number more stable, model parameters are multivariate using triplets of features, imposing a spatial prior for the copy number effect.

### Application

With the proposed R-package __`tigaR`__ for analysis of time-course multilevel molecular continuous (microarray) and count (RNAseq) data, more profound insight may be gained through:
 - Identification of temporal differential gene expression, where the method yields improvements in sensitivity, specificity, and reproducibility compared to existing methods.
 Employing the same spline in modeling up/down-regulated genes may be identified across cell lines, while using different splines may allow more flexibility in capturing the temporal variation over time.
 - Identify temporal differential expression induced by DNA copy number abnormalities and/or miRNA expression levels.
 - Identify miRNA targets for the mRNA gene expression in a time-course fashion.
 
<img src="https://user-images.githubusercontent.com/22052679/148564343-38e60761-cb5e-4e1d-966a-77e541a7d1e1.png" align="top" height="540" width="600">


**Note:** If you can choose to use Windows or Unix/Linux, opt for the latter. __`tigaR`__ runs more efficiently under Unix/Linux than Windows. NOTE:  When running __`tigaR`__ you may see *** WARNINGS ***  from [__`INLA`__](https://www.r-inla.org/) (e.g. on eigenvalues, or convergence, or even something like 18500 Aborted...). They can currently not be suppressed, because C-code produces them. Please ignore them. 

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

## Docker

If your system configuration makes installing __`tigaR`__ natively difficult, a Docker container is an alternative way to get __`tigaR`__ running.

**Note:** Docker Machine has Memory and CPU limits on Mac OS X. To control it, please check instructions either for [CLI](https://stackoverflow.com/questions/32834082/how-to-increase-docker-machine-memory-mac/32834453#32834453) or for [Docker Desktop](https://docs.docker.com/docker-for-mac/#advanced).

For building a Docker image from the Dockerfile, download the [Dockerfile](https://github.com/viktormiok/tigaR/blob/main/Dockerf) (available in this repo) and run the following command to create it:
```
docker build -t tigaR .
```
This will create a __`tigaR`__ docker image on your system (please be patient, as the build could take approximately 30-50 minutes to finish).
You can then run it using the following command:
```
docker run -d -p 8787:8787 -e PASSWORD=pass --name tigaR -it tigaR
```

## Data
All the data required for performing temporal integrative genomics analysis and published in the reference articles have been deposited in the National Center for Biotechnology Information Gene Expression Omnibus (GEO). They are accessible through the GEO Series accession numbers:
| Data type     | GEO number |
| ------------- | ------------- |
| CGH Arrays  | [__`GSE138724`__](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4117045)  |
| mRNA Arrays  | [__`GSE138079`__](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138079)  |
| miRNA Arrays  | [__`GSE78279`__](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78279)  |

To access one of the data sets for instance GSE78279 you need to run the code below. Unpacking the data requires tar and gunzip, which should already be available on most systems.

```
cd ../  #To get to the main GitHub repo folder
mkdir -p data/tigaR_data_analysis/
cd data/tigaR_data_analysis/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78279/suppl/GSE78279_RAW.tar
mkdir GSE78279_RAW
tar -C GSE78279_RAW -xvf GSE78279_RAW.tar
gunzip GSE78279_RAW/*_Regional_*
```

## Tutorials

Please see the following tutorials for detailed examples of how to use __`tigaR`__: 

### tigaR walkthrough:
* [PDF version](https://github.com/viktormiok/tigaR/blob/main/vignettes/tigaR%20vignette/tigaRvignette.pdf)
* [R version - array data](https://github.com/viktormiok/tigaR/blob/main/vignettes/tigaRarray_analysis.R)
* [R version - RNA-seq data](https://github.com/viktormiok/tigaR/blob/main/vignettes/tigaRseq_analysis.R)


## License

__`tigaR`__ is distributed under the MIT license. Please read the license before using __`tigaR`__, distributed in the `LICENSE` file.

## References

Publications related to __`tigaR`__ include:

- **Miok, V.**, Wilting, S.M., Van de Wiel, M.A., Jaspers, A., van Noort, P.I., Brakenhoff, R.H., Snijders, P.J.F., Steenbergen, R.D.M., Van Wieringen, W.N. (2014), "[tigaR: integrative significance analysis of temporal differential gene expression induced by genomic abnormalities](https://doi.org/10.1186/1471-2105-15-327)", *BMC Bioinformatics*, 15, 327.

- Babion, I., **Miok, V.**, Jaspers, A., Huseinovic, A., Steenbergen, R.D., van Wieringen, W.N., Wilting, S.M. (2020), "[Identification of Deregulated Pathways, Key Regulators, and Novel miRNA-mRNA Interactions in HPV-Mediated Transformation](https://doi.org/10.3390/cancers12030700)", *Cancers*, 12(3), 700.

- Babion, I., **Miok, V.**, Jaspers, A., van Wieringen, W.N., Steenbergen, R.D., Wilting, S.M. (2018), "[Comprehensive molecular profiling of HPV-induced transformation over time](https://cancerres.aacrjournals.org/content/78/13_Supplement/5059.short)", *Cancer Research*, 22(3), 200. 


Please cite the relevant publications if you use __`tigaR`__.

