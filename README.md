# Scope

The R-package `tigaR` performs temporal integrative differential expression analysis of omics data, both sequencing counts RNAseq and continuous microarray data.

# Installation

The package `tigaR` depends on [ShrinkBayes](https://github.com/markvdwiel/ShrinkBayes) and on [R >= 3.0.0](https://cran.r-project.org/) and is available from GitHub. This requires the package [devtools](https://cran.r-project.org/web/packages/devtools/index.html):
```
devtools::install_github("viktormiok/tigaR", build_vignettes = TRUE)
```
Please restart R before loading the package and its documentation:
```
library(tigaR)
utils::help(tigaR)
utils::vignette("tigaR")
```
# References

Miok, V., Wilting, S.M., Van de Wiel, M.A., Jaspers, A., van Noort, P.I., Brakenhoff, R.H., Snijders, P.J.F., Steenbergen, R.D.M., Van Wieringen, W.N. (2014), "tigaR: integrative significance analysis of temporal differential gene expression induced by genomic abnormalities", *BMC Bioinformatics*, 15, 327. ([doi.org/10.1186/1471-2105-15-327](https://doi.org/10.1186/1471-2105-15-327)).

Babion, I., Miok, V., Jaspers, A., Huseinovic, A., Steenbergen, R.D., van Wieringen, W.N., Wilting, S.M. (2020), "Identification of Deregulated Pathways, Key Regulators, and Novel miRNA-mRNA Interactions in HPV-Mediated Transformation", *Cancers*, 12(3), 700. ([doi.org/10.3390/cancers12030700](https://doi.org/10.3390/cancers12030700)).