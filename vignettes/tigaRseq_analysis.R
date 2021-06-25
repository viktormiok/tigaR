
rm(list=ls())
library(ShrinkBayes)
library(lattice)

setwd("~/Documents/code/tigaR_1.0.0/")
# load data
load("rnaseq_data.rdata")
datax=dat$datCount

# load code
source("tigaR_functions.R")

# Specify number of cpus to use
numcpus2use <- 4  
# specify the number of genes for fitting
numgen <- 100

###############################################################################
###############################################################################
##           differential gene expression analysis time effect
###############################################################################
###############################################################################

# time points for different splines
timefac <- 1:11
# groups
groupfac <- factor(c(rep(1,6), rep(2,5)))

dsg <- getDesign(timefac=timefac, 
                 groupfac=groupfac,
                 numKnots=3, 
                 deg=2, 
                 diffSpl=TRUE 
)
ZSpline <- dsg$ZSpline
design <- dsg$design

###############################################################################
#
#           shrinkage of parameters of full and reduced model
#
###############################################################################

shrinksimulA <- tigaRshrinkSeq(form=y ~ groupfac + f(timefac,
                                                       model="z",
                                                       Z=ZSpline,
                                                       initial=3,
                                                       prior="loggamma",
                                                       param=c(1, 0.00001)
                              ), 
                              dat=datax,
                              maxiter=3,
                              timefac=timefac,
                              groupfac=groupfac,
                              ZSpline=ZSpline,
                              fams="nb",
                              shrinkfixed="groupfac",
                              shrinkrandom="timefac",
                              ncpus=numcpus2use, 
                              addpackage=c("splines")
)

shrinksimul0 <- tigaRshrinkSeq(form=y ~ groupfac,
                               dat=datax,
                               maxiter=3,
                               timefac=timefac,
                               groupfac=groupfac,
                               ZSpline=ZSpline,
                               fams="nb",
                               shrinkfixed="groupfac",
                               ncpus=numcpus2use, 
                               addpackage=c("splines")
)

###############################################################################
#
#           estimate optimal number of knots 
#
###############################################################################

# optNumKnot <- numKnots(form=y ~ groupfac + f(timefac, 
#                                                model="z",
#                                                Z =ZSpline,
#                                                initial=3,
#                                                prior="loggamma", 
#                                                param=c(1,0.00001)),
#                        dat=datax[1:numgen,],
#                        fams="nb",
#                        timefac=timefac, 
#                        groupfac=groupfac,
#                        deg=2,
#                        shrinksimul=shrinksimulA,
#                        diffSpl=TRUE,
#                        orthogonal=FALSE,
#                        shrink=FALSE,
#                        multivar=FALSE,
#                        rho=0
# )

###############################################################################
#
#           fit the full and reduced model using shrunken parameters
#
###############################################################################
seqFitA <- tigaRshrinkFit(forms=y ~ groupfac + f(timefac, 
                                                   model="z",
                                                   Z =ZSpline,
                                                   initial=3,
                                                   prior="loggamma", 
                                                   param=c(1,0.00001)
                          ),
                          dat=datax[1:numgen,],
                          timefac=timefac, 
                          groupfac=groupfac,
                          ZSpline=ZSpline,
                          fams="nb",
                          shrinksimul=shrinksimulA, 
                          ncpus=numcpus2use
)

seqFit0 <- tigaRshrinkFit(forms=y ~ groupfac,
                           dat=datax[1:numgen,],
                           timefac=timefac, 
                           groupfac=groupfac,
                           ZSpline=ZSpline,
                           fams="nb",
                           shrinksimul=shrinksimul0, 
                           ncpus=numcpus2use
)

###############################################################################
#
#          temporal differential expression test     
#
###############################################################################

Result <- tdge(dat=datax[1:numgen,],
               fitAlt=seqFitA,
               fitNull=seqFit0,
               design=design,
               ZSpline=ZSpline
)
rownames(Result) <- 1:numgen
Result <- Result[order(Result[,3]),]
head(Result)
###############################################################################
#
#           plot fit of the features       
#
###############################################################################
# secify the gene index
i=76
plot_tigaRfit(fit=seqFitA[[1]][[i]],
              fit1=seqFit0[[1]][[i]], 
              timefac=c(1:6, 1:5),
              groupfac=groupfac,
              lattice=TRUE,
              cycle=FALSE,
              title="Gene" # only for the cycle=TRUE
)



###############################################################################
###############################################################################
##           circadian gene expression identification
###############################################################################
###############################################################################


# time points for cycling genes
timefac <- c(1:6, 1:5)
# groups
groupfac <- factor(c(rep(1,6), rep(2,5)))

dsg <- getDesign(timefac=timefac, 
                 groupfac=groupfac,
                 numKnots=3, 
                 deg=2, 
                 diffSpl=FALSE 
)
ZSpline <- dsg$ZSpline
design <- dsg$design

###############################################################################
#
#           shrinkage of parameters of full and reduced model
#
###############################################################################

shrinksimulA <- tigaRshrinkSeq(form=y ~ f(timefac,
                                            model="z",
                                            Z=ZSpline,
                                            initial=3,
                                            prior="loggamma",
                                            param=c(1, 0.00001)
                              ), 
                              dat=datax,
                              maxiter=3,
                              timefac=timefac,
                              groupfac=groupfac,
                              ZSpline=ZSpline,
                              fams="nb",
                              #shrinkfixed="groupfac",
                              shrinkrandom="timefac",
                              ncpus=numcpus2use, 
                              addpackage=c("splines")
)

shrinksimul0 <- tigaRshrinkSeq(form=y ~ 1,
                               dat=datax,
                               maxiter=3,
                               timefac=timefac,
                               groupfac=groupfac,
                               ZSpline=ZSpline,
                               fams="nb",
                               #shrinkfixed="groupfac",
                               ncpus=numcpus2use, 
                               addpackage=c("splines")
)

###############################################################################
#
#           estimate optimal number of knots 
#
###############################################################################

# optNumKnot <- numKnots(form=y ~ f(timefac, 
#                                     model="z",
#                                     Z =ZSpline,
#                                     initial=3,
#                                     prior="loggamma", 
#                                     param=c(1,0.00001)),
#                        dat=datax[1:numgen,],
#                        fams="nb",
#                        timefac=timefac, 
#                        groupfac=groupfac,
#                        deg=2,
#                        shrinksimul=shrinksimulA,
#                        diffSpl=TRUE,
#                        orthogonal=FALSE,
#                        shrink=FALSE,
#                        multivar=FALSE,
#                        rho=0
# )

###############################################################################
#
#           fit the full and reduced model using shrunken parameters
#
###############################################################################
seqFitA <- tigaRshrinkFit(forms=y ~ f(timefac, 
                                        model="z",
                                        Z =ZSpline,
                                        initial=3,
                                        prior="loggamma", 
                                        param=c(1,0.00001)
                          ),
                          dat=datax[1:numgen,],
                          timefac=timefac, 
                          groupfac=groupfac,
                          ZSpline=ZSpline,
                          fams="nb",
                          shrinksimul=shrinksimulA, 
                          ncpus=numcpus2use
)

seqFit0 <- tigaRshrinkFit(forms=y ~ 1,
                          dat=datax[1:numgen,],
                          timefac=timefac, 
                          groupfac=groupfac,
                          ZSpline=ZSpline,
                          fams="nb",
                          shrinksimul=shrinksimul0, 
                          ncpus=numcpus2use
)

###############################################################################
#
#          temporal differential expression test     
#
###############################################################################

Result <- tdge(dat=datax[1:numgen,],
               fitAlt=seqFitA,
               fitNull=seqFit0,
               design=design,
               ZSpline=ZSpline
)
rownames(Result) <- 1:numgen
Result <- Result[order(Result[,3]),]
head(Result)
###############################################################################
#
#           plot fit of the features       
#
###############################################################################
# secify the gene index
i=89
plot_tigaRfit(fit=seqFitA[[1]][[i]],
              fit1=seqFit0[[1]][[i]], 
              timefac=c(1:6, 1:5),
              groupfac=groupfac,
              lattice=FALSE,
              cycle=TRUE,
              title="Gene" # only for the cycle=TRUE
)

