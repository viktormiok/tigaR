rm(list=ls())
library(ShrinkBayes)
library(lattice)

setwd("~/Documents/code/tigaR_1.0.0/")
# load data
load("GEhpv3.RData")
head(GEhpv3)
load("CNhpv3.RData")
head(CNhpv3)

# load the code
source("tigaR_functions.R")

# Specify number of cpus to use
numcpus2use <- 4  
# specify the number of genes for fitting
numgen <- 100

###############################################################################
###############################################################################
##           differential gene expression analysis - copy number effect
###############################################################################
###############################################################################

# time points for different splines
timefac <- 1:32
# groups
groupfac <- factor(rep(1:4, each=8))

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
#           Shrinkage of parameters of full and reduced model
#
###############################################################################

shrinksimulA <- tigaRshrinkGauss(form=y ~ groupfac + x + f(timefac,
                                                             model="z",
                                                             Z=ZSpline,
                                                             initial=3,
                                                             prior="loggamma",
                                                             param=c(1, 10^(-5))
                                  ), 
                                  dat=GEhpv3,
                                  dat1=CNhpv3,
                                  maxiter=3,
                                  timefac=timefac,
                                  groupfac=groupfac,
                                  ZSpline=ZSpline,
                                  shrinkfixed="groupfac",
                                  shrinkaddfixed="x",
                                  shrinkrandom="timefac",
                                  ncpus=numcpus2use, 
                                  addpackage=c("splines"),
                                  orthogonal=TRUE,
                                  shrink=FALSE,
)

shrinksimul0 <- tigaRshrinkGauss(form=y ~ groupfac + x,
                                 dat=GEhpv3,
                                 dat1=CNhpv3,
                                 maxiter=3,
                                 timefac=timefac,
                                 groupfac=groupfac,
                                 ZSpline=ZSpline,
                                 shrinkfixed="groupfac",
                                 shrinkaddfixed="x",
                                 ncpus=numcpus2use, 
                                 addpackage=c("splines"),
                                 orthogonal=TRUE,
                                 shrink=FALSE
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
#                        dat=GEhpv3[1:20,],
#                        fams="gaussian",
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

seqFitA <- tigaRshrinkFit(forms=y ~ groupfac + x + f(timefac, 
                                                       model="z",
                                                       Z =ZSpline,
                                                       initial=3,
                                                       prior="loggamma", 
                                                       param=c(1,10^(-5))
                          ),
                          dat=GEhpv3[1:numgen,],
                          dat1=CNhpv3[1:numgen,],
                          timefac=timefac, 
                          groupfac=groupfac,
                          ZSpline=ZSpline,
                          fams="gaussian",
                          shrinksimul=shrinksimulA, 
                          ncpus=2,
                          orthogonal=TRUE,
                          shrink=FALSE,
                          multivar=TRUE,
                          rho=0.8
)

seqFit0 <- tigaRshrinkFit(forms=y ~ groupfac + x,
                          dat=GEhpv3[1:numgen,],
                          dat1=CNhpv3[1:numgen,],
                          timefac=timefac, 
                          groupfac=groupfac,
                          ZSpline=ZSpline,
                          fams="gaussian",
                          shrinksimul=shrinksimul0, 
                          ncpus=numcpus2use,
                          orthogonal=TRUE,
                          shrink=FALSE,
                          multivar=TRUE,
                          rho=0.8
)

###############################################################################
#
#          temporal differential expression test     
#
###############################################################################
source("tigaR_functions.R")
Result <- tdge(dat=GEhpv3[1:numgen,],
               dat1=CNhpv3[1:numgen,],
               fitAlt=seqFitA,
               fitNull=seqFit0,
               design=design,
               ZSpline=ZSpline,
               multivar=TRUE
)
rownames(Result) <- 1:numgen
Result <- Result[order(Result[,4]),]
head(Result)
###############################################################################
#
#           plot fit of the features       
#
###############################################################################
# secify the gene index
i=5
plot_tigaRfit(fit=seqFitA[[1]][[i]],
              fit1=seqFit0[[1]][[i]], 
              timefac=rep(rep(1:8), 4),
              groupfac=groupfac,
              lattice=TRUE,
              cycle=FALSE,
              title="Gene", # only for the cycle=TRUE
              multivar=T
)
# for multivariate fit
par(mfrow=c(1,2))
i=14
plot(GEhpv3[i,])
lines(seqFitA[[1]][[i]]$multiFit)
lines(seqFitA[[1]][[i]]$multiFitSpline, col="red")
lines(seqFitA[[1]][[i]]$summary.fitted.values$mean, col="blue")

plot(GEhpv3[i,])
lines(seqFit0[[1]][[i]]$multiFit)
lines(seqFit0[[1]][[i]]$multiFitSpline, col="red")
lines(seqFit0[[1]][[i]]$summary.fitted.values$mean, col="blue")
###############################################################################
###############################################################################
##           differential gene expression analysis - time effect
###############################################################################
###############################################################################

# time points for different splines
timefac <- 1:32
# groups
groupfac <- factor(rep(1:4, each=8))

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
#           Shrinkage of parameters of full and reduced model
#
###############################################################################

shrinksimulA <- tigaRshrinkGauss(form=y ~ groupfac + f(timefac,
                                                         model="z",
                                                         Z=ZSpline,
                                                         initial=3,
                                                         prior="loggamma",
                                                         param=c(1, 10^(-5))
                                 ), 
                                 dat=GEhpv3,
                                 dat1=CNhpv3
                                 maxiter=3,
                                 timefac=timefac,
                                 groupfac=groupfac,
                                 ZSpline=ZSpline,
                                 shrinkfixed="groupfac",
                                 shrinkrandom="timefac",
                                 ncpus=numcpus2use, 
                                 addpackage=c("splines")
)

shrinksimul0 <- tigaRshrinkGauss(form=y ~ groupfac,
                                 dat=GEhpv3,
                                 maxiter=3,
                                 timefac=timefac,
                                 groupfac=groupfac,
                                 ZSpline=ZSpline,
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
#                                                param=c(1,10^(-5))),
#                        dat=GEhpv3[1:20,],
#                        fams="gaussian",
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
                                                   param=c(1,10^(-5))
                          ),
                          dat=GEhpv3[1:numgen,],
                          timefac=timefac, 
                          groupfac=groupfac,
                          ZSpline=ZSpline,
                          fams="gaussian",
                          shrinksimul=shrinksimulA, 
                          ncpus=numcpus2use
)

seqFit0 <- tigaRshrinkFit(forms=y ~ groupfac,
                          dat=GEhpv3[1:numgen,],
                          timefac=timefac, 
                          groupfac=groupfac,
                          ZSpline=ZSpline,
                          fams="gaussian",
                          shrinksimul=shrinksimul0, 
                          ncpus=numcpus2use
)

###############################################################################
#
#          temporal differential expression test     
#
###############################################################################

Result <- tdge(dat=GEhpv3[1:numgen,],
               fitAlt=seqFitA,
               fitNull=seqFit0,
               design=design,
               ZSpline=ZSpline
)
rownames(Result) <- 1:numgen
Result <- Result[order(Result[,4]),]
head(Result)
###############################################################################
#
#           plot fit of the features       
#
###############################################################################
# secify the gene index
i=61
plot_tigaRfit(fit=seqFitA[[1]][[i]],
              fit1=seqFit0[[1]][[i]], 
              timefac=rep(rep(1:8), 4),
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
timefac <- rep(rep(1:8), 4)

# groups
groupfac <- factor(rep(1:4, each=8))

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
#           Shrinkage of parameters of full and reduced model
#
###############################################################################

shrinksimulA <- tigaRshrinkGauss(form=y ~ f(timefac,
                                              model="z",
                                              Z=ZSpline,
                                              initial=3,
                                              prior="loggamma",
                                              param=c(1, 10^(-5))
                                 ), 
                                 dat=GEhpv3,
                                 maxiter=3,
                                 timefac=timefac,
                                 groupfac=groupfac,
                                 ZSpline=ZSpline,
                                 shrinkfixed="groupfac",
                                 shrinkrandom="timefac",
                                 ncpus=numcpus2use, 
                                 addpackage=c("splines")
)

shrinksimul0 <- tigaRshrinkGauss(form=y ~ groupfac,
                                 dat=GEhpv3,
                                 maxiter=3,
                                 timefac=timefac,
                                 groupfac=groupfac,
                                 ZSpline=ZSpline,
                                 shrinkfixed="groupfac",
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
#                                     param=c(1,10^(-5))),
#                        dat=GEhpv3[1:20,],
#                        fams="gaussian",
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
                                        param=c(1,10^(-5))
                          ),
                          dat=GEhpv3[1:numgen,],
                          timefac=timefac, 
                          groupfac=groupfac,
                          ZSpline=ZSpline,
                          fams="gaussian",
                          shrinksimul=shrinksimulA, 
                          ncpus=numcpus2use
)

seqFit0 <- tigaRshrinkFit(forms=y ~ 1,
                          dat=GEhpv3[1:numgen,],
                          timefac=timefac, 
                          groupfac=groupfac,
                          ZSpline=ZSpline,
                          fams="gaussian",
                          shrinksimul=shrinksimul0, 
                          ncpus=numcpus2use
)

###############################################################################
#
#          temporal differential expression test     
#
###############################################################################

Result <- tdge(dat=GEhpv3[1:numgen,],
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
i=50
plot_tigaRfit(fit=seqFitA[[1]][[i]],
              fit1=seqFit0[[1]][[i]], 
              timefac=rep(rep(1:8), 4),
              groupfac=groupfac,
              lattice=FALSE,
              cycle=FALSE,
              title="Gene" # only for the cycle=TRUE
)

