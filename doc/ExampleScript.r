##############################################################################################################################################################################
#
# Filename: TDGEtest.r
# Description: Testing for differential expression.
# Author: Viktorian Miok
# Last modified: 29-11-2012
#
############################################################################################################################################################################################

# clear all.
rm(list=ls())
# set memory size to maximum.
memory.size(4000)

library(tigaR)
data(GEhpv3)
data(CNhpv3R)
setwd("N:/Cervical Cancer Data/tigaR analysis")
load("GeneExpr-RelativeCopyNumber_chr3_matching_overlapPlus.RData")

#CNcervix3 <- copynumber(CNdata)
           
# Specify size of the experiment
numtp <- 8       # Number of time points
numgroup <- 4    # Number of cell lines
# Specify the splines
numKnots <- 2    # Number of knots for splines
deg <- 3         # Degree of the splines      
# Specify number of cpus to use
numcpus2use <- 2    
###################################################################################################################################################################################################
#        Function for creating design matrices
###################################################################################################################################################################################################
# numKnots: Number of knots for splines
# deg: Degree of the splines
# numtp: Number of time points
# numgroup: Number of groups (cell lines)
# DiffSpl: if TRUE every group (cell line) is explained with different spline, if FALSE same spline is used for all groups (cell lines)
# groupfac: cell line effect
# timefac: time effect
# ZSpline: Design matrix for low-rank thin-plate splines.
desMat <- GetDesign(numKnots=numKnots, deg=deg, numtp=numtp, numgroup=numgroup, DiffSpl=TRUE)
                           groupfac <- desMat$groupfac
                           timefac <- desMat$timefac
                           DsgGroup <- desMat$DsgGroup
                           ZSpline <- desMat$ZSpline
groupfacHPV <- factor(rep(1:2, each = 16))
setwd("N:/Cervical Cancer Data/tigaR analysis/Fitting")                           
########################################################################################################################################################################################################

###################################################################################################################################################################################################
#        Function of estimation of model hyper-parameters and parameters
###################################################################################################################################################################################################
# form: Defines the model formula, where groupfac, x and f() repesent cell line, copy number and (random) time effect, repsectively.
#       Combining these effects user have possibility to estimate hyper-parameters of  different models.
# dat: Gene expression data
# dat1: DNA copy number data, if model without copy number is considered effect dat1=NULL
# shrinkfixed: Shrink parameter of copy number using parametric prior which combine point mass at zero and Gaussian distribution, if not shrinkfixed=NULL
# shrinkrandom: Shrink parameters of random time effects using Gamma distribution as parametric prior, if not shrinkrandom=NULL
# shrinksigma: Shrink parameters of the errorusing Gamma distribution as parametric prior, if not shrinksigma=NULL
# mixtrand: Shrink random effect with parametri prior which combine point mass at zero and Gamma distribution, if not mixtrand=NULL
# ncpus: Number of cpus used for parallel computation, if ncpus is larger than actual cpus function will still run.
# addpackage: Additional package which will be shared among cpus.
# orthogonal: Orthogonalize spline design matrix onto copy number data, if not orthogonal=FALSE.
# shrink: if TRUE output will be standard INLA output required for ShrinkHyperPar function, if not than output have object necessar for further analysis.
###################################################################################################################################################################################################
#        Function of model fitting using estimated(shrunken) hyper-paramters
###################################################################################################################################################################################################
# shrinksimul: Object whcih contain shrunken estimates of hyper parameters
# fams: Assumed distribution of the data. If the data are microarray gene expression fams="gaussian".
#       For RNA-seq data fams="nb" and for estimating hyper parameters use ShrinkSeq function from ShrinkBayes package.
# multivar: Fit model multivariate
################################################################################
#            Model described by DNA copy number
################################################################################
form <- y ~ groupfac + x-1
shrinksimul <- ShrinkHyperPar(form=form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                              shrinkfixed="x", shrinkrandom=NULL, shrinksigma=TRUE, mixtrand=FALSE,
                              ncpus=numcpus2use, addpackage=c("splines"), orthogonal=FALSE, shrink=TRUE)

CN <- FitIntAllShrink(forms=form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline, 
                      fams="gaussian", shrinksimul=shrinksimul, ncpus=numcpus2use, orthogonal=FALSE, shrink=FALSE, multivar=FALSE, rho=0)

form <- y ~ groupfac-1
CL <- FitIntAllShrink(forms=form, dat=GEhpv3, dat1=NULL, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline, 
                      fams="gaussian", shrinksimul=shrinksimul, ncpus=numcpus2use, orthogonal=FALSE, shrink=FALSE, multivar=FALSE, rho=0)
save(CN, file="CN.RData")
save(CL, file="CL.RData")
rm(CN,CL)
################################################################################
#            Model described by time (splines)
################################################################################
form <- y ~ groupfac + f(timefac, model="z", Z=ZSpline, initial=3, prior="loggamma", param=c(1,0.00001))-1
shrinksimul <- ShrinkHyperPar(form=form, dat=GEhpv3, dat1=NULL,timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                              shrinkfixed=NULL, shrinkrandom="timefac", shrinksigma=TRUE, mixtrand=FALSE,  # !!!!!!!!!! mixtrand=TRUE only when you use mixture of Gamma and point mass
                              ncpus=numcpus2use, addpackage=c("splines"), orthogonal=FALSE, shrink=TRUE)

diffSpl <- FitIntAllShrink(form, dat=GEhpv3, dat1=NULL, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                           fams="gaussian", shrinksimul=shrinksimul, ncpus=numcpus2use, orthogonal=FALSE, shrink=FALSE, multivar=FALSE, rho=0)
save(diffSpl, file="diffSpl.RData")
rm(diffSpl)
################################################################################
#            Model described by time (splines)
################################################################################
form <- y ~ groupfac + f(timefac, model="z", Z=ZSpline, initial=3, prior="loggamma", param=c(1,0.00001))-1
shrinksimul <- ShrinkHyperPar(form=form, dat=GEhpv3, dat1=NULL,timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                              shrinkfixed=NULL, shrinkrandom="timefac", shrinksigma=TRUE, mixtrand=TRUE,  # !!!!!!!!!! mixtrand=TRUE only when you use mixture of Gamma and point mass
                              ncpus=numcpus2use, addpackage=c("splines"), orthogonal=FALSE, shrink=TRUE)

diffSplMix <- FitIntAllShrink(form, dat=GEhpv3, dat1=NULL, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                           fams="gaussian", shrinksimul=shrinksimul, ncpus=numcpus2use, orthogonal=FALSE, shrink=FALSE, multivar=FALSE, rho=0)
save(diffSplMix, file="diffSplMix.RData")
rm(diffSplMix)
################################################################################
#           Model integrating DNA copy number and time (splines)
################################################################################
form <- y ~ groupfac + x + f(timefac, model="z", Z=ZSpline, initial=3, prior="loggamma", param=c(1,0.00001))-1
shrinksimul <- ShrinkHyperPar(form=form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                              shrinkfixed="x", shrinkrandom="timefac", shrinksigma=TRUE, mixtrand=FALSE,
                              ncpus=numcpus2use, addpackage=c("splines"), orthogonal=FALSE, shrink=TRUE)
                                       
CNdiffSpl <- FitIntAllShrink(form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                             fams="gaussian", shrinksimul=shrinksimul, ncpus=numcpus2use, orthogonal=FALSE, shrink=FALSE, multivar=TRUE, rho=0.5)    
save(CNdiffSpl, file="CNdiffSpl.RData")
rm(CNdiffSpl)
################################################################################
#           Model integrating DNA copy number and time (splines)
################################################################################    , shrinkaddfixed="x:timefac"
form <- y ~ groupfac + x + x:timefac + f(timefac, model="z", Z=ZSpline, initial=3, prior="loggamma", param=c(1,0.00001))-1
shrinksimul <- ShrinkHyperPar(form=form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                              shrinkfixed="x", shrinkaddfixed="x:timefac", shrinkrandom="timefac", shrinksigma=TRUE, mixtrand=FALSE,
                              ncpus=numcpus2use, addpackage=c("splines"), orthogonal=FALSE, shrink=TRUE)
                                       
CNdiffIntSpl <- FitIntAllShrink(form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                             fams="gaussian", shrinksimul=shrinksimul, ncpus=numcpus2use, orthogonal=FALSE, shrink=FALSE, multivar=TRUE, rho=0.5)      
save(CNdiffIntSpl, file="CNdiffIntSpl.RData")
rm(CNdiffIntSpl)
################################################################################
#           Model integrating DNA copy number and time (splines)
################################################################################
form <- y ~ groupfac + x + f(timefac, model="z", Z=ZSpline, initial=3, prior="loggamma", param=c(1,0.00001))-1
shrinksimul <- ShrinkHyperPar(form=form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                              shrinkfixed="x", shrinkrandom="timefac", shrinksigma=TRUE, mixtrand=FALSE,
                              ncpus=numcpus2use, addpackage=c("splines"), orthogonal=TRUE, shrink=TRUE)
                                       
CNdiffSplOrt <- FitIntAllShrink(form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                             fams="gaussian", shrinksimul=shrinksimul, ncpus=numcpus2use, orthogonal=TRUE, shrink=FALSE, multivar=TRUE, rho=0.5)  
save(CNdiffSplOrt, file="CNdiffSplOrt.RData")
rm(CNdiffSplOrt)                                                            
################################################################################
#           Model integrating DNA copy number and time (splines)
################################################################################
form <- y ~ groupfac + x + x:timefac + f(timefac, model="z", Z=ZSpline, initial=3, prior="loggamma", param=c(1,0.00001))-1
shrinksimul <- ShrinkHyperPar(form=form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                              shrinkfixed="x", shrinkaddfixed="x:timefac", shrinkrandom="timefac", shrinksigma=TRUE, mixtrand=FALSE,
                              ncpus=numcpus2use, addpackage=c("splines"), orthogonal=TRUE, shrink=TRUE)
                                       
CNsameIntSplOrt <- FitIntAllShrink(form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                             fams="gaussian", shrinksimul=shrinksimul, ncpus=numcpus2use, orthogonal=TRUE, shrink=FALSE, multivar=TRUE, rho=0.5)                                                
save(CNsameIntSplOrt, file="CNsameIntSplOrt.RData")
rm(CNsameIntSplOrt)
################################################################################
#           HPV effect
################################################################################
form <- y ~ x + f(timefac, model="z", Z=ZSpline, initial=3, prior="loggamma", param=c(1,0.00001))-1
shrinksimul <- ShrinkHyperPar(form=form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                              shrinkfixed="x", shrinkrandom="timefac", shrinksigma=TRUE, mixtrand=FALSE,
                              ncpus=numcpus2use, addpackage=c("splines"), orthogonal=TRUE, shrink=TRUE)
                                       
CNdiffSplOrtHPVno <- FitIntAllShrink(form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                             fams="gaussian", shrinksimul=shrinksimul, ncpus=numcpus2use, orthogonal=TRUE, shrink=FALSE, multivar=TRUE, rho=0.5)    
save(CNdiffSplOrtHPVno, file="CNdiffSplOrtHPVno.RData")
rm(CNdiffSplOrtHPVno)

form <- y ~ groupfacHPV + x + f(timefac, model="z", Z=ZSpline, initial=3, prior="loggamma", param=c(1,0.00001))-1
shrinksimul <- ShrinkHyperPar(form=form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfacHPV, ZSpline=ZSpline,
                              shrinkfixed="x", shrinkrandom="timefac", shrinksigma=TRUE, mixtrand=FALSE,
                              ncpus=numcpus2use, addpackage=c("splines"), orthogonal=TRUE, shrink=TRUE)
                                       
CNdiffSplOrtHPV <- FitIntAllShrink(form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfacHPV, ZSpline=ZSpline,
                             fams="gaussian", shrinksimul=shrinksimul, ncpus=numcpus2use, orthogonal=TRUE, shrink=FALSE, multivar=TRUE, rho=0.5)    
save(CNdiffSplOrtHPV, file="CNdiffSplOrtHPV.RData")
rm(CNdiffSplOrtHPV)
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################


form <- y ~ x + f(timefac, model="z", Z=ZSpline, initial=3, prior="loggamma", param=c(1,0.00001))-1
shrinksimul <- ShrinkHyperPar(form=form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                              shrinkfixed="x", shrinkrandom="timefac", shrinksigma=TRUE, mixtrand=FALSE,
                              ncpus=numcpus2use, addpackage=c("splines"), orthogonal=FALSE, shrink=TRUE)
                                       
CNdiffSplHPVno <- FitIntAllShrink(form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                             fams="gaussian", shrinksimul=shrinksimul, ncpus=numcpus2use, orthogonal=FALSE, shrink=FALSE, multivar=TRUE, rho=0.5)    
save(CNdiffSplHPVno, file="CNdiffSplHPVno.RData")
rm(CNdiffSplHPVno)
form <- y ~ groupfacHPV + x + f(timefac, model="z", Z=ZSpline, initial=3, prior="loggamma", param=c(1,0.00001))-1
shrinksimul <- ShrinkHyperPar(form=form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                              shrinkfixed="x", shrinkrandom="timefac", shrinksigma=TRUE, mixtrand=FALSE,
                              ncpus=numcpus2use, addpackage=c("splines"), orthogonal=FALSE, shrink=TRUE)
                                       
CNdiffSplHPV <- FitIntAllShrink(form, dat=GEhpv3, dat1=CNhpv3R, timefac=timefac, groupfac=groupfac, ZSpline=ZSpline,
                             fams="gaussian", shrinksimul=shrinksimul, ncpus=numcpus2use, orthogonal=FALSE, shrink=FALSE, multivar=TRUE, rho=0.5)    
save(CNdiffSplHPV, file="CNdiffSplHPV.RData")
rm(CNdiffSplHPV)

#########################################################################################################################################################################################################
#          Testing
#########################################################################################################################################################################################################

################################################################################
#          Testing for Temporal Differential Expression - Same Spline - Model with copy number
################################################################################
annotation <- cbind(rownames(CNhpv3R),fData(GEdata)[,5])
Result <- TDGE(dat=GEhpv3, dat1=NULL, FitAlt=sameSpl[[1]], FitNull=CL[[1]], DsgGroup=DsgGroup, DsgSpl=ZSpline, annotation=annotation, multivar=FALSE, sortby=5)
head(Result)
p <- as.numeric(as.character(Result[,5]))
length(p[p<0.05])
################################################################################
#          Testing for DNA Copy Number Effect - Same Spline - Full Model
################################################################################
# Univariate
Result <- TDGEindCN(FitAlt=CNsameSpl[[1]], FitNull=diffSpl[[1]], dat=GEhpv3, annotation=annotation, multivar=FALSE)
head(Result)
p <- as.numeric(as.character(Result[,5]))
length(p[p<0.05])

# Multivariate
Result <- TDGEindCN(FitAlt=CNsameSplOrt[[1]], FitNull=sameSpl[[1]], dat=GEhpv3, annotation=annotation, multivar=TRUE)
head(Result)
p <- as.numeric(as.character(Result[,5]))
length(p[p<0.05])

#########################################################################################################################################################################################################

#########################################################################################################################################################################################################
#          Calculate obtimal number of knots
#########################################################################################################################################################################################################
form <- y ~  groupfac + f(timefac, model="z", Z=ZSpline, initial=3, prior="loggamma", param=c(1,0.00001))-1
NumKnots(form=form,dat=GEcervix3,dat1=CNcervix3,numtp=8,numgroup=4,shrinksimul=shrinksimul,DiffSpl=FALSE)

#########################################################################################################################################################################################################
#          Function for ploting model fit
#########################################################################################################################################################################################################
# fit: Fit object of the model
# fit1: Fit object of the another model in the case that you wnat to compare fit of two models
# grpLabel: Labels required for lattice plot
# Around0: Scale data and fit around zero
# lattice: If TRUE make lattice plot

label <- c("CellLine1","CellLine1","CellLine1","CellLine1","CellLine1","CellLine1","CellLine1","CellLine1",
              "CellLine2","CellLine2","CellLine2","CellLine2","CellLine2","CellLine2","CellLine2","CellLine2",
              "CellLine3","CellLine3","CellLine3","CellLine3","CellLine3","CellLine3","CellLine3","CellLine3",
              "CellLine4","CellLine4","CellLine4","CellLine4","CellLine4","CellLine4","CellLine4","CellLine4")
              
Plot_inlafit(fit=sameSpl[[1]][[i]], fit1=group[[1]][[i]]$multiFit, numtp=numtp, numgroup=numgroup, grpLabel=label, Around0=FALSE, lattice=TRUE)
