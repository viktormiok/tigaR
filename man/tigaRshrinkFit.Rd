% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tigaR.R
\name{tigaRshrinkFit}
\alias{tigaRshrinkFit}
\title{Applies INLA to multiple data rows using shrinkage priors}
\usage{
tigaRshrinkFit(
  forms,
  dat,
  dat1 = NULL,
  timefac,
  groupfac,
  ZSpline,
  shrinksimul,
  dispersefixed = 10,
  disperseaddfixed = 1,
  disperserandom = 1,
  maxprecfixed = 4,
  fams = "gaussian",
  ncpus = 2,
  effoutput = TRUE,
  keepmargrand = FALSE,
  keepmarghyper = TRUE,
  setthreads1 = TRUE,
  showupdate = FALSE,
  silentINLA = TRUE,
  updateby = 5000,
  ndigits = 3,
  addpackage = c("splines"),
  safemode = TRUE,
  orthogonal = FALSE,
  shrink = FALSE,
  multivar = FALSE,
  rho = 0,
  ...
)
}
\arguments{
\item{forms}{Formula, or list of formulas the length of which equals the number of data rows. See inla and f for specification of the model formula.}

\item{dat}{Matrix or data frame containing the data. Rows are features, columns are samples.}

\item{dat1}{Matrix of data frame containing explanatory data set. Rows are features, columns are samples.}

\item{timefac}{Vector representing the time points of the design matrix}

\item{groupfac}{Factor representing the groups of the design matrix}

\item{shrinksimul}{A list object resulting from tigaRshrinkSeq or tigaRshrinkGauss.}

\item{dispersefixed}{Numeric. Inflation factor for the variance of the main fixed parameter.}

\item{disperseaddfixed}{Numeric. Inflation factor for the variance of the additional fixed parameter.}

\item{disperserandom}{Numeric. Inflation factor for the variance of the random effects precision.}

\item{maxprecfixed}{Numeric. Maximum precision used for the main fixed effect.}

\item{fams}{Character string. Either equal to "poisson", "zip" (zero-inflated Poisson), "nb" (negative binomial), or "zinb" (zero-inflated negative binomial), "gaussian": likelihood to be used.}

\item{ncpus}{Integer. The number of cpus to use for parallel computations.}

\item{effoutput}{Boolean. If FALSE, all INLA output will be saved. If TRUE, some fields will be deleted.}

\item{keepmargrand}{Boolean. Do you wish to save the marginals of the random effect regression parameters (beta's)?}

\item{keepmarghyper}{Boolean. Do you wish to save the marginals of the hyper-parameters?}

\item{setthreads1}{Boolean. If TRUE, sequential computation is forced within each call to INLA.}

\item{showupdate}{Boolean. Do you wish to see updates on progression of the computation? TRUE may slow down the computations due to less efficient parallel computation.}

\item{silentINLA}{How to silence the output of inla? NOTE: currently not all INLA output can be silenced.}

\item{updateby}{Integer, only relevant when showupdate=TRUE. Show an update for each update by number of features executed.}

\item{ndigits}{Integer. Numerical precision in digits for the output.}

\item{addpackage}{Character string. Additional package that you wish to export to slave nodes when parallel computing.}

\item{safemode}{Boolean. Only relevant for fams="zinb" or fams="nb". If fitting fails, should it re-run in a safe mode?}

\item{orthogonal}{Boolean. If TRUE, ortogonalized predictor data set (DNA copy number) to the respons data set (mRNA or miRNA gene expression)}

\item{shrink}{Boolean. If TRUE, perform the shrinkage of DNA copy number neighborhoods triplets.}

\item{multivar}{Boolean. If TRUE, perform multivariate estimation of the parameters.}

\item{rho}{Number representing the level of the correlation between the triplets.}

\item{...}{Further arguments passed on to inla.}

\item{Zspline}{Design matrix for thin-plate splines}
}
\value{
This function return number of degree of freedom used in the model testing.
}
\description{
This is a wrapper function that facilitates multiple calls to INLA using parallel processing. It currently supports the following univariate likelihoods: Poisson, negative binomial, zero-inflated negative binomial, Gaussian. It uses the output from ShrinkSeq or ShrinkGauss as input for the prior parameters.
}
