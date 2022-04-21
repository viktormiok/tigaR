###############################################################################
#
#         calculate degrees of freedom used in the testing
#
###############################################################################

degFreedom <- function(fit, 
                       design,
                       ZSpline){
  # fit: Fitted object of the full model
  # design: Design matrix for replicates
  # ZSpline:  Design matrix for splines
  if (!is(fit, "list")) {
    stop("Input (fit) is of wrong class.")
  }
  if (!is(design, "matrix")) {
    stop("Input (design) is of wrong class.")
  }
  if (!is(ZSpline, "matrix")) {
    stop("Input (ZSpline) is of wrong class.")
  }
  if (nrow(ZSpline) != nrow(design)) {
    stop("Dimension of the design matrices do not correspond")
  }
  dMat <- cbind(design, ZSpline)
  I <- diag(ncol(dMat))
  m_e <- inla.emarginal(function(x) 1/sqrt(x),
                        fit$marginals.hyper[[2]]
  )
  m_b <- inla.emarginal(function(x) sqrt(x),
                        fit$marginals.hyper[[1]]
  )
  lambda = m_e*m_b
  if(is.null(lambda)){
    print("Problem problem with estimation")}
  else{        
    df <- sum(diag(dMat %*% solve(t(dMat) %*% dMat + 
                                    lambda*I) %*% t(dMat))) -
      (ncol(dMat)-ncol(ZSpline))
    return(df)
  }
}
