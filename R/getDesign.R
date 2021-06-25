###############################################################################
#
#         make design matrix for spline and experiment
#
###############################################################################

getDesign <- function(timefac,
                      groupfac,
                      numKnots,
                      deg,
                      diffSpl){
          # timefac: vector of time points
          # groupfac: factor, represent the desing matrix.
          # numKnots: number of knots for spline.
          # deg: degree of the piecewise polynomial, default is 3 
          # diffspl: if is true different spline for each replicate
          
          if (!is(groupfac, "factor")) {
            stop("Input (groupfac) is of wrong class.")
          }
          if (length(timefac) != length(groupfac)) {
            stop("Dimension of the contrasts do not correspond")
          }
          if (!is(numKnots, "numeric")) {
            stop("Input (numKnots) is of wrong class.")
          }
          if (length(numKnots) != 1) {
            stop("Input (numKnots) is of wrong length.")
          }
          if (is.na(numKnots)) {
            stop("Input (numKnots) is not a positive integer.")
          }
          if (numKnots < 0) {
            stop("Input (numKnots) is not a positive integer.")
          }
          if (!is(deg, "numeric")) {
            stop("Input (numKnots) is of wrong class.")
          }
          if (length(deg) != 1) {
            stop("Input (numKnots) is of wrong length.")
          }
          if (is.na(deg)) {
            stop("Input (numKnots) is not a positive integer.")
          }
          if (deg < 0) {
            stop("Input (numKnots) is not a positive integer.")
          }
          if (!is(diffSpl, "logical")) {
            stop("Input (diffSpl) is of wrong class.")
          }
          # design of the experiment
          dd <- data.frame(a=groupfac)
          design <- as.matrix(sparse.model.matrix(~ -1 + a, dd))
          
          # define knots used for the p-spline at equal spaced
          #quantiles of the covariate
          knots <- quantile(unique(timefac),
                            seq(0, 
                                1,
                                length=(numKnots+2))[-c(1, (numKnots+2))]
          )
          # design matrix for random coefficients using a radial basis
          z_k <- (abs(outer(timefac, knots, "-")))^deg
          omega <- (abs(outer(knots, knots, "-")))^deg
          svd.omega <- svd(omega)
          sqrt.omega <- t(svd.omega$v %*% (t(svd.omega$u)*sqrt(svd.omega$d)))
          ZSpline <- t(solve(sqrt.omega, t(z_k)))
          
          # different spline per replicate
          if(diffSpl){
            p <- split(ZSpline, groupfac)
            x <- list()
            for(i in 1:length(levels(groupfac))){
              x[[i]] <- matrix(p[[i]],,numKnots)
            }
            ZSpline <- as.matrix(bdiag(x))
          }
          return(list(design=design,
                      ZSpline=ZSpline
          )
          )
}
