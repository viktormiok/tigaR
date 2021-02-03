###############################################################################
#
#         estimate optimal number of knots for the splines
#
###############################################################################

numKnots <- function(form,
                     dat,
                     dat1 = NULL,
                     fams,
                     timefac, 
                     groupfac,
                     deg,
                     shrinksimul,
                     diffSpl,
                     orthogonal = FALSE,
                     shrink = FALSE,
                     multivar = FALSE,
                     rho = 0){
          if (!is(dat, "matrix")) {
            stop("Input (dat) is of wrong class.")
          }
          if (!is(form, "formula")) {
            stop("Input (forms) is of wrong class.")
          }
          if (!is(shrinksimul, "list")) {
            stop("Input (shrinksimul) is of wrong class.")
          }
          if (!is(groupfac, "factor")) {
            stop("Input (groupfac) is of wrong class.")
          }
          if (length(timefac) != length(groupfac)) {
            stop("Dimension of the contrasts do not correspond")
          }
          if (!is(fams, "character")) {
            stop("Input (fams) is of wrong class.")
          }
          if (!is(fams, "character")) {
            if (!(fams %in% c("nb",
                              "zinb", 
                              "poisson", 
                              "zip",
                              "gaussian"))) {
              stop("Input (fams) ill-specified.")
            }
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
          if (!is(orthogonal, "logical")) {
            stop("Input (orthogonal) is of wrong class.")
          }
          if (!is(shrink, "logical")) {
            stop("Input (shrink) is of wrong class.")
          }
          if (!is(multivar, "logical")) {
            stop("Input (multivar) is of wrong class.")
          }
          if (!is(rho, "numeric")) {
            stop("Input (rho) is of wrong class.")
          }
          if (length(rho) != 1) {
            stop("Input (rho) is of wrong length.")
          }
          dicEst <- output <- numeric()
          
          for(k in 1:nrow(dat)){
            for(i in 2:((max(timefac)/2) + 1)){
              dsg <- getDesign(timefac = timefac, 
                               groupfac = groupfac,
                               numKnots = i, 
                               deg = deg, 
                               diffSpl = diffSpl  
              )
              ZSpline <- dsg$ZSpline
              design <- dsg$design
              fit <- tigaRshrinkFit(forms = form, 
                                    dat = dat[k,,drop = FALSE],
                                    dat1 = dat1[k,,drop = FALSE],
                                    timefac = timefac,
                                    groupfac = groupfac,
                                    ZSpline = ZSpline, 
                                    fams = fams,
                                    shrinksimul = shrinksimul,
                                    ncpus = 1, 
                                    control.compute = list(dic = TRUE,
                                                           mlik = TRUE,
                                                           cpo = FALSE),
                                    orthogonal = orthogonal,
                                    shrink = shrink, 
                                    multivar = multivar,
                                    rho = rho
              )
              dic <- fit[[1]][[1]]$dic$dic
              numPar <- fit[[1]][[1]]$neffp[1]
              df <- suppressWarnings(degFreedom(fit = fit[[1]][[1]], 
                                                design = design,
                                                ZSpline = ZSpline
              )) 
              dicEst[i-1] <- try(dic - numPar + (ncol(design) + 1 + df))
            }
            output[k] <- which.min(dicEst) + 1
            dicEst <- numeric()
          }    
          res <- as.numeric(names(which.max(table(output))))
          cat("Obtimal number of knots = ", res, "\n")
          barplot(table(output),
                  main = "Optimal number of knots for splines",
                  xlab = "Number of knots",
                  ylab = "Abundance"
          )
          return(output)
}
