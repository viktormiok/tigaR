###############################################################################
#
#         fit of the tigaR models using INLA
#
###############################################################################

fitINLAtigaR <- function (forms, 
                          dat,
                          dat1 = NULL, 
                          timefac,
                          groupfac, 
                          ZSpline, 
                          fams = "zinb",
                          shrinksimul,
                          logdisp = c(0, 0.01),
                          precerr = c(1, 10^(-5)),
                          curvedispfun = NULL, 
                          logitp0 = c(0, 0.01),
                          ncpus = 2, 
                          effoutput = TRUE, 
                          keepmargrand = FALSE,
                          keepmarghyper = TRUE, 
                          setthreads1 = TRUE, 
                          showupdate = FALSE, 
                          silentINLA = 2L, 
                          updateby = 5000,
                          ndigits = 5, 
                          addpackage = NULL,
                          safemode = TRUE, 
                          cf = NULL,
                          orthogonal = FALSE, 
                          shrink = FALSE,
                          multivar = FALSE,
                          rho=0, ...){
          if (!is(forms, "formula")) {
            stop("Input (forms) is of wrong class.")
          }
          if (!is(dat, "matrix")) {
            stop("Input (dat) is of wrong class.")
          }
          if (!is(ZSpline, "matrix")) {
            stop("Input (ZSpline) is of wrong class.")
          }
          if (!is(groupfac, "factor")) {
            stop("Input (groupfac) is of wrong class.")
          }
          if (!is(timefac, "integer")) {
            stop("Input (timefac) is of wrong class.")
          }
          if ((nrow(ZSpline) != length(groupfac)) &
              (nrow(ZSpline) != length(timefac))) {
            stop("Dimension of the design matrices do not correspond")
          }
          #if (!is(shrinksimul, "list")) {
          #   stop("Input (shrinksimul) is of wrong class.")
          #}
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
          if (!is(ncpus, "numeric")) {
            stop("Input (ncpus) is of wrong class.")
          }
          if (length(ncpus) != 1) {
            stop("Input (ncpus) is of wrong length.")
          }
          if (is.na(ncpus)) {
            stop("Input (ncpus) is not a positive integer.")
          }
          if (ncpus < 0) {
            stop("Input (ncpus) is not a positive integer.")
          }
          if (!is(effoutput, "logical")) {
            stop("Input (effoutput) is of wrong class.")
          }
          if (!is(keepmargrand, "logical")) {
            stop("Input (keepmargrand) is of wrong class.")
          }
          if (!is(setthreads1, "logical")) {
            stop("Input (setthreads1) is of wrong class.")
          }
          if (!is(showupdate, "logical")) {
            stop("Input (showupdate) is of wrong class.")
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
          
          cat("NOTE: Warnings from INLA (eigenvalues, convergence, abort) 
              can currently not be surpressed. Please ignore (generally)")
          cat("\n")
          if (setthreads1) 
            inla.setOption("num.threads", 1)
          ngene <- nrow(dat)
          if (mode(forms) == "call") 
            forms <- rep(list(forms), ngene)
          if (length(fams == 1)) 
            fams <- rep(fams, ngene)
          if (is.null(cf)) 
            cf <- list(prec.intercept = 0.001)
          else cf <- c(cf, prec.intercept = 0.001)
          form <- forms[[1]]
          frmchr <- as.character(form)[[3]]
          sp <- strsplit(frmchr, "\\+")
          sp <- as.vector(sapply(sp, function(tt) {
            gsub(" ", "", tt)
          }))
          wht <- which(!is.na(match(sp, 
                                    objects(envir = .GlobalEnv)))
          )
          if (length(wht) > 0) {
            sp2 <- sp[wht]
            dfr <- data.frame(get(sp2[1]))
            if (length(sp2) > 1) {
              for (j in 2:length(sp2)) {
                dfr <- cbind(dfr, get(sp2[j]))
              }
            }
            names(dfr) <- sp2
          }
          else dfr <- data.frame()
          if (dim(dfr)[1] > 0) {
            elmiss <- length(which(is.na(dfr)))
            if (elmiss > 0) {
              cat("Design contains missing fixed effects. 
                  Please use function 'ReDefMiss' first.")
              return(NULL)
            }
          }
          
          fitinlaseqi <- function(i, ...) {
            print(i)
            form <- forms[[i]]
            fam = fams[i]
            di <- as.numeric(dat[i, ])
            
            if (dim(dfr)[1] == 0){ 
              dattag <- data.frame(y = di, 
                                   timefac,
                                   groupfac
              )
            }else{dattag <- cbind(data.frame(y = di),
                                  dfr,
                                  timefac,
                                  groupfac
            )
            }
            # add copy number to the data set
            if(!is.null(dat1)){
              d1i <- as.numeric(dat1[i, ])
              if (dim(dfr)[1] == 0){ 
                dattag <- data.frame(y = di, 
                                     x = d1i,
                                     timefac,
                                     groupfac
                )
              }else{dattag <- cbind(data.frame(y = di),
                                    dfr,
                                    x = d1i,
                                    timefac,
                                    groupfac
              )
              }
            }
            if(orthogonal){
              ZSpline = ZSpline -
                (d1i%*%(solve(t(d1i)%*%d1i))%*%(t(d1i)%*%ZSpline))
              if(shrink){
                param = c(shrinksimul$inputpar$randomprec[1],
                          shrinksimul$inputpar$randomprec[2])
              } else {
                param = c(shrinksimul$pmlist$shaperand,
                          shrinksimul$pmlist$raterand)
              }
              form <- y ~ groupfac + x + f(timefac,
                                           model = "z",
                                           Z = ZSpline,
                                           initial = 3,
                                           prior = "loggamma",
                                           param = param
              )
            }
            if (!is.null(curvedispfun)) {
              mulogdisp <- curvedispfun(log(sum(di)))
              logdispi <- logdisp + c(mulogdisp, 0)
            }
            else {
              logdispi <- logdisp
            }
            if (fam == "zinb") {
              faminla = "zeroinflatednbinomial1"
              cd <- list(prior = c("gaussian", "gaussian"), 
                         param = c(logdispi, logitp0)
              )
            }
            if (fam == "zip") {
              faminla = "zeroinflatedpoisson1"
              cd <- list(prior = "gaussian", 
                         param = logitp0
              )
            }
            if (fam == "nb") {
              faminla = "nbinomial"
              cd <- list(prior = "gaussian",
                         param = logdispi
              )
            }
            if (fam == "poisson") {
              faminla = "poisson"
              cd <- list()
            }
            if (fam == "gaussian") {
              faminla = "gaussian"
              cd = list(hyper = list(prec = list(prior = "loggamma", 
                                                 param = precerr)))
            }
            INLA:::inla.dynload.workaround()
            result <- try(inla(formula = form,
                               family = faminla, 
                               data = dattag, 
                               control.family = cd, 
                               silent = silentINLA, 
                               control.fixed = cf, ...))
            if (class(result) == "try-error") 
              result <- NULL
            
            if (is.null(result$mlik)) {
              if ((faminla == "nbinomial" | faminla == "zeroinflatednbinomial1") & 
                  max(di, na.rm = T) > 10^5) {
                maxval <- max(di, na.rm = T)
                if (maxval > 10^7) 
                  di <- round(di/1000)
                if (maxval > 10^6 & maxval <= 10^7) 
                  di <- round(di/100)
                if (maxval > 10^5 & maxval <= 10^6) 
                  di <- round(di/10)
                if (dim(dfr)[1] == 0) 
                  dattag <- data.frame(y = di)
                else dattag <- cbind(data.frame(y = di), dfr)
                INLA:::inla.dynload.workaround()
                result <- try(inla(formula = form, 
                                   family = faminla, 
                                   data = dattag, 
                                   control.family = cd,
                                   silent = silentINLA, 
                                   control.fixed = cf, ...))
              }
            }
            if (class(result) == "try-error") 
              result <- NULL
            else {
              if (effoutput) {
                nm = names(result)
                todel <- c(".control.defaults", "control.inla", 
                           "model.matrix", "lincomb",
                           "control.expert", "control.mode", 
                           "control.results", "control.lincomb", 
                           "control.predictor", "joint.hyper",
                           "misc", "logfile"
                )
                if (!keepmargrand) 
                  todel <- c(todel, "marginals.random")
                if (!keepmarghyper) 
                  todel <- c(todel, 
                             "marginals.hyperpar", 
                             "internal.marginals.hyperpar"
                  )
                wh <- match(todel, nm)
                wh <- wh[!is.na(wh)]
                result <- result[-wh]
              }
              
              if(multivar) {
                #browser()
                # design matrix for replicate
                numgroup <- as.numeric(max(levels(groupfac)))
                numtp <- ceiling(max(timefac)/numgroup)
                nz <- ncol(ZSpline)
                CL <- kronecker(diag(nz), rep(1, numtp))
                # spline design matrix
                Z <- as.matrix(bdiag(rep(list(ZSpline), 3)))
                if(i==1|i==nrow(dat)){
                  # gene expression
                  y <- rep(dat[i,], 3)
                  # fixed effect
                  X <- cbind(CL, as.matrix(bdiag(rep(list(dat1[i,]), 3))))
                }
                if(i!=1&i!=nrow(dat)){
                  y <- c(dat[i-1,],dat[i,],dat[i+1,]) #c(dat[c(i-1, i, i+1),]) 
                  X <- cbind(CL,as.matrix(bdiag(dat1[i-1,],
                                                dat1[i,],
                                                dat1[i+1,])))
                }
                # precision for error
                e <- result$summary.hyper$mean[1]
                # precision for rando effect
                u <- result$summary.hyper$mean[2]
                # if some hyperparameters are not estimated
                if(is.null(e)|is.null(u)){e <- 1
                u <- 1}
                
                sigma <- matrix(0, nz+3, nz+3)
                if(rho!=0){
                  cor <- matrix(c(1, rho, rho^2, rho, 1, rho, rho^2, rho, 1), 3)
                  # covariance matrix for fixed effect
                  sigma[(nz+1):(nz+3), (nz+1):(nz+3)] <- solve(cor)
                }
                V <- solve(Z%*%diag(ncol(Z))%*%t(Z)*(e/u) + diag(nrow(Z)))
                d <- (shrinksimul$pmlist$precfixed/10)/e
                #estiamation of fixed effect
                fix <- solve(t(X) %*% V %*% X + d*sigma) %*% t(X) %*% V %*% y
                # estimation of random effect
                rand <- solve(t(Z) %*% Z +
                                diag(ncol(Z))*(u/e)) %*% t(Z) %*% (y -
                                                                     X %*% fix)
                #
                par <- rbind(fix, rand)
                dmat <- cbind(X, Z)
                rz1 <- nrow(ZSpline) + 1
                rz2 <- 2*nrow(ZSpline)
                result$multiFitSpline <- as.numeric(Z %*% rand +
                                                    CL %*% fix[1:12])[rz1:rz2]
                
                result$multiFit <- as.numeric(dmat %*% par)[rz1:rz2]
              }
              if (is.list(result)) 
                result <- rapply(result, 
                                 signif, 
                                 digits = ndigits, 
                                 how = "replace", 
                                 classes = "matrix"
                )
            }
            return(list(result))
          }
          if (ncpus == 1 | ngene == 1) {
            results <- sapply(1:ngene, fitinlaseqi, ...)
          }
          else {
            sfInit(parallel = TRUE, cpus = ncpus)
            sfLibrary(INLA)
            if (!is.null(addpackage)) {
              for (i in 1:length(addpackage)) {
                api <- addpackage[i]
                sfLibrary(api, character.only = TRUE)
              }
            }
            print("Exporting to slaves")
            mysfExport(forceexport = c("dat"))
            mysfExport(forceexport = c("dat1"))
            print("Exporting done")
            print("Started fitting")
            if (showupdate & updateby < ngene) {
              results <- as.list(rep(NA, ngene))
              sec <- seq(1, ngene, by = updateby)
              secl <- length(sec)
              if (sec[secl] > ngene - 2) 
                sec[secl] <- ngene + 1
              else sec <- c(sec, ngene + 1)
              secl <- length(sec)
              for (k in 1:(secl - 1)) {
                pmt <- proc.time()
                resultk <- sfSapply((sec[k]):(sec[k + 1] - 1), 
                                    fitinlaseqi, ...)
                print(paste(sec[k + 1] - 1, "data rows done"))
                print(proc.time() - pmt)
                results[(sec[k]):(sec[k + 1] - 1)] <- resultk
              }
            }
            else {
              results <- sfSapply(1:ngene, fitinlaseqi, ...)
            }
          }
          if ((fams[1] == "zinb" | fams[1] == "nb") & safemode) {
            whichna <- which(is.na(mliks(results)))
            if (length(whichna) > 0) {
              logdisp <- c(logdisp[1], 100)
              if (length(whichna) == 1 | ncpus == 1) {
                newresults <- sapply(whichna, fitinlaseqi, ...)
              }
              else {
                sfExport("logdisp")
                newresults <- sfSapply(whichna, fitinlaseqi, 
                                       ...)
              }
              results[whichna] <- newresults
            }
          }
          if (ncpus > 1) {
            sfRemoveAll()
            sfStop()
          }
          return(results)
          if (setthreads1) 
            rm("inla.options")
}
