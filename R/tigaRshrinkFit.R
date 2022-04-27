###############################################################################
#
#         fit of the tigaR models using INLA and shrinkage priors
#
###############################################################################

tigaRshrinkFit <- function (forms, 
                            dat,
                            dat1=NULL,
                            timefac,
                            groupfac,
                            ZSpline,
                            shrinksimul,
                            dispersefixed=10, 
                            disperseaddfixed=1, 
                            disperserandom=1,
                            maxprecfixed=4, 
                            fams="gaussian",
                            ncpus=2,
                            effoutput=TRUE, 
                            keepmargrand=FALSE,
                            keepmarghyper=TRUE,
                            setthreads1=TRUE, 
                            showupdate=FALSE, 
                            silentINLA=TRUE, 
                            updateby=5000,
                            ndigits=3, 
                            addpackage=c("splines"),
                            safemode=TRUE,
                            orthogonal=FALSE, 
                            shrink=FALSE,
                            multivar=FALSE,
                            rho=0, ...){
            if (!is(dat, "matrix")) {
              stop("Input (dat) is of wrong class.")
            }
            if (!is(ZSpline, "matrix")) {
              stop("Input (ZSpline) is of wrong class.")
            }
            if (!is(forms, "formula")) {
              stop("Input (forms) is of wrong class.")
            }
            if (!is(shrinksimul, "list")) {
              stop("Input (shrinksimul) is of wrong class.")
            }
            if (!is(groupfac, "factor")) {
              stop("Input (groupfac) is of wrong class.")
            }
            if ((nrow(ZSpline) != length(groupfac)) & 
                (nrow(ZSpline) != length(timefac))) {
              stop("Dimension of the design matrices do not correspond")
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
            typelik <- shrinksimul$typelik
            ip <- shrinksimul$inputpar
            curvedispfun <- shrinksimul$curvedispfun
            if (typelik == "count") {
              logdisp_shr <- c(mu=shrinksimul$pmlist$mudisp, 
                               prec=shrinksimul$pmlist$precdisp
              )
              logitp0_shr <- c(mu=shrinksimul$pmlist$mup0,
                               prec=shrinksimul$pmlist$precp0
              )
            }
            else {
              logdisp_shr = c(0, 0.01)
              logitp0_shr = c(0, 0.01)
            }
            if (typelik == "gaussian") {
              precerr_shr <- c(shapeerr=shrinksimul$pmlist$shapeerr, 
                               rateerr=shrinksimul$pmlist$rateerr
              )
            }
            else {
              precerr_shr <- c(0.001, 0.001)
            }
            randeff_shr <- c(shape=shrinksimul$pmlist$shaperand,
                             rate=shrinksimul$pmlist$raterand
            )
            randeff_shr <- c(shape=randeff_shr[1]/disperserandom, 
                             rate=randeff_shr[2]/disperserandom
            )
            shrinkrandom <- ip$shrinkrandom
            if (!is.null(shrinkrandom)) {
              if (mode(forms) == "call")
                form_shr <- randreplace(forms,
                                        shrinkrandom, 
                                        randeff_shr
                )
              else form_shr <- lapply(forms, 
                                      randreplace, 
                                      shrinkrandom=shrinkrandom, 
                                      initrandomprec=randeff_shr
              )
            }
            if (is.null(shrinkrandom)) {
              form_shr <- forms
            }
            precfixeddisp <- shrinksimul$pmlist$precfixed/dispersefixed
            precfixeddisp <- min(maxprecfixed, 
                                 precfixeddisp
            )
            mufixed <- shrinksimul$pmlist$mufixed
            precaddfixeddisp <- shrinksimul$pmlist$precaddfixed/disperseaddfixed
            precaddfixeddisp <- min(maxprecfixed,
                                    precaddfixeddisp
            )
            muaddfixed <- shrinksimul$pmlist$muaddfixed
            shrinkfixed <- ip$shrinkfixed
            shrinkaddfixed <- ip$shrinkaddfixed
            addfixed <- shrinksimul$addfixed
            addfixed <- lapply(addfixed,
                               function(af) c(af[1],
                                              min(maxprecfixed,
                                                  af[2]/disperseaddfixed)
                               )
            )
            cf <- list(mean=list(default=0), 
                       prec=list(default=0.01)
            )
            if (!is.null(shrinkfixed)) {
              if (is.factor(try(get(shrinkfixed), silent=TRUE)))
                shrinkfixed <- fact2vec(shrinkfixed)
              nm <- names(cf[[1]])
              for (j1 in 1:length(shrinkfixed)) {
                cf$mean <- c(cf$mean, list(mufixed))
                cf$prec <- c(cf$prec, list(precfixeddisp))
              }
              names(cf[[1]]) <- c(nm, shrinkfixed)
              names(cf[[2]]) <- c(nm, shrinkfixed)
            }
            if (!is.null(shrinkaddfixed)) {
              for (i in 1:length(shrinkaddfixed)) {
                shrinkaddfixedi <- shrinkaddfixed[i]
                addfixedi <- as.numeric(addfixed[[i]])
                nm <- names(cf[[1]])
                if (is.factor(try(get(shrinkaddfixedi), silent=TRUE)))
                  shrinkaddfixedi <- fact2vec(shrinkaddfixedi)
                for (j1 in 1:length(shrinkaddfixedi)) {
                  cf$mean <- c(cf$mean, list(addfixedi[1]))
                  cf$prec <- c(cf$prec, list(addfixedi[2]))
                }
                names(cf[[1]]) <- c(nm, shrinkaddfixedi)
                names(cf[[2]]) <- c(nm, shrinkaddfixedi)
              }
            }
            if ((fams[1] == "poisson" | fams[1] == "zip") & 
                (shrinksimul$pmlist$mixp[1] <0.1)) {
              cat("Estimated proportion to follow (ZI-)
                  Poisson distribution is very low.")
              cat("You may want to abort this computation 
                  and continu with the results from the (ZI-)NB fit")
            }
            
            toret <- fitINLAtigaR(forms=form_shr, 
                                  dat=dat, 
                                  dat1=dat1,
                                  timefac=timefac, 
                                  groupfac=groupfac,
                                  ZSpline=ZSpline,
                                  shrinksimul=shrinksimul,
                                  fams=fams, 
                                  logdisp=logdisp_shr,
                                  precerr=precerr_shr, 
                                  curvedispfun=curvedispfun, 
                                  logitp0=logitp0_shr,
                                  ncpus=ncpus, 
                                  effoutput=effoutput,
                                  keepmargrand=keepmargrand,
                                  keepmarghyper=keepmarghyper, 
                                  setthreads1=setthreads1,
                                  showupdate=showupdate, 
                                  silentINLA=silentINLA, 
                                  updateby=updateby,
                                  ndigits=ndigits, 
                                  addpackage=c("splines"),
                                  cf=cf,
                                  safemode=safemode,
                                  control.predictor=list(compute=TRUE),
                                  orthogonal=orthogonal,
                                  shrink=shrink, 
                                  multivar=multivar, 
                                  rho=rho, ...)
            priors <- list(mufixed=mufixed,
                           precfixed=precfixeddisp,
                           muaddfixed=muaddfixed,
                           precaddfixed=precaddfixeddisp,
                           shaperand=randeff_shr[1], 
                           raterand=randeff_shr[2],
                           mudisp=logdisp_shr[1],
                           precdisp=logdisp_shr[2],
                           mup0=logitp0_shr[1],
                           precp0=logitp0_shr[2],
                           shapeerr=precerr_shr[1],
                           rateerr=precerr_shr[2]
            )
            return(list(res=toret, 
                        priors=priors)
            )
}
