###############################################################################
#
#         simultaneous shrinkage by empirically estimating multiple priors
#                          Gaussian likelihoods
#
###############################################################################

tigaRshrinkGauss <- function(form,
                             dat, 
                             dat1=NULL,
                             timefac, 
                             groupfac,
                             ZSpline,
                             maxiter=10, 
                             shrinkfixed=NULL,
                             shrinkaddfixed=NULL, 
                             shrinkrandom=NULL,
                             shrinkaddrandom=NULL, 
                             shrinksigma=TRUE, 
                             mixtrand=FALSE, 
                             excludefornull=NULL, 
                             fixedmeanzero=FALSE, 
                             addfixedmeanzero=TRUE, 
                             ntag=ifelse(is.null(excludefornull),
                                           c(100, 200, 500, 1000), 
                                           c(1000)), 
                             fixed=c(0, 1/10), 
                             addfixed=c(0, 1/10), 
                             randomprec=c(1, 10^(-5)), 
                             addrandomprec=c(1, 10^(-5)), 
                             precerr=c(1, 10^(-5)), 
                             diracprob0=ifelse((mixtrand | !is.null(excludefornull)),
                                                 0.8, 
                                                 0.2), 
                             fixedseed=TRUE, 
                             ndraw=10000, 
                             safemode=TRUE, 
                             tol=ifelse((mixtrand | !is.null(excludefornull)),
                                          0.005, 
                                          0.01),
                             tolrand=0.02, 
                             mliktol=0.1,
                             orthogonal=FALSE, 
                             shrink=FALSE,
                             multivar=FALSE,
                             rho=0, ...){
              if (!is(form, "formula")) {
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
              if ((nrow(ZSpline) != length(groupfac)) &
                  (nrow(ZSpline) != length(timefac))) {
                stop("Dimension of the design matrices do not correspond")
              }
              if (!is(maxiter, "numeric")) {
                stop("Input (maxiter) is of wrong class.")
              }
              if (length(maxiter) != 1) {
                stop("Input (ncpus) is of wrong length.")
              }
              if (is.na(maxiter)) {
                stop("Input (ncpus) is not a positive integer.")
              }
              if (maxiter < 0) {
                stop("Input (ncpus) is not a positive integer.")
              }
              if (!is(addfixedmeanzero, "logical")) {
                stop("Input (addfixedmeanzero) is of wrong class.")
              }
              if (!is(fixedseed, "logical")) {
                stop("Input (fixedseed) is of wrong class.")
              }
              if (!is(shrinksigma, "logical")) {
                stop("Input (shrinksigma) is of wrong class.")
              }
              if (!is(mixtrand, "logical")) {
                stop("Input (mixtrand) is of wrong class.")
              }
              if (!is(fixedmeanzero, "logical")) {
                stop("Input (fixedmeanzero) is of wrong class.")
              }
              if (!is(safemode, "logical")) {
                stop("Input (safemode) is of wrong class.")
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
              if (mixtrand & !is.null(excludefornull)) {
                cat("Only one mixture prior can be used:
                                    excludefornull is set to NULL")
                exlcludefornull <- NULL
              }
              if (!is.null(excludefornull)) {
                nexcl <- length(excludefornull)
                form0 <- form
                for (j in 1:nexcl){ 
                  form0 <- replacerand0(form0, 
                                        excludefornull[j])
                }  
              }
              if (!is.null(shrinkfixed)) {
                shrinkfixedname <- shrinkfixed
              }
              if (!is.null(shrinkaddfixed)) {
                shrinkaddfixednames <- shrinkaddfixed
              }
              repNA <- function(x) {
                if (is.na(x)) 
                  return(-10^10)
                else return(x)
              }
              iselement2 <- function(elem, set) if (is.null(set) | is.null(elem)) 
                return(FALSE)
              else return(is.element(elem, set))
              els <- length(shrinkaddfixed)
              if (els >= 2) {
                if (mode(addfixedmeanzero) != "list") 
                  addfixedmeanzero <- as.list(rep(addfixedmeanzero, 
                                                  els))
                if (mode(addfixed) != "list") {
                  add0 <- addfixed
                  addfixed <- list()
                  for (i in 1:els) addfixed <- c(addfixed, list(add0))
                }
              }
              else {
                addfixedmeanzero <- list(addfixedmeanzero)
                addfixed <- list(addfixed)
              }
              elsr <- length(shrinkaddrandom)
              if (elsr >= 2) {
                add0 <- addrandomprec
                addrandom <- list()
                for (i in 1:elsr) addrandom <- c(addrandom, list(add0))
              }
              else {
                addrandom <- list(addrandomprec)
              }
              formch <- deparse(form)
              inputpar <- list(formch, shrinkfixed, shrinkaddfixed,
                               shrinkrandom, shrinkaddrandom,
                               shrinksigma, mixtrand, excludefornull, 
                               fixedmeanzero, addfixedmeanzero, 
                               maxiter, tol, tolrand, fixed, 
                               addfixed, randomprec, 
                               addrandomprec, diracprob0, 
                               fixedseed, ndraw, safemode, mliktol
              )
              names(inputpar) <- c("form", "shrinkfixed", 
                                   "shrinkaddfixed", "shrinkrandom",
                                   "shrinkaddrandom", "shrinksigma", 
                                   "mixtrand", "excludefornull", 
                                   "fixedmeanzero", "addfixedmeanzero", 
                                   "maxiter", "tol", "tolrand",
                                   "fixed", "addfixed", "randomprec", 
                                   "addrandomprec", "diracprob0", 
                                   "fixedseed", "ndraw", 
                                   "safemode", "mliktol")
              diracprob <- c(diracprob0, 1 - diracprob0)
              ntagtot <- nrow(dat)
              fams <- rep("gaussian", ntagtot)
              lngene <- length(ntag)
              ntagtotal <- nrow(dat)
              wh2large <- which(ntag > ntagtotal)
              if (length(wh2large) > 0) {
                ntag <- c(ntag[-wh2large], ntagtotal)
              }
              pmlist <- list()
              pmlist <- c(pmlist, list(mufixed=fixed[1], 
                                       precfixed=fixed[2])
              )
              if (!is.null(shrinkaddfixed)) {
                for (i in 1:els) {
                  addf <- as.numeric(addfixed[[i]])
                  name <- c(paste("mu", shrinkaddfixed[i], sep=""), 
                            paste("prec", shrinkaddfixed[i], sep=""))
                  toadd <- list(muaddfixed=addf[1],
                                precaddfixed=addf[2]
                  )
                  names(toadd) <- name
                  pmlist <- c(pmlist, toadd)
                }
              }
              else {
                pmlist <- c(pmlist, 
                            list(muaddfixed=addfixed[[1]][1], 
                                 precaddfixed=addfixed[[1]][2])
                )
              }
              pmlist <- c(pmlist, list(shaperand=randomprec[1],
                                       raterand=randomprec[2], 
                                       mixp=diracprob))
              elsr <- length(shrinkaddrandom)
              if (!is.null(shrinkaddrandom)) {
                for (i in 1:elsr) {
                  addr <- addrandomprec
                  name <- c(paste("shape", shrinkaddrandom[i], sep=""), 
                            paste("rate", shrinkaddrandom[i], sep=""))
                  toadd <- list(shapeaddr=addr[1], rateaddr=addr[2])
                  names(toadd) <- name
                  pmlist <- c(pmlist, toadd)
                }
              }
              else {
                pmlist <- c(pmlist, list(shapeaddr=addrandomprec[1], 
                                         rateaddr=addrandomprec[2]))
              }
              pmlist <- c(pmlist, list(shapeerr=precerr[1], 
                                       rateerr=precerr[2])
              )
              print(pmlist)
              ntagused <- NA
              paraprev <- c(unlist(pmlist), ntagused, meanmlik=NA)
              paraall <- paraprev
              paramidall <- c(NA, NA)
              ksall <- c()
              cf <- list(mean=0, prec=0.01)
              for (j in 1:lngene) {
                ngenej <- ntag[j]
                iter <- 1
                change <- 1
                moreiter <- TRUE
                famlevels <- unique(fams)
                fllen <- length(famlevels)
                if (fixedseed) {
                  set.seed(43545126 + j)
                }
                if (fllen == 1) {
                  sel <- sample(1:ntagtotal, min(ntagtotal, ngenej))
                }
                else {
                  sel <- c()
                  for (j1 in 1:fllen) {
                    wh <- which(fams == famlevels[j1])
                    prop <- length(wh)/ntagtotal
                    nsample <- min(length(wh), round(prop * ngenej))
                    sel <- c(sel, sample(wh, nsample))
                  }
                }
                datshrinkj <- dat[sel, ]
                if (!is.null(dat1)) {
                  CNdatashrinkj <- dat1[sel, ]}
                else {CNdatashrinkj=NULL}
                famsj <- fams[sel]
                mlikmeanall <- c()
                mlikprev <- -Inf
                while (moreiter) {
                  print(paste("iter=", iter))
                  print(paste("ntagsused=", ntag[j]))
                  if (diracprob[1] == 0) {
                    mixtrand <- FALSE
                    excludefornull <- NULL
                    print("Point mass has zero mass...fitting full model only.")
                  }
                  cf <- list(mean=list(default=0),
                             prec=list(default=0.01)
                  )
                  if (!is.null(shrinkfixed)) {
                    if (is.factor(try(get(shrinkfixed), silent=TRUE))) 
                      shrinkfixed <- fact2vec(shrinkfixed)
                    nm <- names(cf[[1]])
                    for (j1 in 1:length(shrinkfixed)) {
                      cf$mean <- c(cf$mean, list(fixed[1]))
                      cf$prec <- c(cf$prec, list(fixed[2]))
                    }
                    names(cf[[1]]) <- c(nm, shrinkfixed)
                    names(cf[[2]]) <- c(nm, shrinkfixed)
                  }
                  if (!is.null(shrinkaddfixed)) {
                    for (i in 1:els) {
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
                  if (!is.null(shrinkrandom)) {
                    form <- randreplace(form, 
                                        shrinkrandom, 
                                        randomprec
                    )
                    formrand0 <- replacerand0(form, shrinkrandom)
                    if (!is.null(excludefornull) & !iselement2(inputpar$shrinkrandom, 
                                                               excludefornull)) {
                      form0 <- randreplace(form0, shrinkrandom, randomprec)
                    }
                  }
                  if (!is.null(shrinkaddrandom)) {
                    for (i in 1:elsr) {
                      form <- randreplace(form, shrinkaddrandom[i], 
                                          addrandom[[i]])
                      if (!is.null(excludefornull) & !iselement2(inputpar$shrinkaddrandom[i], 
                                                                 excludefornull)) {
                        form0 <- randreplace(form0, 
                                             shrinkaddrandom[i], 
                                             addrandom[[i]]
                        )
                      }
                    }
                  }
                  nsamtag <- ceiling(ndraw/ngenej)
                  nsamtagfixed <- ceiling(ndraw/(max(1, 
                                                     length(shrinkfixed)) * 
                                                     ngenej)
                  )
                  repNA <- function(x) {
                    if (is.na(x)) 
                      return(-10^10)
                    else return(x)
                  }
                  fitall <- fitINLAtigaR(forms=form,
                                         dat=datshrinkj, 
                                         dat1=CNdatashrinkj,
                                         timefac=timefac,
                                         groupfac=groupfac,
                                         ZSpline=ZSpline,
                                         fams=famsj, 
                                         precerr=precerr, 
                                         cf=cf, 
                                         control.compute=list(dic=F,
                                                                mlik=T, 
                                                                cpo=F),
                                         orthogonal=FALSE, 
                                         shrink=FALSE,
                                         multivar=FALSE,
                                         rho=0, ...)
                  mliksall <- mliks(fitall)
                  if (mixtrand) {
                    fitallrand0 <- fitINLAtigaR(forms=formrand0, 
                                                dat=datshrinkj, 
                                                dat1=CNdatashrinkj,
                                                timefac=timefac,
                                                groupfac=groupfac,
                                                ZSpline=ZSpline,
                                                precerr=precerr, 
                                                control.compute=list(dic=F, 
                                                                       mlik=T,
                                                                       cpo=F),
                                                cf=cf,
                                                orthogonal=FALSE, 
                                                shrink=FALSE,
                                                multivar=FALSE,
                                                rho=0, ...)
                    mliksall0 <- mliks(fitallrand0)
                    if (!is.null(shrinkfixed)) 
                      postfixed <- fitinlacombine(list(fitallrand0,
                                                       fitall), 
                                                  probs=diracprob,
                                                  modus="fixed", 
                                                  para=shrinkfixed,
                                                  nsam=nsamtagfixed, 
                                                  safemode=safemode
                      )
                    if (!is.null(shrinkaddfixed)) {
                      postaddfixed <- list()
                      for (i in 1:els) {
                        shrinkaddfixedi <- shrinkaddfixed[[i]]
                        postaddfixedi <- fitinlacombine(list(fitallrand0, 
                                                             fitall),
                                                        probs=diracprob,
                                                        modus="fixed", 
                                                        para=shrinkaddfixedi,
                                                        nsam=nsamtagfixed, 
                                                        safemode=safemode
                        )
                        postaddfixed <- c(postaddfixed, list(postaddfixedi))
                      }
                    }
                    if (!is.null(shrinkrandom)) 
                      postrandom <- fitinlacombine(list(fitallrand0,
                                                        fitall),
                                                   probs=diracprob, 
                                                   modus="random", 
                                                   para=shrinkrandom, 
                                                   nsam=nsamtag,
                                                   safemode=safemode
                      )
                    if (!is.null(shrinkaddrandom)) {
                      postaddrandom <- list()
                      for (i in 1:elsr) {
                        shrinkaddrandomi <- shrinkaddrandom[i]
                        postaddrandomi <- fitinlacombine(list(fitallrand0, 
                                                              fitall),
                                                         probs=diracprob,
                                                         modus="random", 
                                                         para=shrinkaddrandomi,
                                                         nsam=nsamtag, 
                                                         safemode=safemode
                        )
                        postaddrandom <- c(postaddrandom, list(postaddrandomi))
                      }
                    }
                    if (shrinksigma) 
                      posterr <- fitinlacombine(list(fitallrand0,
                                                     fitall),
                                                probs=diracprob,
                                                modus="err", 
                                                nsam=nsamtag, 
                                                safemode=safemode
                      )
                  }
                  if (!is.null(excludefornull)) {
                    fitall0 <- fitINLAtigaR(forms=form0, 
                                            dat=datshrinkj, 
                                            dat1=CNdatashrinkj,
                                            timefac=timefac,
                                            groupfac=groupfac,
                                            ZSpline=ZSpline,
                                            fams=famsj, 
                                            precerr=precerr,
                                            cf=cf,
                                            control.compute=list(dic=F,
                                                                   mlik=T, 
                                                                   cpo=F),
                                            orthogonal=FALSE, 
                                            shrink=FALSE,
                                            multivar=FALSE,
                                            rho=0, ...)
                    mliksall0 <- mliks(fitall0)
                    if (!is.null(shrinkfixed)) 
                      postfixed <- fitinlacombine(list(fitall0, 
                                                       fitall), 
                                                  probs=diracprob, 
                                                  modus="fixed", 
                                                  para=shrinkfixed, 
                                                  nsam=nsamtagfixed, 
                                                  safemode=safemode
                      )
                    if (!is.null(shrinkaddfixed)) {
                      postaddfixed <- list()
                      for (i in 1:els) {
                        shrinkaddfixedi <- shrinkaddfixed[[i]]
                        postaddfixedi <- fitinlacombine(list(fitall0, 
                                                             fitall), 
                                                        probs=diracprob,
                                                        modus="fixed", 
                                                        para=shrinkaddfixedi, 
                                                        nsam=nsamtagfixed, 
                                                        safemode=safemode
                        )
                        postaddfixed <- c(postaddfixed, list(postaddfixedi))
                      }
                    }
                    if (!is.null(shrinkrandom)) 
                      postrandom <- fitinlacombine(list(fitall0,
                                                        fitall),
                                                   probs=diracprob, 
                                                   modus="random", 
                                                   para=shrinkrandom,
                                                   nsam=nsamtag,
                                                   safemode=safemode
                      )
                    if (!is.null(shrinkaddrandom)) {
                      postaddrandom <- list()
                      for (i in 1:elsr) {
                        shrinkaddrandomi <- shrinkaddrandom[i]
                        postaddrandomi <- fitinlacombine(list(fitall0, 
                                                              fitall),
                                                         probs=diracprob,
                                                         modus="random", 
                                                         para=shrinkaddrandomi, 
                                                         nsam=nsamtag, 
                                                         safemode=safemode
                        )
                        postaddrandom <- c(postaddrandom, list(postaddrandomi))
                      }
                    }
                    if (shrinksigma) 
                      posterr <- fitinlacombine(list(fitall0, 
                                                     fitall), 
                                                probs=diracprob, 
                                                modus="err", 
                                                nsam=nsamtag, 
                                                safemode=safemode
                      )
                  }
                  if (!mixtrand & is.null(excludefornull)) {
                    if (!is.null(shrinkfixed)) 
                      postfixed <- fitinlacombine(list(fitall), 
                                                  modus="fixed", 
                                                  para=shrinkfixed, 
                                                  nsam=nsamtagfixed
                      )
                    if (!is.null(shrinkaddfixed)) {
                      postaddfixed <- list()
                      for (i in 1:els) {
                        shrinkaddfixedi <- shrinkaddfixed[i]
                        postaddfixedi <- fitinlacombine(list(fitall), 
                                                        modus="fixed",
                                                        para=shrinkaddfixedi, 
                                                        nsam=nsamtagfixed, 
                                                        safemode=safemode
                        )
                        postaddfixed <- c(postaddfixed, list(postaddfixedi))
                      }
                    }
                    if (!is.null(shrinkrandom)) 
                      postrandom <- fitinlacombine(list(fitall), 
                                                   modus="random",
                                                   para=shrinkrandom, 
                                                   nsam=nsamtag
                      )
                    if (!is.null(shrinkaddrandom)) {
                      postaddrandom <- list()
                      for (i in 1:elsr) {
                        shrinkaddrandomi <- shrinkaddrandom[i]
                        postaddrandomi <- fitinlacombine(list(fitall), 
                                                         probs=diracprob,
                                                         modus="random",
                                                         para=shrinkaddrandomi, 
                                                         nsam=nsamtag,
                                                         safemode=safemode
                        )
                        postaddrandom <- c(postaddrandom, list(postaddrandomi))
                      }
                    }
                    if (shrinksigma) 
                      posterr <- fitinlacombine(list(fitall),
                                                modus="err", 
                                                nsam=nsamtag
                      )
                  }
                  if (!mixtrand & is.null(excludefornull)) 
                    mlikmean <- mean(mliksall, na.rm=T)
                  else {
                    if (diracprob[1] == 0) 
                      mlikmean <- mean(mliksall, na.rm=T)
                    if (diracprob[2] == 0) 
                      mlikmean <- mean(mliksall0, na.rm=T)
                    if (diracprob[1] != 0 & diracprob[2] != 0) 
                      mlikmean <- mean(log(diracprob[2]) + 
                                         mliksall + 
                                         log(1 + 
                                             exp((log(diracprob[1]) +
                                                  mliksall0 - 
                                                  log(diracprob[2]) -
                                                  mliksall))), 
                                       na.rm=T
                      )
                  }
                  mlikmeanall <- c(mlikmeanall, mlikmean)
                  print(mlikmeanall)
                  mlikrel <- (mlikmean - mlikprev)/abs(mlikmean)
                  print(mlikrel)
                  mlikconv <- (mlikrel <= (mliktol/100))
                  if (mlikconv) {
                    fixed <- fixedprev
                    addfixed <- addfixedprev
                    addrandom <- addrandomprev
                    randomprec <- randomprecprev
                    precerr <- precerrprev
                    diracprob <- diracprobprev
                  }
                  else {
                    fixedprev <- fixed
                    addfixedprev <- addfixed
                    addrandomprev <- addrandom
                    randomprecprev <- randomprec
                    precerrprev <- precerr
                    diracprobprev <- diracprob
                    if (is.null(excludefornull)) 
                      p0 <- 0
                    else {
                      mlik <- mliks(fitall)
                      mlik0 <- mliks(fitall0)
                      repNA <- function(x) {
                        if (is.na(x)) 
                          return(-10^10)
                        else return(x)
                      }
                      mlik <- sapply(mlik, repNA)
                      mlik0 <- sapply(mlik0, repNA)
                      mlik0 <- mlik + sapply(mlik0 - mlik, 
                                             function(x) min(40, max(x, -40))
                      )
                      maxlik <- as.numeric(apply(cbind(mlik, mlik0), 
                                                 1, max))
                      p0start <- length(which((mlik0 - mlik) > 0))/length(mlik)
                      liktot <- function(p0=0.5) {
                        -sum(log(p0 * exp(mlik0 - maxlik) + 
                                   (1 - p0) * exp(mlik - maxlik)))
                      }
                      res2 <- optimize(liktot, 
                                       lower=0,
                                       upper=1, 
                                       maximum=FALSE)
                      p0 <- res2$minimum
                    }
                    if (!is.null(shrinkfixed)) {
                      if (is.element(inputpar$shrinkfixed, excludefornull)) {
                        wh0 <- which(postfixed == 0)
                        p0 <- length(wh0)/length(postfixed)
                        postfixed <- postfixed[-wh0]
                      }
                      if (length(postfixed) > 1) {
                        if (fixedmeanzero) 
                          fixed <- gaussML(postfixed)
                        else fixed <- gaussMLwithmean(postfixed)
                      }
                    }
                    if (!is.null(shrinkaddfixed)) {
                      addfixed <- list()
                      for (i in 1:els) {
                        addfixedmeanzeroi <- addfixedmeanzero[[i]]
                        postaddfixedi <- unlist(postaddfixed[[i]])
                        if (is.element(inputpar$shrinkaddfixed[i], 
                                       excludefornull)) {
                          wh0 <- which(postaddfixedi == 0)
                          p0 <- length(wh0)/length(postaddfixedi)
                          postaddfixedi <- postaddfixedi[-wh0]
                        }
                        if (length(postaddfixedi) > 1) {
                          if (addfixedmeanzeroi) 
                            addfixedi <- gaussML(postaddfixedi)
                          else addfixedi <- gaussMLwithmean(postaddfixedi)
                        }
                        addfixed <- c(addfixed, list(addfixedi))
                      }
                    }
                    if (!is.null(shrinkrandom)) {
                      if (mixtrand) {
                        mlm <- mixture2loggammaML(postrandom)
                        randomprec <- c(mlm[2], mlm[3])
                        p0 <- mlm[1]
                      }
                      else {
                        if (iselement2(inputpar$shrinkrandom, excludefornull)) {
                          wh0 <- which(postrandom == 0)
                          p0 <- length(wh0)/length(postrandom)
                          postrandom <- postrandom[-wh0]
                        }
                        if (length(postrandom) > 1) {
                          randomprec <- loggammaML(postrandom)
                        }
                      }
                    }
                    if (!is.null(shrinkaddrandom)) {
                      addrandom <- list()
                      for (i in 1:elsr) {
                        postaddrandomi <- unlist(postaddrandom[[i]])
                        if (is.element(inputpar$shrinkaddrandom[i], 
                                       excludefornull)) {
                          wh0 <- which(postaddrandomi == 0)
                          p0 <- length(wh0)/length(postaddrandomi)
                          postaddrandomi <- postaddrandomi[-wh0]
                        }
                        if (length(postaddrandomi) > 1) {
                          addrandomi <- loggammaML(postaddrandomi)
                        }
                        addrandom <- c(addrandom, list(addrandomi))
                      }
                    }
                    diracprob <- c(p0, 1 - p0)
                    
                    if (shrinksigma) {
                      if (length(posterr) > 1) {
                        precerr <- loggammaML(posterr)
                      }
                    }
                    rn <- seq(-8, 8, by=0.01)
                    rnprec <- seq(-8, 1.5, by=0.01)
                    if (!iselement2(inputpar$shrinkfixed, excludefornull)) {
                      KSfixed <- max(abs(pnorm(rn, 
                                               mean=fixedprev[1], 
                                               sd=sqrt(1/fixedprev[2])) - 
                                           pnorm(rn,
                                                 mean=fixed[1], 
                                                 sd=sqrt(1/fixed[2])))
                      )
                    }
                    else {
                      KSfixed <- max(c(abs(diracprobprev[2] * pnorm(rn[rn < 0], 
                                                                    mean=fixedprev[1], 
                                                                    sd=sqrt(1/fixedprev[2])) - 
                                             diracprob[2] * pnorm(rn[rn < 0],
                                                                  mean=fixed[1], 
                                                                  sd=sqrt(1/fixed[2]))), diracprobprev[2] * 
                                         abs(pnorm(rn[rn >= 0],
                                                   mean=fixedprev[1], 
                                                   sd=sqrt(1/fixedprev[2]) + diracprobprev[1]) - 
                                               (diracprob[2] * pnorm(rn[rn >= 0],
                                                                     mean=fixed[1], 
                                                                     sd=sqrt(1/fixed[2])) + +diracprob[1]))))
                    }
                    whi <- which(iselement2(inputpar$shrinkaddfixed, 
                                            excludefornull))
                    KSmaxaddfixed <- 0
                    if (length(whi) > 0) {
                      KSaddfixedall <- sapply(whi, function(i) {
                        addfixedprevi <- addfixedprev[[i]]
                        addfixedi <- addfixed[[i]]
                        KSaddfixedi <- max(c(abs(diracprobprev[2] * pnorm(rn[rn < 0],
                                                                          mean=addfixedprevi[1], 
                                                                          sd=sqrt(1/addfixedprevi[2])) - 
                                                   diracprob[2] * pnorm(rn[rn < 0],
                                                                        mean=addfixedi[1], 
                                                                        sd=sqrt(1/addfixedi[2]))),
                                             diracprobprev[2] * abs(pnorm(rn[rn >= 0],
                                                                          mean=addfixedprevi[1], 
                                                                          sd=sqrt(1/addfixedprevi[2]) + diracprobprev[1]) - 
                                                                      (diracprob[2] * pnorm(rn[rn >= 0], 
                                                                                            mean=addfixedi[1], 
                                                                                            sd=sqrt(1/addfixedi[2])) + +diracprob[1]))))
                      })
                      KSmaxaddfixed <- max(KSaddfixedall)
                    }
                    KSsigma <- max(abs(plgamma(rnprec,
                                               location=0, 
                                               scale=1/precerrprev[2],
                                               precerrprev[1]) - 
                                         plgamma(rnprec, 
                                                 location=0, 
                                                 scale=1/precerr[2], 
                                                 precerr[1])))
                    if (!iselement2(inputpar$shrinkrandom, excludefornull) & 
                        !mixtrand) {
                      KSrandomprec <- max(abs(plgamma(rnprec, 
                                                      location=0,
                                                      scale=1/randomprecprev[2],
                                                      randomprecprev[1]) - 
                                              plgamma(rnprec,
                                                      location=0, 
                                                      scale=1/randomprec[2], 
                                                      randomprec[1]))
                      )
                    }
                    else {
                      KSrandomprec <- max(diracprobprev[2] * abs(plgamma(rnprec, 
                                                                         location=0, 
                                                                         scale=1/randomprecprev[2], 
                                                                         randomprecprev[1]) -
                                           diracprob[2] * plgamma(rnprec, 
                                                                  location=0, 
                                                                  scale=1/randomprec[2],
                                                                  randomprec[1]))
                      )
                    }
                    KS <- c(KSfixed,
                            KSmaxaddfixed, 
                            KSsigma,
                            KSrandomprec
                    )
                    names(KS) <- c("KSfixed", 
                                   "KSmaxaddfixed", 
                                   "KSsigma", 
                                   "KSrandomprec"
                    )
                    ksall <- rbind(ksall, KS)
                    KSmax <- max(KS[-length(KS)], na.rm=T)
                    pmlist <- list()
                    fixed <- as.numeric(fixed)
                    randomprec <- as.numeric(randomprec)
                    precerr <- as.numeric(precerr)
                    pmlist <- c(pmlist, list(mufixed=fixed[1], 
                                             precfixed=fixed[2]))
                    if (!is.null(shrinkaddfixed)) {
                      for (i in 1:els) {
                        addf <- as.numeric(addfixed[[i]])
                        name <- c(paste("mu", shrinkaddfixed[i], sep=""),
                                  paste("prec", shrinkaddfixed[i], sep=""))
                        toadd <- list(muaddfixed=addf[1], 
                                      precaddfixed=addf[2]
                        )
                        names(toadd) <- name
                        pmlist <- c(pmlist, toadd)
                      }
                    }
                    else {
                      pmlist <- c(pmlist, list(muaddfixed=addfixed[[1]][1], 
                                               precaddfixed=addfixed[[1]][2])
                      )
                    }
                    pmlist <- c(pmlist, list(shaperand=randomprec[1], 
                                             raterand=randomprec[2],
                                             mixp=diracprob)
                    )
                    if (!is.null(shrinkaddrandom)) {
                      for (i in 1:elsr) {
                        addr <- as.numeric(addrandom[[i]])
                        name <- c(paste("shape", shrinkaddrandom[i], sep=""),
                                  paste("rate", shrinkaddrandom[i], sep="")
                        )
                        toadd <- list(shapeaddr=addr[1],
                                      rateaddr=addr[2]
                        )
                        names(toadd) <- name
                        pmlist <- c(pmlist, toadd)
                      }
                    }
                    else {
                      pmlist <- c(pmlist, list(shapeaddr=addrandomprec[1], 
                                               rateaddr=addrandomprec[2])
                      )
                    }
                    pmlist <- c(pmlist, list(shapeerr=precerr[1], 
                                             rateerr=precerr[2])
                    )
                    paranew <- c(unlist(pmlist),
                                 nfeat=ngenej, 
                                 meanmlik=mlikmean
                    )
                    paraall <- rbind(paraall, paranew)
                  }
                  print(paranew)
                  print(paranew - paraprev)
                  print(ksall)
                  iter <- iter + 1
                  if (iter > maxiter | (KSmax <= tol & KSrandomprec <= 
                                        tolrand) | mlikconv) 
                    moreiter <- FALSE
                  mlikprev <- mlikmean
                  paraprev <- paranew
                }
              }
              ret <- list(pmlist=pmlist,
                          ksall=ksall,
                          paraall=paraall, 
                          inputpar=inputpar,
                          addfixed=addfixed,
                          addrandom=addrandom, 
                          typelik="gaussian"
              )
              return(ret)
}
