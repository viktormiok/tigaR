###############################################################################
#
#         plot of the tigaR models fit
#
###############################################################################

plot_tigaRfit <- function(fit, 
                          fit1, 
                          timefac,
                          groupfac,
                          lattice=FALSE,
                          cycle=FALSE,
                          title="title",
                          multivar=FALSE){
            # fit, fit1 : Object which contain fit, from FitIntAllShrink
            # Around0: Grafics will be ploted on the scale around 0
            # lattice: Make lattice plot
            
            if (!is(fit, "list")) {
              stop("Input (fit) is of wrong class.")
            }
            if (!is(fit1, "list")) {
              stop("Input (fit1) is of wrong class.")
            }
            if (length(fit[[1]]) != length(fit1[[1]])) {
              stop("Dimension of fitted object do not correspond")
            }
            if (!is(timefac, "integer")) {
              stop("Input (integer) is of wrong class.")
            }
            if (!is(groupfac, "factor")) {
              stop("Input (groupfac) is of wrong class.")
            }
            if (length(timefac) != length(groupfac)) {
              stop("Dimension of the contrasts do not correspond")
            }
            if (!is(lattice, "logical")) {
              stop("Input (lattice) is of wrong class.")
            }
            if (!is(cycle, "logical")) {
              stop("Input (cycle) is of wrong class.")
            }
            if (!is(title, "character")) {
              stop("Input (title) is of wrong class.")
            }
            if (!is(multivar, "logical")) {
              stop("Input (multivar) is of wrong class.")
            }
            numgroup <- as.numeric(max(levels(groupfac))) 
            data <- fit$.args$data$y
            if(multivar){
              fitval <- fit$multiFit
              fitval1 <- fit1$multiFit
            }else{
              fitval <- fit$summary.fitted.values$mean
              fitval1 <- fit1$summary.fitted.values$mean
            }
            
            
            dataSet <- cbind(data,
                             timefac,
                             groupfac
            )
            fitfin <- cbind(fitval,
                            fitval1,
                            groupfac
            )
            if(!cycle){  
              par(mfrow=c(1, numgroup))
              for(j in 1:numgroup){
                datgraf <- subset(dataSet, 
                                  groupfac == j
                )
                plot(datgraf[,2],
                     datgraf[,1],
                     main=paste("group",
                                  j,
                                  sep=" "), 
                     type="p", 
                     ylim=c(min(dataSet[,1]), 
                              max(dataSet[,1])),
                     ylab="Expression",
                     xlab="Time Points"
                )
                
                fitgraf <- subset(fitfin, 
                                  groupfac == j
                )
                lines(fitgraf[,1], 
                      type="l", 
                      col="red",
                      lwd=2
                )
                if(!is.null(fit1)) {
                  lines(fitgraf[,2],
                        col="blue",
                        lty="dashed"
                  )
                }   
              }
            } else {
              par(mfrow=c(1,1))
              fit_cyc <- data.frame(dataSet,
                                    fitfin
              )
              plot(fit_cyc$timefac, 
                   fit_cyc$data,
                   xlab="Time Points",
                   ylab="Expression",
                   main=title
              )
              fit_cyc <- fit_cyc[fit_cyc$groupfac == 1,]
              lines(fit_cyc$timefac,
                    fit_cyc$fitval, 
                    col="red"
              )
              lines(fit_cyc$timefac,
                    fit_cyc$fitval1,
                    col="blue", 
                    lty="dashed"
              )
            }
            if(lattice){     
              p <- data.frame(data, 
                              timefac, 
                              paste("group", groupfac, sep="_")
              )
              colnames(p) <- c("Expression", "TimePoints", "CellLine")
              fram <- xyplot(Expression ~ TimePoints | CellLine,
                             data=p,
                             layout=c(numgroup, 1)
              )                              
              print(fram)
              for(i in 1:numgroup){
                fitgraf <- subset(fitfin, 
                                  groupfac == i
                )
                trellis.focus("panel", i, 1)
                llines(fitgraf[,1],
                       col="red"
                )
                if(!is.null(fit1)){
                  llines(fitgraf[,2],
                         col="blue",
                         lty="dashed"
                  ) 
                }     
              }
              trellis.unfocus()
            }
}
