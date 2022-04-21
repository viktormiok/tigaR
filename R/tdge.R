###############################################################################
#
#         function for testing temporal differential gene expression
#
###############################################################################

tdge <- function(dat, 
                 dat1=NULL,
                 fitAlt,
                 fitNull,
                 design, 
                 ZSpline,
                 multivar=FALSE){
        # fitAlt: Fit object of the alternative hypothesis
        # fitNull: Fit object of the null hypothesis
        # dat: Gene expression data object
        # multivar: Test of multivariate fitted object
        if (!is(dat, "matrix")) {
          stop("Input (dat) is of wrong class.")
        }
        if (!is(fitAlt, "list")) {
          stop("Input (fitAlt) is of wrong class.")
        }
        if (!is(fitNull, "list")) {
          stop("Input (fitNull) is of wrong class.")
        }
        if (length(fitAlt[[1]]) != length(fitNull[[1]])) {
          stop("Dimension of fitted object do not correspond")
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
        if (!is(multivar, "logical")) {
          stop("Input (multivar) is of wrong class.")
        }
        ratio = logLik = numeric()
        for (i in 1:nrow(dat)){
          obj <- fitAlt[[1]][[i]]$.args$data
          if(multivar){
            obj$fitA <- fitAlt[[1]][[i]]$multiFit
            obj$fit0 <- fitNull[[1]][[i]]$multiFit
          } else {
            obj$fitA <- fitAlt[[1]][[i]]$summary.fitted.values$mean
            obj$fit0 <- fitNull[[1]][[i]]$summary.fitted.values$mean
          }
          if(is.null(obj)){
            ratio <- rbind(ratio, 1)
            logLik <- rbind(logLik, 1)
          }
          else{
            lam <- -2*log(sqrt((sum((obj$y - obj$fitA)^2)/
                                  sum((obj$y - obj$fit0)^2))^ncol(dat))
            ) 
            df <- suppressWarnings(degFreedom(fit=fitAlt[[1]][[i]], 
                                              design=design,
                                              ZSpline=ZSpline
            )) 
            pval <- 1- pchisq(lam, df)   
            ratio <- rbind(ratio,
                           pval
            )       
            logLik <- rbind(logLik, 
                            lam
            )      
          }
        }
        cor_ratio <- p.adjust(ratio, 
                              method="fdr", 
                              n=length(ratio)
        )
        output <- data.frame(cbind(rownames(dat), logLik, ratio, cor_ratio))
        colnames(output) <- c("GeneName", "logLik", "pvalue", "padj")
        return(output)
}
