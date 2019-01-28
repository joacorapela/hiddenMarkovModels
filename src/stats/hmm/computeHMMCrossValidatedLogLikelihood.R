computeHMMCrossValidatedLogLikelihood <- function(x, Pi, A, phi,
                                                     getEmissionProbabilities,
                                                     percentageOfDataForInitialStateEstimation=NaN) {
    if(!is.nan(percentageOfDataForInitialStateEstimation)) {
        nSamplesForInitialStateEstimation <-
         ncol(x)*percentageOfDataForInitialStateEstimation
        res <- computeEStep(x=matrix(x[,1:nSamplesForInitialStateEstimation],
                                      nrow=nrow(x)),
                             Pi=Pi, A=A, phi=phi, 
                             getEmissionProbabilities=getEmissionProbabilities)
        Pi <- res$gamma[nrow(res$gamma),]
    }
    p <- getEmissionProbabilities(x=x, phi=phi)
    res <- getAlphaHat(Pi=Pi, A=A, p=p)
    c <- res$c
    logLikelihood <- sum(log(c)) # (13.63)
    return(logLikelihood)
}
