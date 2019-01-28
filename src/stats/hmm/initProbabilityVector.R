
initProbabilityVector <- function(K) {
    probabilityVector <- rep(NA, times=K)
    probabilityVectorBuilt <- FALSE
    while(TRUE) {
        probabilityVector[1:(K-1)] <- runif(n=K-1, min=0, max=2/K)
        sum <- sum(probabilityVector[1:(K-1)])
        if(sum<1) {
            probabilityVector[K] <- 1-sum
            return(probabilityVector)
        }
    }
}
