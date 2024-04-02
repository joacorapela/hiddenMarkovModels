
getARProbabilities <- function(x, phi, xHist) {

   getLogProbXn <- function(xn, xnHist, a, sigmaSq) {
       sum <- 0
       for(i in 1:length(xn)) {
           sum <- sum + getLogProb(xn=xn[i], xnHist=xnHist[i,], a=a,
                                             sigmaSq=sigmaSq)
       }
       return(sum)
   }
   getLogProb <- function(xn, xnHist, a, sigmaSq) { 
       answer <- -.5*(log(2*pi)+log(sigmaSq))-(xn-sum(a*c(1, xnHist)))^2/(2*sigmaSq)
       return(answer)
   }

    N <- ncol(x)
    L <- nrow(phi$a)
    K <- ncol(phi$a)
    p <- matrix(0, nrow=N, ncol=K)
    for(n in 1:N) {
        xnHist <- matrix(xHist[,n,],
                          nrow=dim(xHist)[1],
                          ncol=dim(xHist)[3])
        for(k in 1:K) {
            p[n,k] <- exp(getLogProbXn(xn=x[,n], xnHist=xnHist,
                                                a=phi$a[,k], 
                                                sigmaSq=phi$sigmaSq[k]))
        }
    }
    return(p)
}
