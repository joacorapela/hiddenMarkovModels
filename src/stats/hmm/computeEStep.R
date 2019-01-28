
computeEStep <- function(x, Pi, A, phi, getEmissionProbabilities) {
    p <- getEmissionProbabilities(x=x, phi=phi)
    res <- getAlphaHat(Pi=Pi, A=A, p=p)
    alphaHat <- res$alphaHat
    c <- res$c
    betaHat <- getBetaHat(A=A, p=p, c=c)
    gamma <- alphaHat*betaHat # (13.64)
    answer <- list(gamma=gamma, p=p, c=c, alphaHat=alphaHat, betaHat=betaHat)
    return(answer)
}

