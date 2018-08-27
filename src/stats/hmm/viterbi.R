viterbi <- function(Pi, p, A) {
    N <- nrow(p)
    K <- ncol(p)
    omega <- matrix(NA, nrow=N, ncol=K)
    psi <- matrix(NA, nrow=N, ncol=K)
    vpath <- rep(NA, times=N)

    omega[1,] <- log(Pi) + log(p[1,])
    for(n in 1:(N-1)) {
        for(k in 1:K) {
            vnk <- log(A[,k]) + omega[n,]
            jnk <- which.max(vnk)
            maxVnk <- vnk[jnk]
            psi[n+1,k] <- jnk
            omega[n+1,k] <- log(p[n+1,k])+maxVnk
        }
    }
    vpath[N] <- which.max(omega[N,])
    for(n in seq(from=N-1, to=1, by=-1)) {
        vpath[n] <- psi[n+1, vpath[n+1]]
    }
    return(vpath)
}
