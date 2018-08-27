hiddenMarkovModels
==================

R code to estimate Hidden Markov Models

The main function is src/stats/hmm/emEstimationHMM.R I demonstrate the use of
this function to estimate hidden Markov models with multivariate Gaussian
(test/doTestMGaussian.R) and with Poisson observations. Poisson observations
are modeled without autohistory (test/doTestPoisson.R)
and with autohistory (test/doTestPoissonGLM.R).

doc/hmmPoissonGLM.pdf documents the evaluation of the estimation of Hidden
Markov Models with Poisson observations with autohistory in simulated data and
in multiunit activity count recordings from a subject with epilepsy.

