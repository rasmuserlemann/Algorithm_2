# Code used in the power study and examples in the paper *Conditional Monte Carlo Revisited*
The code is written in R. Professors Bo Henry Lindqvist and Gunnar Taraldsen contributed to writing the code.

## Side ideas
Files Algorithm1exp, Algorithm1trunexp, Algorithm2.5trunexp are side ideas. These methods work for special cases. Like, for the truncated exponential or uniform distributions.

## Main idea
Prior density (function prior(theta) can be adjusted to aid with root finding. Also, some of the parameters in the halving algorithm can be adjusted, depending on the distribution we wish to sample from. When simulating from W_t(u)f(u), we use the Metropolis-Hastings algorithm. The standard uniform distribution is used as the proposal distribution.

This repo has the code for inverse Gaussian, Gamma and truncated exponential distributions.

In the end, we plot the empirical marginal distribution of the conditional samples and compare them with against another method. 

## Gibbs Sampler
File GibbsGamma is based on *Use of the Gibbs sampler to obtain conditional tests, with applications*. There are some numerical problems for some sufficient statistic values. These can be adjusted by varying the range of where we are looking for the root.
