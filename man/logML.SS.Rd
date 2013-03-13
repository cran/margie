\name{logML.SS}
\alias{logML.SS}
\title{Log marginal likelihood of a mixture of Gaussian distributions}
\description{Use the Stepping Stone Sampler of Xie et al. (2011) to calculate the log marginal likelihood of a Bayesian Gaussian Mixture Model}
\usage{logML.SS(final.samples = 10,
		thinning.q = 10,
		num.mix.comp = 2,
		write.posterior.chain = FALSE,
    out.dir = 0,
		mu.prior.Mean = 0,
		mu.prior.Var = 0.5,
		K = 30,
		alpha = 0.3,
		model.Var = c(1, 1),
		tuning.d2 = 0.2,
		epsilon = 0.01,
		NormalMix.observations = rnorm(n = 50, mean = 0, sd = 1),
		omega.prior.DirPar = 1)}
\arguments{
  \item{final.samples}{(>=1) number of MCMC iterations, after thinning, to keep from target distribution g}
  \item{thinning.q}{(>=1) parameter to thin the MCMC chain. The actual number of iterations that the MCMC sampler will produce is thinning.q*final.samples;}
  \item{num.mix.comp}{(>=2) number of mixture components to fit to data points}
  \item{write.posterior.chain}{(TRUE/FALSE) whether the program should write or not the posterior chain of sampled values to files stored within 'out.dir'}
  \item{out.dir}{output directory}
  \item{mu.prior.Mean, mu.prior.Var}{(-inf,inf),(>0) location and scale parameters, respectively, of prior on mu. Data points are fitted with a mixture: Data ~ Normal(mu_1, model.Var_1) + ... + Normal(mu_num.mix.comp, model.Var_num.mix.comp) with prior on mu_j ~ Normal(mu.prior.Mean, mu.prior.Var) and model.Var_1,...,model.Var_num.mix.comp assumed known}
  \item{K}{(>=1) number of factors in which ratio c_{1.0}/c_{0.0} is split into. As a rule of thumb, the greater K the better the estimation of the marginal likelihood. If K=1, the SS sampler will only sample from the posterior (burn-in phase) and the prior (marginal-likelihood estimation phase)}
  \item{alpha}{(>0) positive parameter of the Beta distribution. Xie et al. (2011) recommend setting alpha=0.3}
  \item{model.Var}{(>0) vector of scale parameters of Gaussian mixture model. Data points are fitted with a mixture: Data ~ Normal(mu_1, model.Var) + ... + Normal(mu_num.mix.comp, model.Var) with model.Var assumed known and common to all components}
  \item{tuning.d2}{(>0) tuning parameter of the 'sliding window with Normal proposal' mechanism to update the location parameter of mixture component j, mu_j. The greater tuning.d2, the larger the updating steps throughout the MCMC run}
  \item{epsilon}{(~0) correction parameter of the 'e-shifted Dirichlet proposal' used to update the vector of mixture weights, omega = {omega_1, ..., omega_num.mix.comp}. Set epsilon to a value close to zero}
  \item{NormalMix.observations}{vector of observed data}
  \item{omega.prior.DirPar}{(>=1) parameter of prior on omega ~ Dirichlet(omega.prior.DirPar). omega.prior.DirPar=1 indicates prior ignorance on component sizes. The larger omega.prior.DirPar, the stronger the prior belief that components must be of equal size}

}
\value{The Value that is produced on the R Console is the calculated log marginal likelihood and the sequence of values, B.powers, that are used to power-up the likelihood function along the stepping-stone path from the posterior (i.e. the first stepping stone, B.powers=1) to the prior (i.e. the last stepping stone, B.powers=0).

 If argument \code{write.posterior.chain} is set to TRUE, three output files will be written to the directory specified by argument \code{out.dir}. Each file name will have a random number assigned, to distinguish them (this value is displayed on the R console immediately when run). If two numbers coincide, the new data will be appended onto the older files, so it is best to specify a clean out.dir. The names of the output files also include the number of mixture components.

  The output files are:

  \item{mu.chain}{This file logs the MCMC chain of component means sampled from the posterior distribution (i.e. only during the first stepping stone).}
  \item{omega.chain}{This file logs the MCMC chain of component weights sampled from the posterior distribution.}
  \item{zet.chain}{This file logs the MCMC chain of allocation variables sampled from the posterior distribution.}
}
\section{Warning}{This computation can be quite slow on a poor computer.

  \code{system.time} gave a mean of 81.08 seconds for the \code{logML.SS} run in the example (see below) while, for reference, taking a mean of 6.71 seconds for the R code \code{for(i in 1:5000000){j = i}} on the same computer.

  Unless you have a powerful computer, low values for arguments \code{final.samples} and \code{K} are advisable for large amounts of data. Bear in mind, however, that low values for both of these parameters could make the result less reliable. The higher they are, the more accurate. The example lowers both values, and thus is quicker than most other executions, at the cost of accuracy.
}
\references{
 Kass, R.E. and Raftery, A.E. (1995). Bayes Factors. Journal of the American Statistical Association 90, 773-795.
 
 Xie, W. et al. (2011). Improving Marginal Likelihood Estimation for Bayesian Phylogenetic Model Selection, Systematic Biology, 60: 150-160.
}
\author{Elisa Loza, David Phillips}
\examples{
## Generate 25 Random Samples from a 2-Component Gaussian Mixture 
## Distribution
x = c(rnorm(15, -5, 1), rnorm(10, 5, 1))

## Plot the histogram of the data
hist(x,breaks=20)

## Run the Log-Marginal Likelihood Calculator without writing 
## any output file
   # Fit the data with two components
k.two <- logML.SS(NormalMix.observations = x, model.Var = c(1, 1), 
                  K = 5, final.samples = 10, num.mix.comp = 2, 
		    mu.prior.Mean = mean(x), mu.prior.Var = var(x))   

   # Fit the same data with three components
k.three <- logML.SS(NormalMix.observations = x, 
		      model.Var = c(0.1,0.1,0.1), K = 5, 
		      final.samples = 10, num.mix.comp = 3, 
		      mu.prior.Mean = mean(x), mu.prior.Var = var(x))

   # Compare model fit. Kass and Raftery's (1995) guidelines say 
   # that 2(Log.ML.1 - Log.ML.2) > 10 indicates very strong evidence 
   # in favour of model 1; 0 <= 2(Log.ML.1 - Log.ML.2) <= 2 indicates 
   # no significant difference between the models, and a range of 
   # evidence levels in between these two extremes can be consulted 
   # in Kass and Raftery (1995).
k.two$Log.ML - k.three$Log.ML     
}
\keyword{Marginal Likelihood}
\keyword{Bayesian Gaussian Mixture}
\keyword{Stepping Stone Sampler}