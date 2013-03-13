# Rothamsted Research 'margie' R package

# --------------------------------------------------------------------------------------------

# DESCRIPTION:
# This package calculates the log marginal likelihood of a model using the Stepping Stone (SS) sampler by Xie et al. 2011. 
# Arguments:
#   final.samples:          (>=1) number of MCMC iterations, after thinning, to keep from target distribution g;
#   thinning.q:             (>=1) parameter to thin the MCMC chain. The actual number of iterations that the MCMC sampler will produce is thinning.q*final.samples;
#   *num.mix.comp:          (>=2) number of mixture components to fit to data points;
#    write.posterior.chain  (TRUE/FALSE) whether the program should write or not the posterior chain of sampled values to files stored within 'out.dir' 
#   *out.dir:		            output directory;
#   *mu.prior.Mean:         (-inf,inf) location parameter of prior on mu. Data points are fitted with a mixture data ~ Normal(mu_1, model.Var_1) + ... + Normal(mu_num.mix.comp, model.Var_num.mix.comp) with prior on mu_j ~ Normal(mu.prior.Mean, mu.prior.Var) and model.Var_1,...,model.Var_num.mix.comp assumed known;
#   mu.prior.Var:           (>0) scale parameter of prior on mu. Data points are fitted with a mixture data ~ Normal(mu_1, model.Var_1) + ... + Normal(mu_num.mix.comp, model.Var_num.mix.comp) with prior on mu_j ~ Normal(mu.prior.Mean, mu.prior.Var) and model.Var_1,...,model.Var_num.mix.comp assumed known;
#   K:		                  (>=1) number of factors in which ratio c_{1.0}/c_{0.0} is split into. As a rule of thumb, the greater K the better the estimation of the marginal likelihood. If K=1, the SS sampler will only sample from the posterior (burn-in phase) and the prior (marginal-likelihood estimation phase);
#   alpha:                  (>0)  positive parameter of the Beta distribution. Xie et al. (2011) recommend setting alpha=0.3.
#   *model.Var:             (>0) vector of scale parameters of Normal mixture model; model.Var=(model.Var_1,...,model.Var_num.mix.comp). Data points are fitted with a mixture data ~ Normal(mu_1, model.Var_1) + ... + Normal(mu_num.mix.comp, model.Var_num.mix.comp) with model.Var_1,...,model.Var_num.mix.comp assumed known;
#   tuning.d2:              (>0) tuning parameter of the 'sliding window with Normal proposal' mechanism to update the location parameter of mixture component j, mu_j. The greater tuning.d2, the larger the updating steps throughout the MCMC run;
#   epsilon:                (~0) correction parameter of the 'e-shifted Dirichlet proposal' used to update the vector of mixture weights, omega = {omega_1, ..., omega_num.mix.comp}. Set epsilon to a value close to zero, e.g. 0.0001;
#   *NormalMix.observations: vector of observed data with 'num.observ' data points;
#   omega.prior.DirPar:     (>=1) parameter of prior on omega ~ Dirichlet(omega.prior.DirPar). omega.prior.DirPar=1 indicates prior ignorance on component sizes. The larger omega.prior.DirPar, the stronger the prior belief that components must be of equal size.
#   * indicates arguments that need tuning depending on the problem
#      !! WARNING -- directory into which output will be written must be clean !!

# For more information contact: elisa.loza@rothamsted.ac.uk

# ---------------------------------------------------------------------------------------------

# Copyright (C) Rothamsted Research, 2011-2013.

# Licensed under the GNU General Public License version 3.0 (the "License").

# You may obtain a copy of the License at <http://www.gnu.org/licenses/gpl-3.0.txt>.

# This program is free software: you can redistribute it and/or modify it under the terms 
# of the GNU General Public License as published by the Free Software Foundation, either 
# version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. 
# If not, see <http://www.gnu.org/licenses/gpl-3.0.txt>.

# Date: 13/03/2013

# ---------------------------------------------------------------------------------------------

  
##########################
#      Main function     #
##########################

logML.SS <- function(	final.samples=10,
			thinning.q=10,
			num.mix.comp=2,
			write.posterior.chain=FALSE,
      out.dir=0,
			mu.prior.Mean=0,
			mu.prior.Var=0.5,
			K=30,
			alpha=0.3,
			model.Var=c(1,1),
			tuning.d2=0.2,
			epsilon=0.01,
			NormalMix.observations=rnorm(n=50,mean=0,sd=1),
			omega.prior.DirPar=1				    )
{

  ss.original <- TRUE    # (TRUE/FALSE) if TRUE (Xie et al. 2011) or if FALSE (then, the Generalised SS method is implemented; Fan et al. 2011)

  # Sample identifier for run
  identifier <- sample.int(n=10000, size=1)
  if (write.posterior.chain == TRUE) {
    if (out.dir == 0) {
      cat("ERROR: argument 'out.dir' must be set to the directory where the output will be written\n")
      return()
    } else {
      cat("The identifier of the run is: ", identifier, " and the output files will be written onto directory", out.dir, "\n") 
    }  
  }
  
	# Number of observed data points
	num.observ <- length(NormalMix.observations)

  # Set beta powers: sequence of K+1 values, each one in [0,1], used to power up the likelihood function; betha.powers = {betha_K=1, ..., betha_2=(2/K)^(1/a), betha_1=(1/K)^(1/a), betha_0=0};
  betha.powers <- set.beta(K, alpha)

	# Calculate M given thinning.q and final.samples
	M <- thinning.q * final.samples

	# Starting values for parameters omega, zet and mu:
	omega.0 <- rep(1/num.mix.comp, num.mix.comp)          		# Set all weights to 1/num.mix.comp
	zet.0   <- which(rmultinom(n=num.observ, size=1, prob=omega.0) == 1, arr.ind=TRUE)[,"row"]
                                       					# Allocation i (i=1,...,num.observ) is sampled from a
                                       					# multinomial distribution with the probability of falling into
                                       					# each component given by the starting omega vector, omega.0
	mu.0    <- rnorm(n=num.mix.comp, mean=mu.prior.Mean, sd=sqrt(mu.prior.Var))
                                       					# Mean for component j (j=1,...,num.mix.comp) is independently sampled
                                       					# from the prior for mu_j (equation (G))

#mu.0 <- c(79,520,928,1286)
  
	vector.ln.r  <- c()
	
	# Create progress bar
	pb <- txtProgressBar(min=1, max=(K+1), style=3)
	# Update GUI console
	flush.console()
	
	for(k in 1:(K+1)) 
	{
	    Sys.sleep(0.1)
		    setTxtProgressBar(pb, k)

	    # 1. Set the target distribution as g =(prop.to) likelihood^{b_k} x prior
	    #    (This is done from within the MH.NormalMixture function)
    
	    # 2. Obtain M samples from g
	    returned.list <- MH.NormalMixture(ss.original, betha.powers[k], M, num.mix.comp, omega.0, zet.0, mu.0, out.dir, thinning.q, model.Var, mu.prior.Mean, mu.prior.Var, num.observ, tuning.d2, epsilon, NormalMix.observations, omega.prior.DirPar, identifier, write.posterior.chain)
	    # returned.list[[1]]: chain.ln.likelihood during MCMC run (after thinning it)
	    # returned.list[[2]]: omega^(M/q) = (omega_1^(M/q), ..., omega_k^(M/q))
	    # returned.list[[3]]: zet^(M/q)   = (zet_1^(M/q), ..., zet_N^(M/q))
	    # returned.list[[4]]: mu^(M/q)    = (mu_1^(M/q), ..., mu_k^(M/q))
	    
	    # 3. Set the starting parameter values of the next MCMC run to the last values sampled in the current run
	    omega.0 <- returned.list[[2]]
	    zet.0   <- returned.list[[3]]
	    mu.0    <- returned.list[[4]]
        
	    if(k > 1) {    # If k==1 then g = posterior and we don't require those samples
	        # 4. Factor out the maximum of all likelihood values stored
	        max.ln.lik <- max(returned.list[[1]][-1]) # The first element in chain.ln.likelihood is "NA" and I don't want to consider it
    
   		# 5. Calculate ln.r 
	        #     5.1. First, calculate summation ** Equation (m1) but see page 10 of notes, equation (**) too
	        summation <- 0
	        for(i in 2:length(returned.list[[1]])) {    # For loop runs from 2 because the first element in chain.ln.likelihood (i.e. returned.list[[1]][1]) is "NA"
		        summation <- summation + exp((betha.powers[k-1]-betha.powers[k])*(returned.list[[1]][i]-max.ln.lik))
	        }
    
    		#    5.2. Now calculate ** Equation (M) but see page 10 of notes, equation (**) too
	        ln.r <- ((betha.powers[k-1]-betha.powers[k])*(max.ln.lik))+log((1/length(returned.list[[1]]))*summation)     
  
  	        # 6. Store ln.r
	        vector.ln.r <- c(vector.ln.r, ln.r)
    	    }
	} 
	close(pb)  

	# Return
	   # write(sum(vector.ln.r), file=paste(out.dir,paste("marginal.lik.k",num.mix.comp,".",identifier,".txt",sep=""),sep="/")) 	  
	
  cat("\n", "Log.ML: log marginal likelihood of model", "\n", "B.powers: sequence of power values for likelihood", "\n", "\n")
	
  return(list(Log.ML = sum(vector.ln.r), B.powers = betha.powers))
		# vector.ln.r = {ln.r_{k=K}, ln.r_{k=K-1}, ..., ln.r_{k=1}}, i.e. in reverse order
		# An estimate of the log marginal likelihood is provided by:
		#       sum(vector.ln.r)
}
#########################


  
#########################
#  Auxiliary functions  #
#########################

# To sample from a Dirichlet distribution: draw x1,...,xk from independent gamma distributions
# with common scale and shape parameters a1,...,ak and for each j, let tj=xj/Sum_{i=1}^k xi.
# Then, t ~ Dirichlet(a1, ..., ak). This function returns vector 't' with k elements
# n: number of Dirichlet realisations, i.e. number of t vectors
# k: length of vector t: t=(t1,...,tk) where t ~ Dirichlet(alpha)
# alpha: vector of size k: alpha=(a1,...,ak)
rdirichlet <- function(n, k, alpha)
{
  dirichs <- c()
  for(i in 1:n) {
    gammas  <- c()
    for(i in 1:k) {
      gammas[i] <- rgamma(n=1, shape=alpha[i], scale=alpha[i])
    }
    dirichs <- cbind(dirichs, gammas/sum(gammas), deparse.level=0)
  }
  return(dirichs)
}

# Normal univariate density
ln.normal.univariate <- function(Mean, Var, Observation)
{
  ln.dens <- -((0.5*log(2*pi*Var))+((Observation-Mean)^2/(2*Var)))  # Equation (N)
  return(ln.dens)
}

# Likelihood function ** Equations (D) and (N) **
# model.Mean: vector of k component means mu=(mu_1,...,mu_k)
# model.Var: vector of k component variances sigma2=(sigma2_1,...,sigma2_k)
# allocs: vector of allocation variables z=(z_1,...,z_N)
# num.observ: number of observations N
calc.ln.likelihood <- function(model.Mean, model.Var, num.observ, vector.data, allocs)
{
  ln.likelihood <- 0
  for(i in 1:num.observ) {
    ln.likelihood <- ln.likelihood + ln.normal.univariate(model.Mean[allocs[i]], model.Var[allocs[i]], vector.data[i])
  }

  return(ln.likelihood)
}

# Prior ratio for mu evaluated at "value"
calc.ln.prior.mu <- function(value, mu.prior.Mean, mu.prior.Var)
{
  # Prior probability is assumed N(mu.prior.Mean, mu.prior.Var) *** CHANGE AS REQUIRED ***
  ln.prior <- ln.normal.univariate(mu.prior.Mean, mu.prior.Var, value)
  return(ln.prior)
}

# Log target ratio for zi move, according to Xie et al. 2011 target: g=(prop.to) likelihood^{b_{k-1}}.prior
# current: current allocation zi
# proposed: proposed allocation zi*
# omega: vector of current weights w=(w1,...,wk)
ln.target.ratio.Original.zet <- function(current, proposed, omega, ln.sampling.dist.current, ln.sampling.dist.proposed, betha)
{
  # Equation (J):
  return(log(omega[proposed])-log(omega[current])+betha*(ln.sampling.dist.proposed-ln.sampling.dist.current))
}

# Log target ratio for mu_j move, according to Xie et al. 2011 target: g =(prop.to) likelihood^{b_{k-1}}.prior
ln.target.ratio.Original.mu <- function(current, proposed, betha, ln.lik.proposed, ln.lik.current,mu.prior.Mean,mu.prior.Var)
{
  ln.prior.proposed  <- calc.ln.prior.mu(proposed, mu.prior.Mean, mu.prior.Var)
  ln.prior.current   <- calc.ln.prior.mu(current, mu.prior.Mean, mu.prior.Var)

  # Equation (K):
  return(betha*(ln.lik.proposed-ln.lik.current)+ln.prior.proposed-ln.prior.current)
}

# Log target ratio for omega move, according to Xie et al. 2011 target: g =(prop.to) likelihood^{b_{k-1}}.prior
# N: vector N=(N1,...,Nk) where Nj (j=1,...,k) is the number of observations allocated to component j
# current: vector of current weights w=(w1,...,wk)
# proposed: vector of proposed weights w*=(w1*,...,wk*)
ln.target.ratio.Original.omega <- function(current, proposed, N, k, omega.prior.DirPar)
{
  summation.proposed <- 0
  for(j in 1:k) {
    summation.proposed <- summation.proposed + (N[j]+omega.prior.DirPar-1)*log(proposed[j])
  }
  
  summation.current  <- 0
  for(j in 1:k) {
    summation.current <- summation.current + (N[j]+omega.prior.DirPar-1)*log(current[j])
  }
  
  # Equation (I1):
  return(summation.proposed - summation.current)
}

# Log proposal ratio for omega move
ln.proposal.ratio.omega <- function(current, proposed, N, k, epsilon,omega.prior.DirPar)
{
  summation.proposed <- 0
  for(j in 1:k) {
    summation.proposed <- summation.proposed + (N[j]+omega.prior.DirPar+epsilon-1)*log(proposed[j])
  }

  summation.current  <- 0
  for(j in 1:k) {
    summation.current <- summation.current + (N[j]+omega.prior.DirPar+epsilon-1)*log(current[j])
  }

  # Equation (I2):
  return(summation.current - summation.proposed)
}

# Acceptance probability
calc.ln.acceptance <- function(ln.target.ratio, ln.proposal.ratio)
{
  ln.alpha <- min(0, (ln.target.ratio + ln.proposal.ratio))
  return(ln.alpha)
}

# ss.original: TRUE (Xie et al. 2011) or FALSE (then, the Generalised SS method is implemented; Fan et al. 2011)
# k: number of mixture components
# betha: power of power-posterior; betha=0: the target is the prior; betha=1: the target is the posterior
# (starting.omega, starting.zet, starting.mu): starting parameter values for MCMC run
# out.dir: name of directory to write output to. Default: the current working directory in R
# Parameters of inferential interest: omega (vector of k mixture proportions), zet (vector of 'num.observ' allocations) and mu (vector of k Normal means)
MH.NormalMixture <- function(   ss.original,
				betha, 
				MCMC.its, 
				k,
				starting.omega, 
				starting.zet, 
				starting.mu, 
				out.dir, 
				thinning.q, 
				model.Var,
				mu.prior.Mean,
				mu.prior.Var, 
				num.observ,
				tuning.d2,
				epsilon,
				NormalMix.observations,
				omega.prior.DirPar,
        identifier,
        write.posterior.chain 	)
{
  # Declare parameters and likelihood chains
  omega.chain  <- list(starting.omega)
  zet.chain    <- list(starting.zet)
  mu.chain     <- list(starting.mu)
  ln.lik.chain <- c(NA)      # Likelihood of starting sampled value is missing (i.e. set to "NA") --cannot see how to avoid this

  for(n in 2:(MCMC.its+1)) {   # Start MCMC loop
    #print(paste("Iteration ", n, sep="")) 
    # 1. Set (omega.curr, zet.curr, mu.curr)=(omega, zet, mu)^(n-1)
    omega.curr <- omega.chain[[n-1]]
    zet.curr   <- zet.chain[[n-1]]
    mu.curr    <- mu.chain[[n-1]]

    #cat("Mu current: ", mu.curr, "\n")
    #cat("Omega current: ", omega.curr, "\n")
    #cat("Zet current: ", zet.curr, "\n")
    
    # 2. Generate candidate parameter values and accept/reject.
    #    ** UPDATE MEANS **
    #    2.1. Update all mu_j (j=1,...,k) using the 'sliding window with Normal proposal' mechanism:
    #             mu_j* ~ N(mu.curr_j, d2), where d2 is a tuning parameter
    #component  <- sample.int(n=k, size=1)
    for(component in 1:k) {
      #print(paste("Update mean for component ", component, sep=""))
      mu.star    <- mu.curr   # Copy first, then substitute only the jth component (given by 'component' above)
      mu.star[component] <- rnorm(1, mean=mu.curr[component], sd=sqrt(tuning.d2))

      # Accept mu_j*, that is set mu_j^(n)=mu_j*, with probability alpha (equation (K)):
      ln.lik.proposed <- calc.ln.likelihood(mu.star, model.Var, num.observ, NormalMix.observations, zet.curr)
      ln.lik.current  <- calc.ln.likelihood(mu.curr, model.Var, num.observ, NormalMix.observations, zet.curr)

      if(identical(ss.original, TRUE)){                 # Original SS method (Xie et al. 2011)
        ln.target.ratio <- ln.target.ratio.Original.mu(mu.curr, mu.star, betha, ln.lik.proposed, ln.lik.current, mu.prior.Mean, mu.prior.Var)
      } else {
        #ln.target.ratio <- ln.target.ratio.Generalised()  # Generalised SS method (Fan et al. 2011)
      }
      ln.proposal.ratio <- log(1)


      ln.alpha <- calc.ln.acceptance(ln.target.ratio, ln.proposal.ratio)
      if(ln.alpha > 0) {
        ln.alpha <- 0
      }
      r <- runif(1, min=0, max=1)
  	  if (log(r) < ln.alpha) {  # i.e. alpha > 1
        #print("Accept proposed")
  		  mu.curr <- mu.star
      }      # Otherwise reject proposed
    }

    #    ** UPDATE ALL WEIGHTS **
    #    2.2. Update omega using an e-shifted Dirichlet proposal:
    #                   omega* ~ Dir_k (r+N1+e, ..., r+Nk+e),
    #       where e>0 and Nj=Sum_{i=1}^N I[zi=j] is the number of observations
    #       allocated to component j, and I[.] is the indicator function
    #print("Update weights")
    N <- rep.int(0, times=k)
    membership <- table(zet.curr)  
    for(entry in 1:length(membership)) {
      N[as.numeric(names(membership[entry]))]<-membership[entry]
    }     # N=(N1, ..., Nk)
    omega.star <- as.vector(rdirichlet(n=1, k, alpha=N+omega.prior.DirPar+epsilon))
    
    # Accept omega*, that is set omega^(n)=omega*, with probability alpha (equation (I)):
    if(identical(ss.original, TRUE)){                 # Original SS method (Xie et al. 2011)
      ln.target.ratio <- ln.target.ratio.Original.omega(omega.curr, omega.star, N, k, omega.prior.DirPar)
    } else {
      #ln.target.ratio <- ln.target.ratio.Generalised()  # Generalised SS method (Fan et al. 2011)
    }


    ln.proposal.ratio <- ln.proposal.ratio.omega(omega.curr, omega.star, N, k, epsilon, omega.prior.DirPar)  # Equation (I2)

    ln.alpha <- calc.ln.acceptance(ln.target.ratio, ln.proposal.ratio)
    if(ln.alpha > 0) {
      ln.alpha <- 0
    }
    r <- runif(1, min=0, max=1)
  	if (log(r) < ln.alpha) {  # i.e. alpha > 1
      #print("Accept proposed")
  		omega.curr <- omega.star
	  }   # Otherwise reject proposed

    #    ** UPDATE ALLOCATIONS **
    #    2.3. Update all allocations zi (i=1,...,num.observ) by drawing from the set
    #         {1,...,k}\{zi}, with all elements equally likely to be drawn.
    #observation  <- sample.int(n=num.observ, size=1)
    for(observation in 1:num.observ) {
      #print(paste("Update allocation for observation ", observation, sep=""))
      zet.star     <- zet.curr   # Copy first, then substitute only the ith element (given by 'observation' above)
      zet.star[observation] <- sample(c(1:k)[-zet.curr[observation]], size=1)
    
      # Accept zi*, that is set zi^(n)=zi*, with probability alpha (equation (J)):
      ln.normal.proposed <- ln.normal.univariate(mu.curr[zet.star[observation]], model.Var[zet.star[observation]], NormalMix.observations[observation])
      ln.normal.current  <- ln.normal.univariate(mu.curr[zet.curr[observation]], model.Var[zet.curr[observation]], NormalMix.observations[observation])

      if(identical(ss.original, TRUE)){                 # Original SS method (Xie et al. 2011)
        ln.target.ratio <- ln.target.ratio.Original.zet(zet.curr[observation], zet.star[observation], omega.curr, ln.normal.current, ln.normal.proposed, betha)
      } else {
        #ln.target.ratio <- ln.target.ratio.Generalised()  # Generalised SS method (Fan et al. 2011)
      }
      ln.proposal.ratio <- log(1)

      ln.alpha <- calc.ln.acceptance(ln.target.ratio, ln.proposal.ratio)
      if(ln.alpha > 0) {
        ln.alpha <- 0
      }
      r <- runif(1, min=0, max=1)
      if (log(r) < ln.alpha) {  # i.e. alpha > 1
        #print("Accept proposed")
  		  zet.curr <- zet.star
	    }      # Otherwise reject proposed
    }

    # Write values to chains
    omega.chain[[n]]  <- omega.curr
    zet.chain[[n]]    <- zet.curr
    mu.chain[[n]]     <- mu.curr
    ln.lik.chain[n]   <- ln.lik.current  #NOTE: Likelihood of starting sampled value is missing (i.e. set to NA)

    if ((betha != 1) && ((n-1)%%thinning.q == 0)) {
        # Calculate likelihood given current parameter values
        ln.lik.current <- calc.ln.likelihood(mu.curr, model.Var, num.observ, NormalMix.observations, zet.curr)
    }

    if ((write.posterior.chain == TRUE) && (betha == 1) && ((n-1)%%thinning.q == 0)) {
        # Write chains to files
        write(omega.chain[[n]], file=paste(out.dir,paste("omega.chain.k",k,".",identifier,".txt",sep=""),sep="/"), append=TRUE, ncolumns=k)
        write(zet.chain[[n]], file=paste(out.dir,paste("zet.chain.k",k,".",identifier,".txt",sep=""),sep="/"), append=TRUE, ncolumns=num.observ)
        write(mu.chain[[n]], file=paste(out.dir,paste("mu.chain.k",k,".",identifier,".txt",sep=""),sep="/"), append=TRUE, ncolumns=k)
    }
    
  }  # End MCMC loop

  last.omega <- omega.chain[[n]]
  last.zet   <- zet.chain[[n]]
  last.mu    <- mu.chain[[n]]

  return(list(ln.lik.chain, last.omega, last.zet, last.mu))  ## DOUBLE CHECK!!
    # [[1]]: chain.ln.likelihood during MCMC run (after thinning it)
    # [[2]]: last.omega: omega^(M/q) = (omega_1^(M/q), ..., omega_k^(M/q))
    # [[3]]: last.zet:   zet^(M/q)   = (zet_1^(M/q), ..., zet_N^(M/q))
    # [[4]]: last.mu:    mu^(M/q)    = (mu_1^(M/q), ..., mu_k^(M/q))

}  # End function


# Set beta powers: This function sets the values of beta (to power up the likelihood function in the Stepping Stone Sampler; see Xie et al. 2011) according to uniform quantiles of a Beta(alpha,1) distribution.
#			      In practice, this entails setting beta_n=(n/K)^3.33 (n=1,...,K). Xie et al. (2011) show that the accuracy of the SS sampler is optimal in a Gaussian model example when beta values are set in this way.
# Arguments:
#   K:       (>=1) NOTE: This value has to COINCIDE with the K value in logML.SS
#   alpha:   (>0)  positive parameter of the Beta distribution. Xie et al. (2011) recommend setting alpha=0.3.

set.beta <- function(K=100, alpha=0.3) 
{
	vector.betha <- c(0)
	for(i in 2:(K+1)) {
		vector.betha[i] <- ((i-1)/K)^(1/alpha)
	}

	# vector.betha = {betha_0=0, betha_1=(1/K)^(1/a), betha_2=(2/K)^(1/a),..., betha_K=1}
	# So, this vector has K+1 elements

	# Note: I want to return the reversed version of vector.betha
	return(rev(vector.betha))
}

#########################

