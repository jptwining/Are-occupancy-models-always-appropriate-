# This script contains functions and code for a simulation which explores impact of model type on inference from count data where a species is widespread and abundant (naive occupancy ~1)
# Here we simulate data from a Poisson-Poisson N-mixture model to represent sampling of a widely occuring population using passive detectors e.g., camera traps
# where individuals are unidentifiable and can leave more than 1 detection each over the sampling period.

# Author: Joshua P. Twining
# Date: 08/22/2025


# load packages
library('nimble')
library('coda')

# Function to simulate count data from a Poisson-Poisson N-mixture model
data.fn <- function(I = 200, J = 3, beta0 = 3.5, beta1 = 1, beta2 = 2,
                   alpha0 = -1, alpha1 = 1, alpha2 = -1){
  
  # set coefficents for parameters (on logs-odd scale, use plogis(x) to see value on probability-scale, e.g. for beta0.psi/mean occupancy)
  # I = 100
  # J = 3
  # beta0 = 3
  # beta1 = 1 
  # beta2 = 2 
  # alpha0 = 0 
  # alpha1 = 1 
  # alpha2 = -1
  
# simulate some site covariates
forest <- rnorm(I) # simulates canopy cover values at I sites from a normal distribution with mean=0 and sd =1 (scaled)
grassland <- rnorm(I)
road <- rnorm(I)

# simulate some site-by-observation covariates
temperature <-array(runif(n=I*J, -2, 2), dim =c(I,J)) # simulates temp values at I sites at J plots from a uniform distribution from -2 to 2
rainfall <- array(runif(n=I*J, -2, 2), dim =c(I,J)) 

# calculate lambda using a log-link
lambda <- matrix(ncol = I)  
for (i in 1:I){  
  lambda[i] <- exp(beta0 + beta1 * forest[i] + beta2 * grassland[i])
}


# calculate N states (abundance N, drawn from Poisson distribution)
N <- matrix(ncol = I)
for (i in 1:I){
  N[i] <- rpois(1, lambda[i])
}

# calculate detection rate (theta)
#3D matrix (i=site, j=plot, k=occasion)
theta = array(0, dim = c(I, J))
for (i in 1:I){
  for (j in 1:J){
      theta[i,j] <- exp(alpha0 + alpha1 * temperature[i,j] + alpha2 * rainfall[i,j])
    }
  }


# simulate observations (sampling) at site i, plot j, and occasion k, as a bernoulli/binomial trail of probability p, conditional on use of the plot 
y = array(0, dim = c(I, J))
for (i in 1:I){
  for (j in 1:J){
      y[i,j] <- (rpois(1, theta[i,j]*N[i]))
    }
  }

# quantize counts into detection/non-detection data for occupancy model
yoccu <- y
yoccu[yoccu>1]<- 1

return(list(I = I, J = J, lambda = lambda, beta0 = beta0, beta1 = beta1, beta2 = beta2, theta = theta, alpha0 = alpha0, alpha1 = alpha1,
            alpha2 = alpha2, N = N, y = y, yoccu=yoccu, forest = forest, grassland = grassland, road = road, temperature = temperature, 
            rainfall = rainfall))
}

# use function to simulate dataset with default parameter values
data <- data.fn()


# Prerun of all models to produce arrays of appropriate dimensions to save out estimates into

# # # # # # # # # # # # # # # # # # Poisson GLM # # # # # # # # # # # # # # # # # # 
PoissonGLMconstants <- list(I=data$I,J=data$J,forest=data$forest, grass = data$grassland)
PoissonGLMNimdata <- list(y=data$y)
PoissonGLMNiminits <- list(beta0 = runif(1),  beta1 = runif(1), beta2 = runif(1))
#set parameters to monitor
PoissonGLMparameters <- c("beta0","beta1", "beta2")

PoissonGLMNimModel <- nimbleCode({
  # priors
  beta0 ~ dnorm(0, sd= 5)
  beta1 ~ dnorm(0, sd = 5)
  beta2 ~ dnorm(0, sd = 5)

  # Likelihood
  for(i in 1:I){
    log(lambda[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i] 
    for(j in 1:J){
      y[i,j] ~ dpois(lambda[i]) 
    }
  }
})# end mode

# Build the model, configure the mcmc, and compileConfigure
start.time <- Sys.time()
PoissonGLMRmodel <- nimbleModel(code=PoissonGLMNimModel, constants=PoissonGLMconstants, data=PoissonGLMNimdata,check=FALSE,inits=PoissonGLMNiminits)
PoissonGLMconf <- configureMCMC(PoissonGLMRmodel,monitors=PoissonGLMparameters, thin=3, useConjugacy=FALSE)

# Build and compile
PoissonGLMRmcmc <- buildMCMC(PoissonGLMconf)
PoissonGLMCmodel <- compileNimble(PoissonGLMRmodel)
PoissonGLMCmcmc <- compileNimble(PoissonGLMRmcmc, project = PoissonGLMRmodel)

# Run the model
start.time2 <- Sys.time()
PoissonGLMCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

#get the chains
PoissonGLMmvSamples = as.matrix(PoissonGLMCmcmc$mvSamples)

burnin = 1000
iter = 10000

# # # # # # # # # # # # Negative binomial GLM # # # # # # # # # # # # # # # # # # 
NegbinGLMconstants <- list(I=data$I,J=data$J,forest=data$forest, grass = data$grassland)
NegbinGLMNimdata <- list(y=data$y)
NegbinGLMNiminits <- list(beta0 = runif(1),  beta1 = runif(1), beta2 = runif(1), r = runif(1))
#set parameters to monitor
NegbinGLMparameters <- c("beta0","beta1", "beta2")

NegbinGLMNimModel <- nimbleCode({
  # priors
  beta0 ~ dnorm(0, sd= 5)
  beta1 ~ dnorm(0, sd = 5)
  beta2 ~ dnorm(0, sd = 5)
  r ~ dgamma(0.01, 0.01)
  # Likelihood
  for(i in 1:I){
    log(mu[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i] 
    p[i] <- r / (r + mu[i])
    for(j in 1:J){
      y[i,j] ~ dnegbin(prob = p[i], size =r) 
    }
  }
})# end mode

# Build the model, configure the mcmc, and compileConfigure
start.time <- Sys.time()
NegbinGLMRmodel <- nimbleModel(code=NegbinGLMNimModel, constants=NegbinGLMconstants, data=NegbinGLMNimdata,check=FALSE,inits=NegbinGLMNiminits)
NegbinGLMconf <- configureMCMC(NegbinGLMRmodel,monitors=NegbinGLMparameters, thin=3, useConjugacy=FALSE)

# Build and compile
NegbinGLMRmcmc <- buildMCMC(NegbinGLMconf)
NegbinGLMCmodel <- compileNimble(NegbinGLMRmodel)
NegbinGLMCmcmc <- compileNimble(NegbinGLMRmcmc, project = NegbinGLMRmodel)

# Run the model
start.time2 <- Sys.time()
NegbinGLMCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

#get the chains
NegbinGLMmvSamples = as.matrix(NegbinGLMCmcmc$mvSamples)

burnin = 1000
iter = 10000

# # # # # # Binomial N-mixture model # # # # 
BinomialNmixconstants <- list(I=data$I,J=data$J,forest=data$forest, grass = data$grassland,  rain = data$rainfall, temp= data$temperature)
BinomialNmixNimdata <- list(y=data$y)

# create initial values for N
Nst <- matrix(NA, nrow = nrow(data$y))
for (i in 1:nrow(data$y)){
  Nst[i] <- max(data$y[i])
}

Nst <- as.vector(Nst)
BinomialNmixNiminits <- list(N = Nst, beta0 = runif(1),  beta1 = runif(1), beta2 = runif(1), alpha0 = runif(1), alpha1 = runif(1), alpha2 = runif(1))
#set parameters to monitor
BinomialNmixparameters <- c("beta0", "beta1", "beta2", "alpha0", "alpha1", "alpha2")
#
BinomialNmixNimModel <- nimbleCode({
  # priors
  beta0 ~ dnorm(0, sd= 5)
  beta1 ~ dnorm(0, sd = 5)
  beta2 ~ dnorm(0, sd = 5)
  alpha0 ~ dlogis(0, 1)
  alpha1 ~ dnorm(0, sd = 5)
  alpha2 ~ dnorm(0, sd = 5)

  # Binomial N mixture model
  # State model for N at site i
  for(i in 1:I){
    # estimate lambda as function  of covs
    log(lambda[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i]
    # estimate N from Poisson distribution
    N[i] ~ dpois(lambda[i])

  # Binomial observation model for detection of individuals at site i on onccasion j
    for(j in 1:J){
      y[i,j] ~ dbinom(p[i,j], N[i])
      logit(p[i,j]) <- alpha0 + alpha1 * rain[i,j] + alpha2 * temp[i,j]
    }
  }
})# end mode

## Build the model, configure the mcmc, and compileConfigure
start.time <- Sys.time()
BinomialNmixRmodel <- nimbleModel(code=BinomialNmixNimModel, constants=BinomialNmixconstants, data=BinomialNmixNimdata,check=FALSE,inits=BinomialNmixNiminits)
BinomialNmixconf <- configureMCMC(BinomialNmixRmodel,monitors=BinomialNmixparameters, thin=3, useConjugacy=FALSE)

# Build and compile
BinomialNmixRmcmc <- buildMCMC(BinomialNmixconf)
BinomialNmixCmodel <- compileNimble(BinomialNmixRmodel)
BinomialNmixCmcmc <- compileNimble(BinomialNmixRmcmc, project = BinomialNmixRmodel)

# Run the model
start.time2 <- Sys.time()
BinomialNmixCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

#get the chains
BinomialNmixmvSamples = as.matrix(BinomialNmixCmcmc$mvSamples)

burnin = 1000
iter = 10000


# # # # # # # # # # # # # # # # # # Poisson-Poisson N-mixture model # # # # # # # # # # # # # # # # # # 
Nst <- matrix(NA, nrow = nrow(data$y))
for (i in 1:nrow(data$y)){
  Nst[i] <- max(data$y[i])>1
}
Nst[Nst=="TRUE"] <- 1
Nst[Nst=="FALSE"] <-0
Nst <- as.vector(Nst)

PoissonNmixconstants <- list(I=data$I,J=data$J,forest=data$forest, grass = data$grassland,  rain = data$rainfall, temp= data$temperature)
PoissonNmixNimdata <- list(y=data$y)
PoissonNmixNiminits <- list(N = Nst, beta0 = runif(1),  beta1 = runif(1), beta2 = runif(1), alpha0 = runif(1), alpha1 = runif(1), alpha2 = runif(1))

#set parameters to monitor
PoissonNmixparameters <- c("beta0", "beta1", "beta2", "alpha0", "alpha1", "alpha2")

# Poisson-Poisson N-mixture model

PoissonNmixNimModel <- nimbleCode({
  # priors
  beta0 ~ dnorm(0, sd= 5)
  beta1 ~ dnorm(0, sd = 5)
  beta2 ~ dnorm(0, sd = 5)
  alpha0 ~ dnorm(0, sd= 5)
  alpha1 ~ dnorm(0, sd = 5)
  alpha2 ~ dnorm(0, sd = 5)
  
  # Likelihood
  # State model for N at site i
  for(i in 1:I){
    # estimate lambda as function  of covs
    log(lambda[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i]
    # estimate N from Poisson distribution
    N[i] ~ dpois(lambda[i])
    
    # Poisson observation model for detection of individuals at site i on on occasion j
    for(j in 1:J){
      y[i,j] ~ dpois(theta[i,j]*N[i])
      # estimate detection rate (theta) as a function of covariates
      log(theta[i,j]) <- alpha0 + alpha1 * rain[i,j] + alpha2 * temp[i,j]
    }
  }
})# end mode


# Build the model, configure the mcmc, and compileConfigure
start.time <- Sys.time()
PoissonNmixRmodel <- nimbleModel(code=PoissonNmixNimModel, constants=PoissonNmixconstants, data=PoissonNmixNimdata,check=FALSE,inits=PoissonNmixNiminits)
PoissonNmixconf <- configureMCMC(PoissonNmixRmodel,monitors=PoissonNmixparameters, thin=3, useConjugacy=FALSE)

# Build and compile
PoissonNmixRmcmc <- buildMCMC(PoissonNmixconf)
PoissonNmixCmodel <- compileNimble(PoissonNmixRmodel)
PoissonNmixCmcmc <- compileNimble(PoissonNmixRmcmc, project = PoissonNmixRmodel)

# Run the model
start.time2 <- Sys.time()
PoissonNmixCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

#get the chains
PoissonNmixmvSamples = as.matrix(PoissonNmixCmcmc$mvSamples)

burnin = 1000
iter = 10000


# # # # # # # # # # # # # # # # # # Occupancy model # # # # # # # # # # # # # # # # # # 
  occuconstants <- list(I=data$I,J=data$J,forest=data$forest, grass=data$grassland, rain = data$rainfall, temp= data$temperature)
  occuNimdata <- list(y=data$yoccu)
  
  # create initial values for Z state of site i
  Zst <- matrix(NA, nrow = nrow(data$yoccu))
  for (i in 1:nrow(data$yoccu)){
    Zst[i] <- max(data$yoccu[i])
  }
 
  Zst <- as.vector(Zst)
 
  # bundle up initial values
   occuNiminits <- list(beta0 = runif(1), alpha0=runif(1),  beta1 = runif(1), beta2 = runif(1), alpha1 = runif(1), alpha2= runif(1))
  
  #set parameters to monitor
  occuparameters <- c("beta0", "beta1", "beta2", "alpha0" ,"alpha1", 'alpha2')
  
  # Occupancy model
  occuNimModel <- nimbleCode({
    # priors
    beta0 ~ dlogis(0,1)
    alpha0 ~ dlogis(0,1)
    beta1 ~ dnorm(0, sd = 5)
    beta2 ~ dnorm(0, sd = 5)
    alpha1 ~ dnorm(0, sd = 5)
    alpha2 ~ dnorm(0, sd = 5)
    
    # likelihood for state model
    for(i in 1:I){
      # estimate psi as function  of covs
      logit(psi[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i]
      z[i] ~ dbern(psi[i]) #is site occupied

        # likelihood for observation model
        for(j in 1:J){
          logit(p[i,j]) <- alpha0 + alpha1 * temp[i,j] + alpha2 * rain[i,j]
          y[i,j] ~ dbern(p[i,j]*z[i]) #detection|site occupied
        }
      }
  })# end model
  
  # Build the model, configure the mcmc, and compileConfigure
  start.time <- Sys.time()
  occuRmodel <- nimbleModel(code=occuNimModel, constants=occuconstants, data=occuNimdata,check=FALSE,inits=occuNiminits)
  occuconf <- configureMCMC(occuRmodel,monitors=occuparameters, thin=3, useConjugacy=FALSE)
  
  # Build and compile
  occuRmcmc <- buildMCMC(occuconf)
  occuCmodel <- compileNimble(occuRmodel)
  occuCmcmc <- compileNimble(occuRmcmc, project = occuRmodel)
  
  # Run the model
  start.time2 <- Sys.time()
  occuCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
  end.time <- Sys.time()
  end.time - start.time  # total time for compilation, replacing samplers, and fitting
  end.time - start.time2 # post-compilation run time
  
  #get the chains
  occumvSamples = as.matrix(occuCmcmc$mvSamples)

  burnin = 5000
  iter = 10000


# Set up the simulations
  setwd('C:/...')
  
  i <- 250
  # set the number of cores to run in parrelel
  library(snow)
  library(doSNOW)
  library(foreach)
  cores=5
  cl.tmp = makeCluster(rep("localhost",cores), type="SOCK")
  registerDoSNOW(cl.tmp)
  
  # set number of sim reps and create arrays to stuff estimates into
  simrep <- 250   # Number of simulation replicates
  esti_occu <- array(NA, dim = c(4, (dim(occumvSamples)[2])))
  esti_PoissonGLM <- array(NA, dim = c(4, (dim(PoissonGLMmvSamples)[2])))
  esti_NegbinGLM <- array(NA, dim = c(4, (dim(NegbinGLMmvSamples)[2])))
  esti_PoissonNmix <- array(NA, dim = c(4, (dim(PoissonNmixmvSamples)[2])))
  esti_BinomialNmix <- array(NA, dim = c(4, (dim(BinomialNmixmvSamples)[2])))

  # name and label arrays
  colnames(esti_occu) <- colnames(occumvSamples)
  rownames(esti_occu)<- c("mean", "lowerCI","upperCI", "sd")
  colnames(esti_PoissonGLM) <- colnames(PoissonGLMmvSamples)
  rownames(esti_PoissonGLM)<- c("mean", "lowerCI","upperCI", "sd")
  colnames(esti_NegbinGLM) <- colnames(NegbinGLMmvSamples)
  rownames(esti_NegbinGLM)<- c("mean", "lowerCI","upperCI", "sd")
  colnames(esti_PoissonNmix) <- colnames(PoissonNmixmvSamples)
  rownames(esti_PoissonNmix)<- c("mean", "lowerCI","upperCI", "sd")
  colnames(esti_BinomialNmix) <- colnames(BinomialNmixmvSamples)
  rownames(esti_BinomialNmix)<- c("mean", "lowerCI","upperCI", "sd")

  # create array for true values
  true.vals <- array(NA, dim = c(6, 1))
  rownames(true.vals) <- c("beta0", "beta1", "beta2", "alpha0",   "alpha1",   "alpha2")
  
  # create function for searching for already run simulations
  `%!in%` = Negate(`%in%`)
  
  # Launch simulation
  out = foreach(rep=1:simrep) %dopar% {
    setwd('C:/...')
    files = list.files('C:/...')
    
    if(paste("count_data_model_sim",rep, ".RData", sep="") %!in% files){
      library(nimble)
      library(coda)
      cat(paste("\n\n*** Simrep Nr.", i, "***\n\n"))26
      
      # simulate data with coefs from uniform distribution from -1.5 to 1.5
      data <- data.fn(beta1 = runif(1, -1.5, 1.5),
                      beta2 = runif(1, -1.5, 1.5))

      # save out true vals from each sim rep
      true.vals[1,] <- data$beta0
      true.vals[2,] <- data$beta1
      true.vals[3,] <- data$beta2
      true.vals[4,] <- data$alpha0
      true.vals[5,] <- data$alpha1
      true.vals[6,] <- data$alpha2

      #model run parameters
      n.iter <- 30000

      #sampling parameters
      I = data$I
      J = data$J
      
      
      # run the chains - occupancy model
      n.chains = 3
      occuchains = vector("list", n.chains)
      # chains2 = vector("list", n.chains)
      for(chain in 1:n.chains){
        
        occuconstants <- list(I=data$I,J=data$J,forest=data$forest, grass=data$grassland, rain = data$rainfall, temp= data$temperature)
        occuNimdata <- list(y=data$yoccu)
        
        Zst <- matrix(NA, nrow = nrow(data$yoccu))
        for (i in 1:nrow(data$yoccu)){
          Zst[i] <- max(data$yoccu[i])
        }
        
        Zst <- as.vector(Zst)
        
        occuNiminits <- list(beta0 = runif(1), alpha0=runif(1),  beta1 = runif(1), beta2 = runif(1), alpha1 = runif(1), alpha2= runif(1))
        
        
        #set parameters to monitor
        occuparameters <- c("beta0", "beta1", "beta2", "alpha0" ,"alpha1", 'alpha2')
        
        # The model
        occuNimModel <- nimbleCode({
          # priors
          beta0 ~ dlogis(0,1)
          alpha0 ~ dlogis(0,1)
          beta1 ~ dnorm(0, sd = 5)
          beta2 ~ dnorm(0, sd = 5)
          alpha1 ~ dnorm(0, sd = 5)
          alpha2 ~ dnorm(0, sd = 5)
          
          # likelihood for state model
          for(i in 1:I){
            # estimate psi as function  of covs
            logit(psi[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i]
            z[i] ~ dbern(psi[i]) #is site occupied
            
            # likelihood for observation model
            for(j in 1:J){
              logit(p[i,j]) <- alpha0 + alpha1 * temp[i,j] + alpha2 * rain[i,j]
              y[i,j] ~ dbern(p[i,j]*z[i]) #detection|site occupied
            }
          }
        })# end model
        
        # Build the model, configure the mcmc, and compileConfigure
        start.time <- Sys.time()
        occuRmodel <- nimbleModel(code=occuNimModel, constants=occuconstants, data=occuNimdata,check=FALSE,inits=occuNiminits)
        occuconf <- configureMCMC(occuRmodel,monitors=occuparameters, thin=3, useConjugacy=FALSE)
        
        # Build and compile
        occuRmcmc <- buildMCMC(occuconf)
        occuCmodel <- compileNimble(occuRmodel)
        occuCmcmc <- compileNimble(occuRmcmc, project = occuRmodel)
        
        # Run the model
        start.time2 <- Sys.time()
        occuCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
        end.time <- Sys.time()
        end.time - start.time  # total time for compilation, replacing samplers, and fitting
        end.time - start.time2 # post-compilation run time
  
        #get the chains
       occumvSamples = as.matrix(occuCmcmc$mvSamples)
       occuchains[[chain]]=occumvSamples
  
}

        
# run the chains - Poisson GLM
n.chains = 3
PoissonGLMchains = vector("list", n.chains)
# chains2 = vector("list", n.chains)
for(chain in 1:n.chains){
  
  # fit the model
  PoissonGLMconstants <- list(I=data$I,J=data$J,forest=data$forest, grass = data$grassland)
  PoissonGLMNimdata <- list(y=data$y)
  PoissonGLMNiminits <- list(beta0 = runif(1),  beta1 = runif(1), beta2 = runif(1))
  
  #set parameters to monitor
  PoissonGLMparameters <- c("beta0","beta1", "beta2")
  
  PoissonGLMNimModel <- nimbleCode({
    # priors
    beta0 ~ dnorm(0, sd= 5)
    beta1 ~ dnorm(0, sd = 5)
    beta2 ~ dnorm(0, sd = 5)
    # Likelihood
    for(i in 1:I){
      log(lambda[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i] 
      for(j in 1:J){
        y[i,j] ~ dpois(lambda[i]) 
      }
    }
  })# end mode
  
  # Build the model, configure the mcmc, and compileConfigure
  start.time <- Sys.time()
  PoissonGLMRmodel <- nimbleModel(code=PoissonGLMNimModel, constants=PoissonGLMconstants, data=PoissonGLMNimdata,check=FALSE,inits=PoissonGLMNiminits)
  PoissonGLMconf <- configureMCMC(PoissonGLMRmodel,monitors=PoissonGLMparameters, thin=3, useConjugacy=FALSE)
  
  # Build and compile
  PoissonGLMRmcmc <- buildMCMC(PoissonGLMconf)
  PoissonGLMCmodel <- compileNimble(PoissonGLMRmodel)
  PoissonGLMCmcmc <- compileNimble(PoissonGLMRmcmc, project = PoissonGLMRmodel)
  
  # Run the model
  start.time2 <- Sys.time()
  PoissonGLMCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
  end.time <- Sys.time()
  end.time - start.time  # total time for compilation, replacing samplers, and fitting
  end.time - start.time2 # post-compilation run time
  
  #get the chains
  PoissonGLMmvSamples = as.matrix(PoissonGLMCmcmc$mvSamples)

  PoissonGLMchains[[chain]]=PoissonGLMmvSamples

}

# run the chains - Negbin GLM
n.chains = 3
NegbinGLMchains = vector("list", n.chains)
for(chain in 1:n.chains){
  
  # bundle the data
  NegbinGLMconstants <- list(I=data$I,J=data$J,forest=data$forest, grass = data$grassland)
  NegbinGLMNimdata <- list(y=data$y)
  
  # set initial values
  NegbinGLMNiminits <- list(beta0 = runif(1),  beta1 = runif(1), beta2 = runif(1), r = runif(1))
  
  #set parameters to monitor
  NegbinGLMparameters <- c("beta0","beta1", "beta2")
  
  # the model
  NegbinGLMNimModel <- nimbleCode({
    # priors
    beta0 ~ dnorm(0, sd= 5)
    beta1 ~ dnorm(0, sd = 5)
    beta2 ~ dnorm(0, sd = 5)
    r ~ dgamma(0.01, 0.01)
    # Likelihood
    for(i in 1:I){
      log(mu[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i] 
      p[i] <- r / (r + mu[i])
      for(j in 1:J){
        y[i,j] ~ dnegbin(prob = p[i], size =r) 
      }
    }
  })# end mode
  
  # Build the model, configure the mcmc, and compileConfigure
  start.time <- Sys.time()
  NegbinGLMRmodel <- nimbleModel(code=NegbinGLMNimModel, constants=NegbinGLMconstants, data=NegbinGLMNimdata,check=FALSE,inits=NegbinGLMNiminits)
  NegbinGLMconf <- configureMCMC(NegbinGLMRmodel,monitors=NegbinGLMparameters, thin=3, useConjugacy=FALSE)
  
  # Build and compile
  NegbinGLMRmcmc <- buildMCMC(NegbinGLMconf)
  NegbinGLMCmodel <- compileNimble(NegbinGLMRmodel)
  NegbinGLMCmcmc <- compileNimble(NegbinGLMRmcmc, project = NegbinGLMRmodel)
  
  # Run the model
  start.time2 <- Sys.time()
  NegbinGLMCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
  end.time <- Sys.time()
  end.time - start.time  # total time for compilation, replacing samplers, and fitting
  end.time - start.time2 # post-compilation run time
  
  #get the chains
  NegbinGLMmvSamples = as.matrix(NegbinGLMCmcmc$mvSamples)
  NegbinGLMchains[[chain]]=NegbinGLMmvSamples

}

# run the chains - Poisson N-mixture
n.chains = 3
PoissonNmixchains = vector("list", n.chains)
# chains2 = vector("list", n.chains)
for(chain in 1:n.chains){

  # bundle the data
  PoissonNmixconstants <- list(I=data$I,J=data$J,forest=data$forest, grass = data$grassland,  rain = data$rainfall, temp= data$temperature)
  PoissonNmixNimdata <- list(y=data$y)
  
  # set initial values
  Nst <- matrix(NA, nrow = nrow(data$y))
  for (i in 1:nrow(data$y)){
    Nst[i] <- max(data$y[i])>1
  }
  Nst[Nst=="TRUE"] <- 1
  Nst[Nst=="FALSE"] <-0
  Nst <- as.vector(Nst)
  
  PoissonNmixNiminits <- list(N = Nst, beta0 = runif(1),  beta1 = runif(1), beta2 = runif(1), alpha0 = runif(1), alpha1 = runif(1), alpha2 = runif(1))
  
  #set parameters to monitor
  PoissonNmixparameters <- c("beta0", "beta1", "beta2", "alpha0" ,"alpha1", 'alpha2')
  
  # the model
  PoissonNmixNimModel <- nimbleCode({
    # priors
    beta0 ~ dnorm(0, sd= 5)
    beta1 ~ dnorm(0, sd = 5)
    beta2 ~ dnorm(0, sd = 5)
    alpha0 ~ dnorm(0, sd= 5)
    alpha1 ~ dnorm(0, sd = 5)
    alpha2 ~ dnorm(0, sd = 5)
    
    # Likelihood
    # State model for N at site i
    for(i in 1:I){
      # estimate lambda as function  of covs
      log(lambda[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i]
      # estimate N from Poisson distribution
      N[i] ~ dpois(lambda[i])
      
      # Poisson observation model for detection of individuals at site i on on occasion j
      for(j in 1:J){
        y[i,j] ~ dpois(theta[i,j]*N[i])
        # estimate detection rate theta as a function of covariates
        log(theta[i,j]) <- alpha0 + alpha1 * rain[i,j] + alpha2 * temp[i,j]
      }
    }
  })# end mode
  
  
  # Build the model, configure the mcmc, and compileConfigure
  start.time <- Sys.time()
  PoissonNmixRmodel <- nimbleModel(code=PoissonNmixNimModel, constants=PoissonNmixconstants, data=PoissonNmixNimdata,check=FALSE,inits=PoissonNmixNiminits)
  PoissonNmixconf <- configureMCMC(PoissonNmixRmodel,monitors=PoissonNmixparameters, thin=6, useConjugacy=FALSE)
  
  # Build and compile
  PoissonNmixRmcmc <- buildMCMC(PoissonNmixconf)
  PoissonNmixCmodel <- compileNimble(PoissonNmixRmodel)
  PoissonNmixCmcmc <- compileNimble(PoissonNmixRmcmc, project = PoissonNmixRmodel)
  
  # Run the model
  start.time2 <- Sys.time()
  PoissonNmixCmcmc$run(, ,reset=FALSE) #Can keep extending the run by rerunning this line
  end.time <- Sys.time()
  end.time - start.time  # total time for compilation, replacing samplers, and fitting
  end.time - start.time2 # post-compilation run time
  
  #get the chains
  PoissonNmixmvSamples = as.matrix(PoissonNmixCmcmc$mvSamples)
  PoissonNmixchains[[chain]]=PoissonNmixmvSamples
}


# run the chains - Binomial N-mixture
n.chains = 3
BinomialNmixchains = vector("list", n.chains)
# chains2 = vector("list", n.chains)
for(chain in 1:n.chains){
  # bundle the data
  BinomialNmixconstants <- list(I=data$I,J=data$J,forest=data$forest, grass = data$grassland,  rain = data$rainfall, temp= data$temperature)
  BinomialNmixNimdata <- list(y=data$y)
  
  # set initial vals
  BinomialNmixNiminits <- list(N = Nst, beta0 = runif(1),  beta1 = runif(1), beta2 = runif(1), alpha0 = runif(1), alpha1 = runif(1), alpha2 = runif(1))
  
  #set parameters to monitor
  BinomialNmixparameters <- c("beta0", "beta1", "beta2", "alpha0" ,"alpha1", 'alpha2')
  
  # the model
  BinomialNmixNimModel <- nimbleCode({
    # priors
    beta0 ~ dnorm(0, sd= 5)
    beta1 ~ dnorm(0, sd = 5)
    beta2 ~ dnorm(0, sd = 5)
    alpha0 ~ dlogis(0, 1)
    alpha1 ~ dnorm(0, sd = 5)
    alpha2 ~ dnorm(0, sd = 5)
    
    # Likelihood
    # State model for N at site i
    for(i in 1:I){
      # estimate lambda as function  of covs
      log(lambda[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i]
      # estimate N from Poisson distribution
      N[i] ~ dpois(lambda[i])
      
      # Binomial observation model for detection of individuals at site i on onccasion j
      for(j in 1:J){
        y[i,j] ~ dbinom(p[i,j], N[i])
        logit(p[i,j]) <- alpha0 + alpha1 * rain[i,j] + alpha2 * temp[i,j]
      }
    }
  })# end mode
  
  ## Build the model, configure the mcmc, and compileConfigure
  start.time <- Sys.time()
  BinomialNmixRmodel <- nimbleModel(code=BinomialNmixNimModel, constants=BinomialNmixconstants, data=BinomialNmixNimdata,check=FALSE,inits=BinomialNmixNiminits)
  BinomialNmixconf <- configureMCMC(BinomialNmixRmodel,monitors=BinomialNmixparameters, thin=6, useConjugacy=FALSE)
  
  # Build and compile
  BinomialNmixRmcmc <- buildMCMC(BinomialNmixconf)
  BinomialNmixCmodel <- compileNimble(BinomialNmixRmodel)
  BinomialNmixCmcmc <- compileNimble(BinomialNmixRmcmc, project = BinomialNmixRmodel)
  
  # Run the model
  start.time2 <- Sys.time()
  BinomialNmixCmcmc$run(60000,reset=FALSE) #Can keep extending the run by rerunning this line
  end.time <- Sys.time()
  end.time - start.time  # total time for compilation, replacing samplers, and fitting
  end.time - start.time2 # post-compilation run time
  
  #get the chains
  BinomialNmixmvSamples = as.matrix(BinomialNmixCmcmc$mvSamples)
  BinomialNmixchains[[chain]]=BinomialNmixmvSamples
}
  
  
n.iter = 10000
n.burn = 5000

#combine the chains and burn
occumod=mcmc.list(mcmc(occuchains[[1]][n.burn:n.iter,]),
                     mcmc(occuchains[[2]][n.burn:n.iter,]),
                     mcmc(occuchains[[3]][n.burn:n.iter,]))

PoissonGLM=mcmc.list(mcmc(PoissonGLMchains[[1]][n.burn:n.iter,]),
                      mcmc(PoissonGLMchains[[2]][n.burn:n.iter,]),
                      mcmc(PoissonGLMchains[[3]][n.burn:n.iter,]))

NegbinGLM=mcmc.list(mcmc(NegbinGLMchains[[1]][n.burn:n.iter,]),
                     mcmc(NegbinGLMchains[[2]][n.burn:n.iter,]),
                     mcmc(NegbinGLMchains[[3]][n.burn:n.iter,]))

PoissonNmix=mcmc.list(mcmc(PoissonNmixchains[[1]][n.burn:n.iter,]),
                     mcmc(PoissonNmixchains[[2]][n.burn:n.iter,]),
                     mcmc(PoissonNmixchains[[3]][n.burn:n.iter,]))

BinomialNmix=mcmc.list(mcmc(BinomialNmixchains[[1]][n.burn:n.iter,]),
                      mcmc(BinomialNmixchains[[2]][n.burn:n.iter,]),
                      mcmc(BinomialNmixchains[[3]][n.burn:n.iter,]))

# gelman rubin diags
occugelman <- gelman.diag(occumod)
PoissonGLMgelman <- gelman.diag(PoissonGLM)
NegbinGLMgelman <- gelman.diag(NegbinGLM)
PoissonNmixgelman <- gelman.diag(PoissonNmix)
BinomialNmixgelman <- gelman.diag(BinomialNmix)

# combine the chains
occumod=runjags::combine.mcmc(occumod)
PoissonGLM=runjags::combine.mcmc(PoissonGLM)
NegbinGLM=runjags::combine.mcmc(NegbinGLM)
PoissonNmix=runjags::combine.mcmc(PoissonNmix)
BinomialNmix=runjags::combine.mcmc(BinomialNmix)

# save out estimates
esti_occu[1,1] <- mean(occumod[,"beta0"])
esti_occu[c(2,3),1]<- quantile(occumod[,"beta0"], probs = c(2.5,97.5)/100)
esti_occu[4,1]<- sd(occumod[,"beta0"])
esti_occu[1,2] <- mean(occumod[,"beta1"])
esti_occu[c(2,3),2]<- quantile(occumod[,"beta1"], probs = c(2.5,97.5)/100)
esti_occu[4,2]<- sd(occumod[,"beta1"])
esti_occu[1,3] <- mean(occumod[,"beta2"])
esti_occu[c(2,3),3]<- quantile(occumod[,"beta2"], probs = c(2.5,97.5)/100)
esti_occu[4,3]<- sd(occumod[,"beta2"])
esti_occu[1,4] <- mean(occumod[,"alpha0"])
esti_occu[c(2,3),4]<- quantile(occumod[,"alpha0"], probs = c(2.5,97.5)/100)
esti_occu[4,4]<- sd(occumod[,"alpha0"])
esti_occu[1,5] <- mean(occumod[,"alpha1"])
esti_occu[c(2,3),5]<- quantile(occumod[,"alpha1"], probs = c(2.5,97.5)/100)
esti_occu[4,5]<- sd(occumod[,"alpha1"])
esti_occu[1,6] <- mean(occumod[,"alpha2"])
esti_occu[c(2,3),6]<- quantile(occumod[,"alpha2"], probs = c(2.5,97.5)/100)
esti_occu[4,6]<- sd(occumod[,"alpha2"])

esti_PoissonGLM[1,1] <- mean(PoissonGLM[,"beta0"])
esti_PoissonGLM[c(2,3),1]<- quantile(PoissonGLM[,"beta0"], probs = c(2.5,97.5)/100)
esti_PoissonGLM[4,1]<- sd(PoissonGLM[,"beta0"])
esti_PoissonGLM[1,2] <- mean(PoissonGLM[,"beta1"])
esti_PoissonGLM[c(2,3),2]<- quantile(PoissonGLM[,"beta1"], probs = c(2.5,97.5)/100)
esti_PoissonGLM[4,2]<- sd(PoissonGLM[,"beta1"])
esti_PoissonGLM[1,3] <- mean(PoissonGLM[,"beta2"])
esti_PoissonGLM[c(2,3),3]<- quantile(PoissonGLM[,"beta2"], probs = c(2.5,97.5)/100)
esti_PoissonGLM[4,3]<- sd(PoissonGLM[,"beta2"])

esti_NegbinGLM[1,1] <- mean(NegbinGLM[,"beta0"])
esti_NegbinGLM[c(2,3),1]<- quantile(NegbinGLM[,"beta0"], probs = c(2.5,97.5)/100)
esti_NegbinGLM[4,1]<- sd(NegbinGLM[,"beta0"])
esti_NegbinGLM[1,2] <- mean(NegbinGLM[,"beta1"])
esti_NegbinGLM[c(2,3),2]<- quantile(NegbinGLM[,"beta1"], probs = c(2.5,97.5)/100)
esti_NegbinGLM[4,2]<- sd(NegbinGLM[,"beta1"])
esti_NegbinGLM[1,3] <- mean(NegbinGLM[,"beta2"])
esti_NegbinGLM[c(2,3),3]<- quantile(NegbinGLM[,"beta2"], probs = c(2.5,97.5)/100)
esti_NegbinGLM[4,3]<- sd(NegbinGLM[,"beta2"])

esti_PoissonNmix[1,1] <- mean(PoissonNmix[,"beta0"])
esti_PoissonNmix[c(2,3),1]<- quantile(PoissonNmix[,"beta0"], probs = c(2.5,97.5)/100)
esti_PoissonNmix[4,1]<- sd(PoissonNmix[,"beta0"])
esti_PoissonNmix[1,2] <- mean(PoissonNmix[,"beta1"])
esti_PoissonNmix[c(2,3),2]<- quantile(PoissonNmix[,"beta1"], probs = c(2.5,97.5)/100)
esti_PoissonNmix[4,2]<- sd(PoissonNmix[,"beta1"])
esti_PoissonNmix[1,3] <- mean(PoissonNmix[,"beta2"])
esti_PoissonNmix[c(2,3),3]<- quantile(PoissonNmix[,"beta2"], probs = c(2.5,97.5)/100)
esti_PoissonNmix[4,3]<- sd(PoissonNmix[,"beta2"])
esti_PoissonNmix[1,4] <- mean(PoissonNmix[,"alpha0"])
esti_PoissonNmix[c(2,3),4]<- quantile(PoissonNmix[,"alpha0"], probs = c(2.5,97.5)/100)
esti_PoissonNmix[4,4]<- sd(PoissonNmix[,"alpha0"])
esti_PoissonNmix[1,5] <- mean(PoissonNmix[,"alpha1"])
esti_PoissonNmix[c(2,3),5]<- quantile(PoissonNmix[,"alpha1"], probs = c(2.5,97.5)/100)
esti_PoissonNmix[4,5]<- sd(PoissonNmix[,"alpha1"])
esti_PoissonNmix[1,6] <- mean(PoissonNmix[,"alpha2"])
esti_PoissonNmix[c(2,3),6]<- quantile(PoissonNmix[,"alpha2"], probs = c(2.5,97.5)/100)
esti_PoissonNmix[4,6]<- sd(PoissonNmix[,"alpha2"])

esti_BinomialNmix[1,1] <- mean(BinomialNmix[,"beta0"])
esti_BinomialNmix[c(2,3),1]<- quantile(BinomialNmix[,"beta0"], probs = c(2.5,97.5)/100)
esti_BinomialNmix[4,1]<- sd(BinomialNmix[,"beta0"])
esti_BinomialNmix[1,2] <- mean(BinomialNmix[,"beta1"])
esti_BinomialNmix[c(2,3),2]<- quantile(BinomialNmix[,"beta1"], probs = c(2.5,97.5)/100)
esti_BinomialNmix[4,2]<- sd(BinomialNmix[,"beta1"])
esti_BinomialNmix[1,3] <- mean(BinomialNmix[,"beta2"])
esti_BinomialNmix[c(2,3),3]<- quantile(BinomialNmix[,"beta2"], probs = c(2.5,97.5)/100)
esti_BinomialNmix[4,3]<- sd(BinomialNmix[,"beta2"])
esti_BinomialNmix[1,4] <- mean(BinomialNmix[,"alpha0"])
esti_BinomialNmix[c(2,3),4]<- quantile(BinomialNmix[,"alpha0"], probs = c(2.5,97.5)/100)
esti_BinomialNmix[4,4]<- sd(BinomialNmix[,"alpha0"])
esti_BinomialNmix[1,5] <- mean(BinomialNmix[,"alpha1"])
esti_BinomialNmix[c(2,3),5]<- quantile(BinomialNmix[,"alpha1"], probs = c(2.5,97.5)/100)
esti_BinomialNmix[4,5]<- sd(BinomialNmix[,"alpha1"])
esti_BinomialNmix[1,6] <- mean(BinomialNmix[,"alpha2"])
esti_BinomialNmix[c(2,3),6]<- quantile(BinomialNmix[,"alpha2"], probs = c(2.5,97.5)/100)
esti_BinomialNmix[4,6]<- sd(BinomialNmix[,"alpha2"])


#put all the sim outputs together and save it.
setwd('C:/Users/twininjo/Documents/R/Modelling_count_data_simulations/I=200')
allthebits <- list(esti_occu, esti_PoissonGLM, esti_NegbinGLM, esti_PoissonNmix, esti_BinomialNmix, true.vals, occugelman,
                  PoissonGLMgelman, NegbinGLMgelman, PoissonNmixgelman, BinomialNmixgelman)
save(allthebits, file = paste("count_data_model_sim",rep, ".RData", sep=""))

rm(chains)
gc()
    }
  }
  stopCluster(cl.tmp)

  

