# This script contains functions and code for a simulation which explores impact on sample size on inference from applying an occupancy model vs. a Bernoulli GLM to simulated 
# detection/non-detection data

# Author: Joshua P. Twining
# Date: 08/22/2025


# This function simulates detection/non-detection data subject to imperfect detection at I sites for J occasions

data.fn <- function(I = 50, J = 3, beta0 = 0, beta1 = 1, beta2 = -1, beta3 = -1,
                   alpha0 = 0, alpha1 = 1, alpha2 = -1){
  # set coefficents for parameters (on logs-odd scale, use plogis(x) to see value on probability-scale, e.g. for beta0.psi/mean occupancy)

# simulate some site covariates
forest <- rnorm(I) # simulates canopy cover values at I sites from a normal distribution with mean=0 and sd =1 (scaled)
grassland <- rnorm(I)
road_density <- rnorm(I)

# simulate some site-by-observation covariates
temperature <-array(runif(n=I*J, -2, 2), dim =c(I,J)) # simulates temp values at I sites at J plots from a uniform distribution from -2 to 2
rainfall <- array(runif(n=I*J, -2, 2), dim =c(I,J)) 

# calculate psi using a logit-link
psi <- matrix(ncol = I)  
for (i in 1:I){  
  psi[i] <- plogis(beta0 + beta1 * forest[i] + beta2 * grassland[i] + beta3 * road_density[i])
}

# calculate true z states (bernoulli/binomial trail of probability psi)
z <- matrix(ncol = I)
for (i in 1:I){
  z[i] <- rbinom(1, 1, psi[i])
}

# calculate detection probability (p)
#3D matrix (i=site, j=plot, k=occasion)
p = array(0, dim = c(I, J))
for (i in 1:I){
  for (j in 1:J){
      p[i,j] <- plogis(alpha0 + alpha1 * temperature[i,j] + alpha2 * rainfall[i,j])
    }
  }


# simulate observations (sampling) at site i, plot j, and occasion k, as a bernoulli/binomial trail of probability p, conditional on use of the plot 
y = array(0, dim = c(I, J))
for (i in 1:I){
  for (j in 1:J){
      y[i,j] <- (rbinom(1, 1, p[i,j])*z[i])
    }
  }

str(y)

zsum <- sum(z)

return(list(I = I, J = J, psi = psi, beta0 = beta0, beta1 = beta1, beta2 = beta2, beta3=beta3, p = p, alpha0 = alpha0, alpha1 = alpha1,
            alpha2 = alpha2, z = z, y = y, forest = forest, grassland = grassland, road_density = road_density, temperature = temperature, 
            rainfall = rainfall))
}


# run the function with default values
data <- data.fn()


#prerun simulated dataset to get dimensions for estimates for each rep 
library(nimble)
library(coda)


  #fit full occupancy model
  fulloccuconstants <- list(I=data$I,J=data$J,forest=data$forest, grass=data$grassland, road = data$road_density, rain = data$rainfall, temp= data$temperature)
  fulloccuNimdata <- list(y=data$y)
  
  Zst <- matrix(NA, nrow = nrow(data$y))
  for (i in 1:nrow(data$y)){
    Zst[i] <- max(data$y[i])
  }
 
  Zst <- as.vector(Zst)
  fulloccuNiminits <- list(beta0 = runif(1), alpha0=runif(1),  beta1 = runif(1), beta2 = runif(1), beta3 = runif(1), alpha1 = runif(1), alpha2= runif(1))
  
  
  #set parameters to monitor
  fulloccuparameters <- c("beta0", "beta1", "beta2","beta3","alpha0" ,"alpha1", 'alpha2')
  
  
  fulloccuNimModel <- nimbleCode({
    # priors
    beta0 ~ dlogis(0,1)
    alpha0 ~ dlogis(0,1)
    beta1 ~ dnorm(0, sd = 5)
    beta2 ~ dnorm(0, sd = 5)
    beta3 ~ dnorm(0, sd = 5)
    alpha1 ~ dnorm(0, sd = 5)
    alpha2 ~ dnorm(0, sd = 5)
    
    # likelihood for state model
    for(i in 1:I){
      # estimate psi as function  of covs
      logit(psi[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i] + beta3 * road[i]
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
  fulloccuRmodel <- nimbleModel(code=fulloccuNimModel, constants=fulloccuconstants, data=fulloccuNimdata,check=FALSE,inits=fulloccuNiminits)
  fulloccuconf <- configureMCMC(fulloccuRmodel,monitors=fulloccuparameters, thin=3, useConjugacy=FALSE)
  
  # Build and compile
  fulloccuRmcmc <- buildMCMC(fulloccuconf)
  fulloccuCmodel <- compileNimble(fulloccuRmodel)
  fulloccuCmcmc <- compileNimble(fulloccuRmcmc, project = fulloccuRmodel)
  
  # Run the model
  start.time2 <- Sys.time()
  fulloccuCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
  end.time <- Sys.time()
  end.time - start.time  # total time for compilation, replacing samplers, and fitting
  end.time - start.time2 # post-compilation run time
  
  #get the chains
  fulloccumvSamples = as.matrix(fulloccuCmcmc$mvSamples)


  #fit null detection occupancy model
  pdotoccuconstants <- list(I=data$I,J=data$J,forest=data$forest, grass=data$grassland, road = data$road_density)
  pdotoccuNimdata <- list(y=data$y)
  pdotoccuNiminits <- list(beta0 = runif(1), alpha0=runif(1),  beta1 = runif(1), beta2 = runif(1), beta3 = runif(1))

  #set parameters to monitor
  pdotoccuparameters <- c("beta0", "beta1", "beta2","beta3","alpha0")


  # Null detection occupancy model
  pdotoccuNimModel <- nimbleCode({
    # priors
    beta0 ~ dlogis(0,1)
    alpha0 ~ dlogis(0,1)
    beta1 ~ dnorm(0, sd = 5)
    beta2 ~ dnorm(0, sd = 5)
    beta3 ~ dnorm(0, sd = 5)

    # likelihood for state model
    for(i in 1:I){
      # estimate psi as function  of covs
      logit(psi[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i] + beta3 * road[i]
      z[i] ~ dbern(psi[i]) #is site occupied

      # likelihood for observation model
      for(j in 1:J){
        logit(p[i,j]) <- alpha0
        y[i,j] ~ dbern(p[i,j]*z[i]) #detection|site occupied
      }
    }
  })# end model

  # Build the model, configure the mcmc, and compileConfigure
  start.time <- Sys.time()
  pdotoccuRmodel <- nimbleModel(code=pdotoccuNimModel, constants=pdotoccuconstants, data=pdotoccuNimdata,check=FALSE,inits=pdotoccuNiminits)
  pdotoccuconf <- configureMCMC(pdotoccuRmodel,monitors=pdotoccuparameters, thin=3, useConjugacy=FALSE)

  # Build and compile
  pdotoccuRmcmc <- buildMCMC(pdotoccuconf)
  pdotoccuCmodel <- compileNimble(pdotoccuRmodel)
  pdotoccuCmcmc <- compileNimble(pdotoccuRmcmc, project = pdotoccuRmodel)

  # Run the model
  start.time2 <- Sys.time()
  pdotoccuCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
  end.time <- Sys.time()
  end.time - start.time  # total time for compilation, replacing samplers, and fitting
  end.time - start.time2 # post-compilation run time

  #get the chains
  pdotoccumvSamples = as.matrix(pdotoccuCmcmc$mvSamples)

# run the chains - logistic regression model
  #fit model
  logregconstants <- list(I=data$I,J=data$J,forest=data$forest, grass=data$grassland, road = data$road_density)
  logregNiminits <- list(beta0 = runif(1),  beta1 = runif(1), beta2 = runif(1), beta3 = runif(1))
  logregNimdata <- list(y=data$y)
  
  #set parameters to monitor
  logregparameters <- c("beta0","beta1", "beta2","beta3")
  
  # logistic regression model
  logregNimModel <- nimbleCode({
    # priors
    beta0 ~ dlogis(0,1)
    beta1 ~ dnorm(0, sd = 5)
    beta2 ~ dnorm(0, sd = 5)
    beta3 ~ dnorm(0, sd = 5)

    # Bernoulli GLM
    for(i in 1:I){
      # estimate p as function  of covs
      logit(p[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i] + beta3 * road[i]
      for(j in 1:J){
      y[i,j] ~ dbern(p[i]) 
      }
    }
  })# end model
  
  # Build the model, configure the mcmc, and compileConfigure
  start.time <- Sys.time()
  logregRmodel <- nimbleModel(code=logregNimModel, constants=logregconstants, data=logregNimdata,check=FALSE,inits=logregNiminits)
  logregconf <- configureMCMC(logregRmodel,monitors=logregparameters, thin=3, useConjugacy=FALSE)
  
  # Build and compile
  logregRmcmc <- buildMCMC(logregconf)
  logregCmodel <- compileNimble(logregRmodel)
  logregCmcmc <- compileNimble(logregRmcmc, project = logregRmodel)
  
  # Run the model
  start.time2 <- Sys.time()
  logregCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
  end.time <- Sys.time()
  end.time - start.time  # total time for compilation, replacing samplers, and fitting
  end.time - start.time2 # post-compilation run time
  
  #get the chains
  logregmvSamples = as.matrix(logregCmcmc$mvSamples)

# Now we set up the simulations to run a bunch of times in paralell (change depending on number of cores available)
  library(snow)
  library(doSNOW)
  library(foreach)
  
  setwd('C:/...')

# set number of cores
  cores=5
  cl.tmp = makeCluster(rep("localhost",cores), type="SOCK")
  registerDoSNOW(cl.tmp)
  
# set Number of simulation replicates
  simrep <- 250 
  
  # create arrays to stuff estimates into each rep
  esti_fulloccu <- array(NA, dim = c(4, (dim(fulloccumvSamples)[2])))
  esti_pdotoccu <- array(NA, dim = c(4, (dim(pdotoccumvSamples)[2])))
  esti_logreg <- array(NA, dim = c(4, (dim(logregmvSamples)[2])))
  
  # provide column and rownames so things make sense on otherside
  colnames(esti_fulloccu) <- colnames(fulloccumvSamples)
  rownames(esti_fulloccu)<- c("mean", "lowerCI","upperCI", "sd")
  colnames(esti_pdotoccu) <- colnames(pdotoccumvSamples)
  rownames(esti_pdotoccu)<- c("mean", "lowerCI","upperCI", "sd")
  colnames(esti_logreg) <- colnames(esti_logreg)
  rownames(esti_logreg)<- c("mean", "lowerCI","upperCI", "sd")
  true.vals <- array(NA, dim = c(7, 1))
  rownames(true.vals) <- c("beta0", "beta1", "beta2", "beta3", "alpha0",   "alpha1",   "alpha2")
  
  # 'Not in' function to check folder to not reproduce already run simulations (for if sim run crashes and need to rerun)
  `%!in%` = Negate(`%in%`)
  
  # Launch simulation
  out = foreach(rep=1:simrep) %dopar% {
    setwd('C:/...')
    files = list.files('C:/...')
    
    if(paste("occu_log_reg_sim",rep, ".RData", sep="") %!in% files){
      library(nimble)
      library(coda)
      cat(paste("\n\n*** Simrep Nr.", i, "***\n\n"))
      
      # simulate data with slope coefficents drawn from values between -1.5 and 1.5
      data <- data.fn(beta1 = runif(1, -1.5, 1.5),
                      beta2 = runif(1, -1.5, 1.5),
                      beta3 = runif(1, -1.5, 1.5))
     
      # save out true values for each sim rep
      true.vals[1,] <- data$beta0
      true.vals[2,] <- data$beta1
      true.vals[3,] <- data$beta2
      true.vals[4,] <- data$beta3
      true.vals[5,] <- data$alpha0
      true.vals[6,] <- data$alpha1
      true.vals[7,] <- data$alpha2

      #model run parameters
      n.iter <- 30000

      #sampling parameters
      I = data$I
      J = data$J
      

# run the chains - full occupancy model
n.chains = 3
fulloccuchains = vector("list", n.chains)
# chains2 = vector("list", n.chains)
for(chain in 1:n.chains){
  
  #fit model
  fulloccuconstants <- list(I=data$I,J=data$J,forest=data$forest, grass=data$grassland, road = data$road_density, rain = data$rainfall, temp= data$temperature)
  fulloccuNimdata <- list(y=data$y)
  
  fulloccuNiminits <- list(beta0 = runif(1), alpha0=runif(1),  beta1 = runif(1), beta2 = runif(1), beta3 = runif(1), alpha1 = runif(1), alpha2= runif(1))
  
  
  #set parameters to monitor
  fulloccuparameters <- c("beta0", "beta1", "beta2","beta3","alpha0" ,"alpha1", 'alpha2')
  
  
  fulloccuNimModel <- nimbleCode({
    # priors
    beta0 ~ dlogis(0,1)
    alpha0 ~ dlogis(0,1)
    beta1 ~ dnorm(0, sd = 5)
    beta2 ~ dnorm(0, sd = 5)
    beta3 ~ dnorm(0, sd = 5)
    alpha1 ~ dnorm(0, sd = 5)
    alpha2 ~ dnorm(0, sd = 5)
    
    # likelihood for state model
    for(i in 1:I){
      # estimate psi as function  of covs
      logit(psi[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i] + beta3 * road[i]
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
  fulloccuRmodel <- nimbleModel(code=fulloccuNimModel, constants=fulloccuconstants, data=fulloccuNimdata,check=FALSE,inits=fulloccuNiminits)
  fulloccuconf <- configureMCMC(fulloccuRmodel,monitors=fulloccuparameters, thin=3, useConjugacy=FALSE)
  
  # Build and compile
  fulloccuRmcmc <- buildMCMC(fulloccuconf)
  fulloccuCmodel <- compileNimble(fulloccuRmodel)
  fulloccuCmcmc <- compileNimble(fulloccuRmcmc, project = fulloccuRmodel)
  
  # Run the model
  start.time2 <- Sys.time()
  fulloccuCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
  end.time <- Sys.time()
  end.time - start.time  # total time for compilation, replacing samplers, and fitting
  end.time - start.time2 # post-compilation run time
  
  #get the chains
  fulloccumvSamples = as.matrix(fulloccuCmcmc$mvSamples)
  
  fulloccuchains[[chain]]=fulloccumvSamples
  
}


# run the chains - pdot occupancy model
n.chains = 3
pdotoccuchains = vector("list", n.chains)
# chains2 = vector("list", n.chains)
for(chain in 1:n.chains){
  
  #fit model
  pdotoccuconstants <- list(I=data$I,J=data$J,forest=data$forest, grass=data$grassland, road = data$road_density)
  pdotoccuNimdata <- list(y=data$y)
  
  pdotoccuNiminits <- list(beta0 = runif(1), alpha0=runif(1),  beta1 = runif(1), beta2 = runif(1), beta3 = runif(1))
  
  
  #set parameters to monitor
  pdotoccuparameters <- c("beta0", "beta1", "beta2","beta3","alpha0")
  
  
  pdotoccuNimModel <- nimbleCode({
    # priors
    beta0 ~ dlogis(0,1)
    alpha0 ~ dlogis(0,1)
    beta1 ~ dnorm(0, sd = 5)
    beta2 ~ dnorm(0, sd = 5)
    beta3 ~ dnorm(0, sd = 5)
    
    
    # likelihood for state model
    for(i in 1:I){
      # estimate psi as function  of covs
      logit(psi[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i] + beta3 * road[i]
      z[i] ~ dbern(psi[i]) #is site occupied
      
      # likelihood for observation model
      for(j in 1:J){
        logit(p[i,j]) <- alpha0
        y[i,j] ~ dbern(p[i,j]*z[i]) #detection|site occupied
      }
    }
  })# end model
  
  # Build the model, configure the mcmc, and compileConfigure
  start.time <- Sys.time()
  pdotoccuRmodel <- nimbleModel(code=pdotoccuNimModel, constants=pdotoccuconstants, data=pdotoccuNimdata,check=FALSE,inits=pdotoccuNiminits)
  pdotoccuconf <- configureMCMC(pdotoccuRmodel,monitors=pdotoccuparameters, thin=3, useConjugacy=FALSE)
  
  # Build and compile
  pdotoccuRmcmc <- buildMCMC(pdotoccuconf)
  pdotoccuCmodel <- compileNimble(pdotoccuRmodel)
  pdotoccuCmcmc <- compileNimble(pdotoccuRmcmc, project = pdotoccuRmodel)
  
  # Run the model
  start.time2 <- Sys.time()
  pdotoccuCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
  end.time <- Sys.time()
  end.time - start.time  # total time for compilation, replacing samplers, and fitting
  end.time - start.time2 # post-compilation run time
  
  #get the chains
  pdotoccumvSamples = as.matrix(pdotoccuCmcmc$mvSamples)
  # mvSamples2 = as.matrix(Cmcmc$mvSamples2)
  
  pdotoccuchains[[chain]]=pdotoccumvSamples
  # chains2[[chain]]=mvSamples2
  
}

# run the chains - logistic regression model
n.chains = 3
logregchains = vector("list", n.chains)
# chains2 = vector("list", n.chains)
for(chain in 1:n.chains){
  
  #fit model
  logregconstants <- list(I=data$I,J=data$J,forest=data$forest, grass=data$grassland, road = data$road_density)
  logregNiminits <- list(beta0 = runif(1),  beta1 = runif(1), beta2 = runif(1), beta3 = runif(1))
  logregNimdata <- list(y=data$y)
  
  #set parameters to monitor
  logregparameters <- c("beta0","beta1", "beta2","beta3")
  
  
  logregNimModel <- nimbleCode({
    # priors
    beta0 ~ dlogis(0,1)
    beta1 ~ dnorm(0, sd = 5)
    beta2 ~ dnorm(0, sd = 5)
    beta3 ~ dnorm(0, sd = 5)
    
    # Bernoulli GLM
    for(i in 1:I){
      # estimate p as function  of covs
      logit(p[i]) <- beta0 + beta1 * forest[i] + beta2 * grass[i] + beta3 * road[i]
      for(j in 1:J){
        y[i,j] ~ dbern(p[i]) 
      }
    }
  })# end model
  
  # Build the model, configure the mcmc, and compileConfigure
  start.time <- Sys.time()
  logregRmodel <- nimbleModel(code=logregNimModel, constants=logregconstants, data=logregNimdata,check=FALSE,inits=logregNiminits)
  logregconf <- configureMCMC(logregRmodel,monitors=logregparameters, thin=3, useConjugacy=FALSE)
  
  # Build and compile
  logregRmcmc <- buildMCMC(logregconf)
  logregCmodel <- compileNimble(logregRmodel)
  logregCmcmc <- compileNimble(logregRmcmc, project = logregRmodel)
  
  # Run the model
  start.time2 <- Sys.time()
  logregCmcmc$run(30000,reset=FALSE) #Can keep extending the run by rerunning this line
  end.time <- Sys.time()
  end.time - start.time  # total time for compilation, replacing samplers, and fitting
  end.time - start.time2 # post-compilation run time
  
  #get the chains
  logregmvSamples = as.matrix(logregCmcmc$mvSamples)
  # mvSamples2 = as.matrix(Cmcmc$mvSamples2)
  
  logregchains[[chain]]=logregmvSamples
  # chains2[[chain]]=mvSamples2
  
}

n.iter = 10000
n.burn = 5000

#combine the chains and burn
fulloccumod=mcmc.list(mcmc(fulloccuchains[[1]][n.burn:n.iter,]),
                     mcmc(fulloccuchains[[2]][n.burn:n.iter,]),
                     mcmc(fulloccuchains[[3]][n.burn:n.iter,]))

pdotoccumod=mcmc.list(mcmc(pdotoccuchains[[1]][n.burn:n.iter,]),
                      mcmc(pdotoccuchains[[2]][n.burn:n.iter,]),
                      mcmc(pdotoccuchains[[3]][n.burn:n.iter,]))

logregmod=mcmc.list(mcmc(logregchains[[1]][n.burn:n.iter,]),
                      mcmc(logregchains[[2]][n.burn:n.iter,]),
                      mcmc(logregchains[[3]][n.burn:n.iter,]))
summary(fulloccumod)

fulloccugelman <- gelman.diag(fulloccumod)
pdotoccugelman <- gelman.diag(pdotoccumod)
logreggelman <- gelman.diag(logregmod)

fulloccumod=runjags::combine.mcmc(fulloccumod)
pdotoccumod=runjags::combine.mcmc(pdotoccumod)
logregmod=runjags::combine.mcmc(logregmod)

esti_fulloccu[1,1] <- mean(fulloccumod[,"beta0"])
esti_fulloccu[c(2,3),1]<- quantile(fulloccumod[,"beta0"], probs = c(2.5,97.5)/100)
esti_fulloccu[4,1]<- sd(fulloccumod[,"beta0"])
esti_fulloccu[1,2] <- mean(fulloccumod[,"beta1"])
esti_fulloccu[c(2,3),2]<- quantile(fulloccumod[,"beta1"], probs = c(2.5,97.5)/100)
esti_fulloccu[4,2]<- sd(fulloccumod[,"beta1"])
esti_fulloccu[1,3] <- mean(fulloccumod[,"beta2"])
esti_fulloccu[c(2,3),3]<- quantile(fulloccumod[,"beta2"], probs = c(2.5,97.5)/100)
esti_fulloccu[4,3]<- sd(fulloccumod[,"beta2"])
esti_fulloccu[1,4] <- mean(fulloccumod[,"beta3"])
esti_fulloccu[c(2,3),4]<- quantile(fulloccumod[,"beta3"], probs = c(2.5,97.5)/100)
esti_fulloccu[4,4]<- sd(fulloccumod[,"beta3"])
esti_fulloccu[1,5] <- mean(fulloccumod[,"alpha0"])
esti_fulloccu[c(2,3),5]<- quantile(fulloccumod[,"alpha0"], probs = c(2.5,97.5)/100)
esti_fulloccu[4,5]<- sd(fulloccumod[,"alpha0"])
esti_fulloccu[1,6] <- mean(fulloccumod[,"alpha1"])
esti_fulloccu[c(2,3),6]<- quantile(fulloccumod[,"alpha1"], probs = c(2.5,97.5)/100)
esti_fulloccu[4,6]<- sd(fulloccumod[,"alpha1"])
esti_fulloccu[1,7] <- mean(fulloccumod[,"alpha2"])
esti_fulloccu[c(2,3),7]<- quantile(fulloccumod[,"alpha2"], probs = c(2.5,97.5)/100)
esti_fulloccu[4,7]<- sd(fulloccumod[,"alpha2"])


esti_pdotoccu[1,1] <- mean(pdotoccumod[,"beta0"])
esti_pdotoccu[c(2,3),1]<- quantile(pdotoccumod[,"beta0"], probs = c(2.5,97.5)/100)
esti_pdotoccu[4,1]<- sd(pdotoccumod[,"beta0"])
esti_pdotoccu[1,2] <- mean(pdotoccumod[,"beta1"])
esti_pdotoccu[c(2,3),2]<- quantile(pdotoccumod[,"beta1"], probs = c(2.5,97.5)/100)
esti_pdotoccu[4,2]<- sd(pdotoccumod[,"beta1"])
esti_pdotoccu[1,3] <- mean(pdotoccumod[,"beta2"])
esti_pdotoccu[c(2,3),3]<- quantile(pdotoccumod[,"beta2"], probs = c(2.5,97.5)/100)
esti_pdotoccu[4,3]<- sd(pdotoccumod[,"beta2"])
esti_pdotoccu[1,4] <- mean(pdotoccumod[,"beta3"])
esti_pdotoccu[c(2,3),4]<- quantile(pdotoccumod[,"beta3"], probs = c(2.5,97.5)/100)
esti_pdotoccu[4,4]<- sd(pdotoccumod[,"beta3"])
esti_pdotoccu[1,5] <- mean(pdotoccumod[,"alpha0"])
esti_pdotoccu[c(2,3),5]<- quantile(pdotoccumod[,"alpha0"], probs = c(2.5,97.5)/100)
esti_pdotoccu[4,5]<- sd(pdotoccumod[,"alpha0"])


esti_logreg[1,1] <- mean(logregmod[,"beta0"])
esti_logreg[c(2,3),1]<- quantile(logregmod[,"beta0"], probs = c(2.5,97.5)/100)
esti_logreg[4,1]<- sd(logregmod[,"beta0"])
esti_logreg[1,2] <- mean(logregmod[,"beta1"])
esti_logreg[c(2,3),2]<- quantile(logregmod[,"beta1"], probs = c(2.5,97.5)/100)
esti_logreg[4,2]<- sd(logregmod[,"beta1"])
esti_logreg[1,3] <- mean(logregmod[,"beta2"])
esti_logreg[c(2,3),3]<- quantile(logregmod[,"beta2"], probs = c(2.5,97.5)/100)
esti_logreg[4,3]<- sd(logregmod[,"beta2"])
esti_logreg[1,4] <- mean(logregmod[,"beta3"])
esti_logreg[c(2,3),4]<- quantile(logregmod[,"beta3"], probs = c(2.5,97.5)/100)
esti_logreg[4,4]<- sd(logregmod[,"beta3"])


#put all the sim outputs together and save it.
setwd('C:/...')
allthebits <- list(esti_fulloccu, esti_pdotoccu, esti_logreg, true.vals, fulloccugelman, pdotoccugelman, logreggelman)
save(allthebits, file = paste("occu_log_reg_sim",rep, ".RData", sep=""))

rm(chains)
gc()
    }
  }
  stopCluster(cl.tmp)

