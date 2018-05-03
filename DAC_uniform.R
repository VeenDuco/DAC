# Copyright statement:
# The GNU GENERAL PUBLIC LICENSE v3.0 applies.

# Author comment
# no comment at this time

# File description comment
# This file is meant to calculate the Data Agreement Criterion (DAC) (Bousquet, 2008)
# for a 2 parameter (mean / sd) model. The DAC can be calculated for multiple priors 
# simultaniously such that we can rank and compare these priors. 
# Priors can be elicited experts' beliefs such that a ranking of experts can be obtained.
# In this specific function we calculate the DAC whilst using a Uniform benchmark prior. 

# loading required packages
library(blavaan)
library(flexmix)
library(sfsmisc)

# Function definitions
# Function to calculate the DAC for a 2 parameter model mean / sd. 
DAC.uniform <- function(from, to, by, data, priors, lb.bench, ub.bench, n.iter=10000){
  #
  # Args:
  #          from: Lower bound of the parameter space that is to be evaluated, as in the seq function of the base package
  #            to: Upper bound of the parameter space that is to be evaluated, as in the seq function of the base package  
  #            by: Step size by which the defined parameter space is mapped out, as in the seq function of the base package  
  #          data: A vector of your data points. 
  #        priors: A matrix of densities with in each column a density of a specific prior mapped on the paramater space that
  #                is equal to the parameter space that is supplied using the from, to, by statements. E.g. the parameter space
  #                runs from -10 to 10 in steps of 0.01 than your density of a standard normal distribution should be obtained
  #                using dnorm(x = seq(from = -10, to = 10, by = 0.01), mean = 0, sd = 1). The first column will thus describe this
  #                density using 2001 rows and all other columns should use the same density mapping to the parameter space.
  #      lb.bench: Lower boundary of the benchmark prior.
  #      ub.bench: Upper boundary of the benchmark prior.
  #        n.iter: The number of iterations that is used to obtain the posterior distribution of the data and the benchmark prior
  #                note that only half of these iterations will be used to obtain samples, the other half is used for adaptation and
  #                burnin.
  #
  # Returns: 
  #   
  #
  # Warning messages 
  # We will check if the parameter space defined matches the length of the densities of the priors
  if(length(seq(from = from,  to = to, by = by)) != nrow(priors) ){
    stop("The length of your defined parameter space does not match the length of the densities supplies in the priors input.")
  }
  
  # We will check if all distributions are propper and integrate to one
  for(npriors in 1:ncol(priors)){
    if(round(integrate.xy(x = seq(from = from,  to = to, by = by), fx = priors[, npriors]), 2) != 1){
      stop("One or more of your defined priors is not propper in the sense that the density does not integrate to one over
           the defined parameter space. You can use the integrate.xy function of the sfsmisc package to check this.")
    }
  } # end for loop
  
  # now checking if the benchmark prior integrates to one
  if(round(integrate.xy(x = seq(from = from,  to = to, by = by), 
                        fx = dunif(x = seq(from = from,  to = to, by = by), min=lb.bench, max=ub.bench) ), 2) != 1){
    stop("Your benchmark prior is not propper in the sense that the density does not integrate to one over
         the defined parameter space. You can use the integrate.xy function of the sfsmisc package to check this.")
  }
  
  # Create space to store the results 
  
  output.data <- matrix(NA, nrow=1, ncol=6)
  colnames(output.data) <- c("Sample mean", "Sample SE", "Posterior mean", "Posterior SD", "Abs. mean diff", "Abs. SD diff") 
  
  # Sample mean en SE of mean
    output.data[1,1] <- mean(data)
    output.data[1,2] <- sd(data)/sqrt(length(data))
    
    #-----------------------------------------------------------------
    
    # Posterior calculation with blavaan
    
    # Defining the uniform prior for blavaan
    prior <- paste("dunif(",lb.bench,",",ub.bench,")") 
    
    get.posterior <- function(data, prior){ 
      # Creating a matrix of the data with column name y by which the model specified below is generally applicable for different datasets.
      data <- as.matrix(data)
      colnames(data) <- c("y")
      
      #Defining the model, which is an intercept only model in case of the mean 
      model <- '#Intercept
      y ~ 1 
      
      #variance
      y ~~ prior("dgamma(1, 1)")*y
      '
      
      fit.model <- blavaan(model, data=data, n.chains = 2, burnin = n.iter/2, sample = n.iter, adapt=n.iter/2, dp=dpriors(nu=prior))
      return(fit.model)
    }
    
    #surpress blavaan output during DAC_Uniform run
    capture.output(output.blavaan <- get.posterior(data, prior))
    
    #Save posterior mean and posterior sd of the mean
    output.data[1,3] <- as.numeric(blavInspect(output.blavaan, what="postmean")[1])
    output.data[1,4] <- as.numeric(blavInspect(output.blavaan, what="se")$nu)
    
    #Save the absolute difference between sample and posterior estimates
    output.data[1,5] <- abs(output.data[1,3] - output.data[1,1])
    output.data[1,6] <- abs(output.data[1,4] - output.data[1,2])
    
    #------------------------------------------------------------------
    # Kullback-Leibler calculation 
    
    # Set range of x-axis 
    x.axis <- seq(from = from, to=to, by = by)
    
    # Specify posterior density 
    post <- dtruncnorm(x.axis, a=lb.bench, b=ub.bench, output[1,3], output[1,4])
    
    #density(hierdechainssamengevoegd, n = length( seq(from = from, to = to, by = by) ), from = from, to = to)
    
    
    # Specify benchmark density
    benchmark <- dunif(x.axis, lb.bench, ub.bench)
    
    #Creating a matrix of the posterior and benchprior for the use of the KLdiv function
    matrix1 <- cbind(post, benchmark)
    
    #Kullback-Leibler Divergence for posterior and benchmark prior
    KL.bench <- KLdiv(matrix1, method= c("continuous"), eps=10^-250)[1,2]
    
    
    # Priors 
    # Matrix with posterior density, and expert densities
    matrix2 <- cbind(post, priors)
    
    # Create space to store the KL expert values, and DAC scores 
    KL.experts <- matrix(NA, nrow=ncol(priors), ncol=2)
    colnames(KL.experts) <- c("KL expert", "DAC")
    
    # Get KL divergences for the priors 
    for(i in 1:ncol(priors)){
      KL.experts[i,1] <- KLdiv(cbind(matrix2[,1],matrix2[,i+1]),method= c("continuous"), eps=10^-250)[1,2]
      # get KL divergence between posterior (column 1 of matrix 2) for each expert (expert 1 = column 2 in matrx 2 etc.)
    }

    
    # DAC scores for each expert 
    for(i in 1:ncol(priors)){
      KL.experts[i,2] <- KL.experts[i,1]/KL.bench
    }
    
    out <- list(output.data = output.data, KL.DAC = KL.experts)
    
    print(KL.experts)
    
    # Return output
    return(out)
  
}

