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

DAC.uniform <- function(from, to, by, data, priors, lb.bench, ub.bench, n.iter){

  # Create space to store the results 
  
  output.data <- matrix(NA, nrow=1, ncol=6)
  colnames(output) <- c("Sample mean", "Sample SE", "Posterior mean", "Posterior SD", "Abs. mean diff", "Abs. SD diff") 
  
  # Sample mean en SE
    output.data[1,1] <- mean(data)
    output[1,2] <- sd(data)/sqrt(length(data))
    
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
    
    # Set range of x-axis to largest uniform prior used, with specified length.out
    x.axis <- seq(from = from, to=to, by = by)
    
    # Specify posterior density 
    post <- dtruncnorm(x.axis, a=lb.bench, b=ub.bench, output[1,3], output[1,4])
    
    # Specify benchmark density
    benchmark <- dunif(x.axis, lb.bench, ub.bench)
    
    #Creating a matrix of the posterior and benchprior for the use of the KL D function
    matrix1 <- cbind(post, benchmark)
    
    #Kullback-Leibler Divergence for posterior and benchmark prior
    KL.bench <- KLdiv(matrix1, method= c("continuous"), eps=10^-250)[1,2]
    
    # Priors 
    
    # Matrix with priors 
    
    # Create space to store the KL expert values 
    
    # Matrix with posterior density, and expert densities
    matrix2 <- cbind(post, priors)
    
    KL.experts <- matrix(NA,nrow=ncol(priors), ncol=2)
    colnames(KL.experts) <- c("KL expert", "DAC")
    
    
    for(i in 1:ncol(priors)){
      KL.experts[i,1] <- KLdiv(cbind(matrix2[,1],matrix2[,i+1]),method= c("continuous"), eps=10^-250)[1,2]
      # get KL divergence between posterior (column 1 of matrix 2) for each expert (expert 1 = column 2 in matrx 2 etc.)
    }

    
    # DAC scores for each expert 
    for(i in 1:ncol(priors)){
      KL.experts[i,2] <- KL.experts[i,1]/KL.bench
    }
    
    out <- list(output.data = output.data, KL.experts = KL.experts)
    
    print(list(Means = Means, Variances = Variances, Ranges = ranges))
    
    # Return output
    return(out)
  
}