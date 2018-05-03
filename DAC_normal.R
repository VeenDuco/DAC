# Copyright statement:
# The GNU GENERAL PUBLIC LICENSE v3.0 applies.

# Author comment
# no comment at this time

# File description comment
# This file is meant to calculate the Data Agreement Criterion (DAC) (Bousquet, 2008)
# for a 2 parameter (mean / sd) model. The DAC can be calculated for multiple priors 
# simultaniously such that we can rank and compare these priors. 
# Priors can be elicited experts' beliefs such that a ranking of experts can be obtained.
# In this specific function we calculate the DAC whilst using a Normal benchmark prior. 

# loading required packages
library(blavaan)
library(flexmix)
library(sfsmisc)

# Function defenitions
# Function to calculate the DAC for a 2 parameter model mean / sd. 
DAC.normal <- function (from, to, by, data, priors, mean.bench, sd.bench, n.itter = 10000) {
  #
  # Args:
  #          from: Lower bound of the parameter space that is to be evaluated, as in the seq function of the base package
  #            to: Upper bound of the parameter space that is to be evaluated, as in the seq function of the base package  
  #            by: Step size by which the defined parameter space is mapped out, as in the seq function of the base package  
  #          data: A vector of your data points. 
  #        priors: A matrix of densities with in each column a density of a specific prior mapped on the paramater space that
  #                is equal to the parameter space that is supplied using the from, to, by statements. E.g. the parameter space
  #                runs from -10 to 10 in steps of 0.01 than your density of a standard normal distribution shoudl be obtained
  #                using dnorm(x = seq(from = -10, to = 10, by = 0.01), mean = 0, sd = 1). The first column will thus describe this
  #                density using 2001 rows and all other columns should use the same density mapping to the parameter space.
  #    mean.bench: Mean of the benchmark prior.
  #      sd.bench: sd of the benchmark prior.
  #       n.itter: The number of itterations that is used to obtain the posterior distribution of the data and the benchmark prior
  #                note that only half of these itterations will be used to obtain samples, the other half is used for adaptation and
  #                burnin.
  #
  # Returns: 
  #   
  #
  # Error handling:
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
                        fx = dnorm(x = seq(from = from,  to = to, by = by), mean = mean.bench, sd = sd.bench) ), 2) != 1){
    stop("Your benchmark prior is not propper in the sense that the density does not integrate to one over
         the defined parameter space. You can use the integrate.xy function of the sfsmisc package to check this.")
  }
  
  
} # end function