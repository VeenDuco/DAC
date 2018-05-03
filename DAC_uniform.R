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
