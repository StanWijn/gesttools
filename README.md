# gesttools: General Purpose G estimation in R
Here we document gesttools, a series of general purpose, user friendly functions with which to perform
g-estimation of structural nested mean models (SNMMs), with the intent to provide g-estimation techniques for datasets 
with a wide variety of exposure and outcome types in R.

The suite of functions provided in this repository are designed for use with longitudinal, long format datasets, with data organised into
a set of discrete time periods t=1,...,T. The datasets exposure must be time varying, but may be continuous, binary or even 
categorical in nature. The datasets outcome may be measured at the end of the study, or repeatedly measured at every time period, 
and may be a continous, binary or count outcome. 

In addition gesttools provides users with a choice of four types of SNMMs to fit, which allows for a variety of causal effects 
to be estimated from the data. In particular gesttools permits users to fit time-specific causal effects, as well as allowing causal effects
to be modified by some other covariate.

The functions also support censoring weights for missing data and confidence intervals are available via a bootstrapping procedure.

In order to run gesttools, first run "Loadpackages.R" to obtain the required dependent packages, then run the 5 major functions
that perform g-estimation, "gest.R", "gestmult.R", "gest.cat", "gestmult.cat", and the bootstrap function "gest.boot". 

A user guide for using gesttools can be found in this repository under "gesttoolsimplementation.pdf" which uses R code found in 
"SimulatedExamples.R"

An associated paper is in development to provide a comprehensive introduction to gesttools including an explanation of the SNMM types,
the g-estimation algorithm, instructions to set up the users dataset, explanation of each function, and a turotial to perform g-estimation,

For any queries please contact Dr Daniel Tompsett at d.tompsett@ucl.ac.uk.

