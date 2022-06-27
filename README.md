# CPUEspatial
A R-package that I developed during my PhD. It is geostatistical model that was developed to explore preferential sampling from fisheries when analysing fishery-dependent catch and effort data, also known as a Catch Per Unit Effort (CPUE) analysis.
It was built to explore preferential sampling, in an attempt to improve CPUE methods by modeling the spatial distribution of both catch rates and fishing locations using shared spatial covariates and Gaussian Fields. The package exploits recent advancements in Software TMB (Kristensen et al., [2016](https://www.jstatsoft.org/article/view/v070i05)) and Random Markov Gaussian Field theory (Lindgren et al., [2011](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2011.00777.x)). The package is very similar to [geostatistical_delta-GLMM](https://github.com/nwfsc-assess/geostatistical_delta-GLMM) and was a source of inspiration. It differs in a few key ways firstly CPUEspatial is only univariate, and secondly CPUEspatial allows the preferential sampling model. The structure of the Package was created by TMBtools [highly recommended](https://github.com/mlysy/TMBtools) for those who want to create packages that depend on TMB. 

## The Model
Focuses on a single species (univariate response variable), compared with other packages such as VAST that extend to multispecies cases and the delta method. This is a GLMM model with forced temporal model structure, and optional spatial structure. The idea of the package was to code up a GLMM which had all the helper summary plotting functions that are easy to use for users. The systematic component allows for linear terms (classis GLM) and spline-based smoothers (for continuous covariates). As well as time-invariant spatial Gaussian Fields (omega) and a time-varying Gaussian Fields (epsilon). The structure of the systematic component is similar to that of VAST where spatial, temporal and catchability covariates are seperated out in the linear predictor term. This was done to investigate accounting for preferential sampling in a CPUE analysis context. The model matrix for catchability covariates is consistent with a normal GLM model matrix with intercept and contrasts (for factors) or slope for continuous variables. However, the time coefficients are constrained to sum = 0 and spatial factor coeffecients are constrained to sum to zero. This is for identifiability, users should see the documentation in the folder [here](https://github.com/Craig44/CPUEspatial/tree/master/inst/examples).

## Features to test and validate
- Simulation (self test)
- Presence/Absence 
- Comprehensive documentation....
- Example case
- Diagnostics currently based on DHARMa R package [see here](https://github.com/florianhartig/DHARMa)
- Merge nearest neighbor and triangulation spatial approach currently in two separate source files that share 90% the same code
- Variable selection help - an issue I have found with these models are they can be slow to evaluate MLE estimates, so variable selection can be a tiresome approach if you try an stepwise approach, i.e. forward with some deviance acceptance criteria


During my PhD I identified the next most important model developments should be 
- Explore alternative point process likelihoods for fishing locations
- Allow point process to have different covariates and Gaussian Fields (GF) to the catch rate (more flexibility)
- Add R functions to conduct the predictive joint log-likelihood diagnostic

## Installation
System requirements include R version `> 3.5` so you will need a newish R version to install this R package. It also depend on in the R package INLA, which is not on CRAN. See [here](https://www.r-inla.org/download-install) for install instructions. Most other packages should be on CRAN.
```
devtools::install_github("Craig44/CPUEspatial")
```
