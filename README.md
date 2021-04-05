# CPUEspatial
An R-package that I started developing during my PhD. It applies geostatistical models mainly for analysising Commercial Fishing data, also known as Catch Per Unit Effort (CPUE) analysis.
It was built to use extend current practices and exploit recent advancements in Software TMB (add reference) and geostatistical techniques (link to INLA and lindgren et. al). The package is very similar to https://github.com/nwfsc-assess/geostatistical_delta-GLMM, the reason I have gone throught all this is work. Is sometimes you want to add stuff or do silly things that packages owners don't want. The structure of the Package was created by TMBtools [highly recommended](https://github.com/mlysy/TMBtools) for those who want to create packages that depend on TMB. 

## The Model
Focuses on a single species (univariate response variable), compared with other packages such as VAST that extend to multispecies cases. Bascially this is a GLMM model with forced temporal structure, and optional spatial structure. The idea of the package was to code up a GLMM which had all the helper summary plotting functions that are easy to use for users. The systematic component allows for linear terms (classis GLM) and Penalised cubic splines (for continuous covariates).
As well as time-invariant spatial GMRF (omega) and a time-varying spatial GMRF (epsilon). The structure of the systematic component is similar to that of VAST where spatially varying covariates and catchability covariates are seperated. This was done to investigate accounting for preferential sampling in a CPUE analysis context. The model matrix for catchability covariates is consistent with a normal GLM model matrix with intercept and contrasts (for factors) or slope for continuous variables. But the time coeffecients are constrained to sum = 0 and spatial factor coeffecients are constrained to sum to zero. This is for identifiability, users should see the documentation in the folder [here](https://github.com/Craig44/CPUEspatial/tree/master/inst/examples).

## Features to test and validate
- Simulation (self test)
- Precense/Absence 
- Comprehensive documentation....
- Example case
- Diagnostics currently based on DHARMa R package [see here](https://github.com/florianhartig/DHARMa)
- Merge nearest neighbour and triangulation spatial approachs currently in two seperate source files that share 90% the same code
- Variable selection help - an issue I have found with these models are they can be slow to evalueate MLE estimates, so variable selection can be a tiresome approach if you try an stepwise approach, i.e. forward with some deviance acceptance criteria
## Installation
System requirements include R version `> 3.5` so you will need a newish R version to install this R package. It also depend on in the R package INLA, which is not on CRAN. See [here](https://www.r-inla.org/download-install) for install instructions. Most other packages should be on CRAN.
```
devtools::install_github("Craig44/CPUEspatial")
```
