# Mixed-Discrete-Continuous-Gaussian-Copula

This github repository provides the data and code to reproduce the results from my doctoral dissertation entitled: *Beyond a Company of Soliders: Exploring Phenotypic Integration across the Multivariate Human Growth and Development Phenotype* submitted to the University of Nevada, Reno in May 2023. 

This dissertation fits a Mixed Discrete-Continuous Gaussian copula in Stan to address the underlying dependency structure of human growth traits. The **Data** folder contains the dataset used to fit the models, the mean posterior correlation matrices extracted from the individual Stan models, and additional files put together to complete the results. The **Code** folder contains both the Stan model and downstream code utilzed in all analyses separated by chapter in the dissertation. Note, the code for the *discussion* section needs additional data directly from individual Stan files that can be recreated with the Stan code and "model_fit" / "model_prep" scripts.

Note, the copula code is modified from two sources
1. https://spinkney.github.io/helpful_stan_functions/group__normal.html#ga53b0ac15ddf226fa86c96800014d89fc. 
2. https://github.com/stan-dev/math/blob/develop/stan/math/prim/prob/ordered_probit_lpmf.hpp. 
