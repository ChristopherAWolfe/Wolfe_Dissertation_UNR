
#################################################################
#####   The following file includes the final code to       ##### 
#####   fit each stan model in cmdstanr. Note, this         #####
#####   file is to be after the model prep code to          #####
#####   prepare the stan data.                              #####
#################################################################

####################################BEGIN####################################### 


library(cmdstanr)

stanmodel <- cmdstan_model(stan_file = "mixed_gaussian_copula_FINAL2.stan")

fit <-stanmodel$sample(data = stan_dat,
                 seed = 1991,
                 refresh = 500,
                 init = 0,
                 chains = 4,
                 parallel_chains = 4,
                 iter_warmup = 3000,
                 iter_sampling = 6000,
                 max_treedepth = 12,
                 adapt_delta = 0.99)

##############################END##############################################