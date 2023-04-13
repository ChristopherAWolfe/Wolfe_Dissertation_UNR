
#################################################################
#####   The following file includes code to recreate        ##### 
#####   all figures and tables in Chapter 4 Materials       ##### 
#####   and Chapter 5 Methodology. Figures are reproducible ##### 
#####   with the files in the data folder.                  #####
#################################################################

####################################BEGIN#######################################        

## Set seed for reproducibility
set.seed(1991)
      
## Load in all necessary packages
library(tidyverse)
library(ggpubr)
library(ggExtra)
library(magrittr)
library(MASS)

## Load in all necessary files (located in the data folder)
observed_data <- read.csv("us_diss.csv") # full dataset 1,316 x 54
ppd_parameters <- read.csv("ppd_parameters.csv") # parameters for Fig. 5.3-5.6 
cor_matrix_FULL <- as.matrix(read.csv("allvars_matrix.csv")) # full corr matrix
rownames(cor_matrix_FULL) <- colnames(cor_matrix_FULL) # tidy up row names


####################
##   FIGURE 4.1   ##
###################

# Age and biological sex distribution of the study sample.

observed_data %>% dplyr::select(agey, SEX) %>% na.omit() %>% ggplot() +
geom_histogram(aes(agey, fill = SEX), position = "dodge") +
scale_x_continuous(name = "Age [years]", breaks = seq(0,22,1)) +
scale_y_continuous(name = "Number of Individuals", breaks = seq(0,95,5)) +
labs(fill = "Biological Sex") + scale_fill_manual(labels = c("Female (n = 533)",
    "Male (n = 783)"), values = c("black","grey")) + theme_classic() +  
  theme(legend.position = c(0.91,0.85), legend.title.align = 0.5, 
        legend.box.background = element_rect(fill='transparent'), 
        legend.background = element_rect(fill='transparent'))


####################
##   TABLE 4.2   ##
###################

# Count of individuals (pooled sex) across each developmental stage.

observed_data %>% dplyr::select(agey) %>% mutate(stage = ifelse(agey < 3, 
  "infancy", ifelse(agey >= 3.0 & agey <= 6.99, "childhood","juv_adol"))) %>% 
  group_by(stage) %>% summarise(count = n())
  

####################
##   FIGURE 5.1   ##
###################

# Demonstration of copula transformations.

## Implicit Gaussian Copula = MVNorm(0, CorMat)

randnums <- mvrnorm(n = 10000, mu = rep(0,54), Sigma = cor_matrix_FULL, 
                          empirical = T)

## Probability integral transform (uniformize the marginals)

p_out <- pnorm(randnums)

## Marginals are assumed Std Norm(0,1) - need to do inverse of the above

### Note technically this function will give similar results to the initial
### mvrnorm() call. I show here to demonstrate the full process if the marginal
### distributions were not std_normal - i.e., N(3,1). 

q_out <- qnorm(p_out)

## Observed Data FDL x RDL (top)

orig <- observed_data %>% ggplot() + geom_point(aes(FDL_L, RDL_L)) + 
  xlab("Femur Length [mm]") + ylab("Radius Length [mm]") + theme_classic()

orig <- ggMarginal(orig, type = "histogram") 

## Uniformized data (middle)

unif <- as.data.frame(p_out) %>% ggplot() + geom_point(aes(FDL, RDL)) + 
  xlab("Femur Length") + ylab("Radius Length") + theme_classic()

unif <- ggMarginal(unif, type = "histogram") 

## Inverse probability transform - std_normal()
cop <- as.data.frame(q_out) %>% ggplot() + geom_point(aes(FDL, RDL)) + 
  xlab("Femur Length") + ylab("Radius Length") + theme_classic()

cop <- ggMarginal(cop, type = "histogram") 

## Full plot

ggarrange(orig, unif, cop, nrow=3, ncol=1)

#####################
##   FIGURE 5.2   ##
####################

# Trace plots demonstrating convergence

## Note, the following code is commented out because the plots were directly 
## created from the samples - a large (>10gb) file not included. This is to
## demonstrate the steps taken to complete the plot. 

#library(bayesplot)
#color_scheme_set("blue")

#mcmc_trace(allvars_cors$post_warmup_draws, pars = c(...), facet_args = list())


#####################
##   FIGURE 5.3   ##
####################

# PPD check of FDL and RDL.

## Empty containers for calculations

fdl_mean <- NULL
fdl_sd <- NULL
rdl_mean <- NULL
rdl_sd <- NULL
fdl_z <- NULL
rdl_z <- NULL

# Loop through each values of x (age) computing the means and sds and then 
# standardizing each variable - this is the same step in the Stan code in the
# normal_marginal() function

for(i in 1:nrow(observed_data)){
  rdl_mean[i] <-  ppd_parameters$posterior.mean[6] * observed_data$agey[i] ^ 
    ppd_parameters$posterior.mean[7] + ppd_parameters$posterior.mean[8]
  rdl_sd[i] <- ppd_parameters$posterior.mean[9] * 
    (1 + ppd_parameters$posterior.mean[10] * observed_data$agey[i])
  fdl_mean[i] <- ppd_parameters$posterior.mean[1] * observed_data$agey[i] ^ 
    ppd_parameters$posterior.mean[2] + ppd_parameters$posterior.mean[3]
  fdl_sd[i] <- ppd_parameters$posterior.mean[4] * 
    (1 + ppd_parameters$posterior.mean[5] * observed_data$agey[i])
  fdl_z[i] <- (observed_data$FDL_L[i] - (fdl_mean[i])) / (fdl_sd[i])
  rdl_z[i] <- (observed_data$RDL_L[i] - (rdl_mean[i])) / (rdl_sd[i])
}

# put standardized values into same matrix for ease

z <- cbind(fdl_z, rdl_z)
fdl_ppd <- as.data.frame(cbind(fdl_mean, observed_data$FDL_L))

## Figure 5.3, top

fdl_ppd %>% na.omit() %>% ggplot() +geom_density(aes(fdl_mean, 
  fill = "Predicted"), alpha = 0.25) + geom_density(aes(V2, fill = "Observed"), 
  alpha = 0.25) + xlab("Femur Diaphyseal Length [mm]") +ylab("Density") + 
  scale_fill_manual(name ="", breaks = c("Predicted", "Observed"), 
  values = c("steelblue", "tomato")) + theme_classic() + 
  theme(legend.position = c(0.8, 0.8))

# Figure 5.3, bottom

ggplot() + geom_point(data = as.data.frame(q_out), aes(FDL, RDL, 
  color = "Predicted")) + geom_point(data = as.data.frame(z), aes(fdl_z, rdl_z, 
  color = "Observed")) + xlab("Femur Diaphyseal Length") + 
  ylab("Radius Diaphyseal Length") + scale_color_manual(name ="", 
  breaks = c("Predicted", "Observed"), values = c("steelblue", "tomato")) + 
  theme_classic() + theme(legend.position = c(0.2, 0.8))


#####################
##   FIGURE 5.4   ##
####################

# PPD of M1.

## Empty container
z_m1 <- NULL

## Distill parameters into vectors for ease
thresholds <- ppd_parameters$posterior.mean[12:22]
beta <- ppd_parameters$posterior.mean[11]
u <- ppd_parameters$posterior.mean[23:1338]

## Here I complete the same procedure from the probit_marginal function 
## that transforms the ordinal data to the latent scale --> z_m1

for(i in 1:nrow(observed_data)){
  if(is.na(observed_data$max_M1_L[i])){
    z_m1[i] <- NA
  }else if (observed_data$max_M1_L[i] == 1){
    bound <- pnorm((thresholds[1] - observed_data$agey[i] * beta))
    z_m1[i] <- qnorm((bound * u[i]))
  } else if(observed_data$max_M1_L[i] == 12){
    bound <- pnorm((thresholds[11] - observed_data$agey[i] * beta))
    z_m1[i] <- qnorm((bound + (1 - bound) * u[i]))
  } else{
    ub <- pnorm((thresholds[observed_data$max_M1_L[i]] - 
                   observed_data$agey[i] * beta))
    lb <- pnorm((thresholds[observed_data$max_M1_L[i] - 1] - 
                   observed_data$agey[i] * beta))
    z_m1[i] <- qnorm((lb + (ub - lb) * u[i]))
  }
}

## combining standardized / latent values for ease
z_both <- cbind(fdl_z, z_m1)


## Figure 5.4, top
ggplot() + geom_density(data = as.data.frame(q_out), aes(max_M1, 
  fill = "Predicted"), alpha = 0.25) + geom_density(data = as.data.frame(z_m1), 
  aes(z_m1, fill = "Observed"), alpha = 0.25) + 
  xlab("Maxillary M1 Latent Value") + ylab("Density") + 
  scale_fill_manual(name ="", breaks = c("Predicted", "Observed"), 
  values = c("steelblue", "tomato")) + theme_classic() + 
  theme(legend.position = c(0.8, 0.8))

## Figure 5.4, bottom
ggplot() + geom_point(data = as.data.frame(q_out), aes(FDL, max_M1, 
  color = "Predicted")) + geom_point(data = as.data.frame(z_both), 
  aes(fdl_z, z_m1, color = "Observed")) + xlab("Femur Diaphyseal Length") + 
  ylab("Maxillary M1 Latent Value") + scale_color_manual(name ="", 
  breaks = c("Predicted", "Observed"), values = c("steelblue", "tomato")) + 
  theme_classic() + theme(legend.position = c(0.15, 0.9))


#####################
##   FIGURE 5.5   ##
####################

# PPD of M1 on discrete scale

## Here I turn the simulated "pseudo-observations" from q_out into discrete 
## values as a means to complete a posterior check on the discrete scale

## Cumulative probabilities across all 12 dentition categories
post.p1 <- pnorm(thresholds[1] - observed_data$agey * beta)
post.p2 <- pnorm(thresholds[2] - observed_data$agey * beta) - post.p1
post.p3 <- pnorm(thresholds[3] - observed_data$agey * beta) - post.p1 - post.p2
post.p4 <- pnorm(thresholds[4] - observed_data$agey * beta) - post.p1 - 
        post.p2 - post.p3
post.p5 <- pnorm(thresholds[5] - observed_data$agey * beta) - post.p1 - 
        post.p2 - post.p3 - post.p4
post.p6 <- pnorm(thresholds[6] - observed_data$agey * beta) - post.p1 - 
        post.p2 - post.p3 - post.p4 - post.p5
post.p7 <- pnorm(thresholds[7] - observed_data$agey * beta) - post.p1 - 
        post.p2 - post.p3 - post.p4 - post.p5 - post.p6
post.p8 <- pnorm(thresholds[8] - observed_data$agey * beta) - post.p1 - 
        post.p2 - post.p3 - post.p4 - post.p5 - post.p6 - post.p7
post.p9 <- pnorm(thresholds[9] - observed_data$agey * beta) - post.p1 - 
        post.p2 - post.p3 - post.p4 - post.p5 - post.p6 - post.p7 - post.p8
post.p10 <- pnorm(thresholds[10] - observed_data$agey * beta) - post.p1 - 
        post.p2 - post.p3 - post.p4 - post.p5 - post.p6 - post.p7 - post.p8 - 
        post.p9
post.p11 <- pnorm(thresholds[11] - observed_data$agey * beta) - post.p1 - 
        post.p2 - post.p3 - post.p4 - post.p5 - post.p6 - post.p7 - post.p8 - 
        post.p9 - post.p10
post.p12 <- 1 - post.p1 - post.p2 - post.p3 - post.p4 - post.p5 - post.p6 - 
        post.p7 - post.p8 - post.p9 - post.p10 - post.p11

iter <- 5000

l <- length(observed_data$agey)

## put into same array
post.column.prob <- array(c(post.p1, post.p2, post.p3, post.p4, post.p5, 
    post.p6, post.p7, post.p8, post.p9, post.p10, post.p11, post.p12), c(l, 12))

## empty containers
simulated.category=array(0,c(length(observed_data$agey),12))

cats <- NULL

# Draw from the multinomial distribution based based on cumulative probs above
for(i in 1:l){
  for(m in 1:iter){
    simulated.category[i,] <- rmultinom(1, 1, post.column.prob[i,])
    ColNum = which(simulated.category[i,] == 1)
    cats[i] <- ColNum
  }
}

## combined for ease
m1_ppd <- cbind(cats, observed_data$max_M1_L)


## Figure 5.5
ggplot(as.data.frame(m1_ppd)) + geom_histogram(aes(cats, fill = "Predicted"), 
  alpha = 0.45) + geom_histogram(aes(V2, fill = "Observed"), alpha = 0.65) + 
  scale_fill_manual(name ="", breaks = c("Predicted", "Observed"), 
  values = c("steelblue", "tomato")) + 
  scale_x_continuous(breaks = seq(1,12,1)) + xlab("Maxillary M1 Score") + 
  theme_classic() + theme(legend.position = c(0.2, 0.8))

#####################################END########################################
