
#################################################################
#####   The following file includes code to recreate        ##### 
#####   all figures and tables in Chapter 7 Discussion      ##### 
#####   Figures are reproducible with the files in          ##### 
#####   the data folder.                                    #####
#################################################################

## Note, data was read in from each model to complete the plots in the 
## discussion. Further, the data for MCP comparisons are part of another project
## and can be found at the github ElaineYChu/aaba_2022 For full reproducibility,
## the copula models and MCP models can be requested. 

####################################BEGIN####################################### 

## set seed for reproducibility

set.seed(1991)

## load in all necessary packages

library(tidyverse)
library(magrittr)
library(condMVNorm)
library(cmdstanr)
library(posterior)

#####################
# Copula Imputation #
#####################

# Load in the correlation matrix from the full model
cop_corr <- as.matrix(read.csv("allvars_matrix.csv"))
rownames(cop_corr) <- colnames(cop_corr)

# Load in initial US dataset
dat <- read.csv("us_diss.csv")

# I want to use a subset of traits to predict femoral growth: 
# (FDL), FBDL, HPB, max_M1, man_C, TC_Oss, RDE_EF

fdl_pars <- read_cmdstan_csv(files = 
  c("allvars/mixed_gaussian_copula_FINAL-202302060556-1-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-2-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-3-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-4-94c541.csv"),
  variables = c("constant[1]", "exponent[1]", "offset[1]","noise1[1]",
  "noise2[1]"), sampler_diagnostics = "")

fdl_summary <- summarise_draws(fdl_pars$post_warmup_draws, 
                               .num_args = list(sigfig = 4, notation = "dec"))

fbdl_pars <- read_cmdstan_csv(files = 
  c("allvars/mixed_gaussian_copula_FINAL-202302060556-1-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-2-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-3-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-4-94c541.csv"), 
  variables = c("constant[8]", "exponent[8]", "offset[8]","noise1[8]",
  "noise2[8]"), sampler_diagnostics = "")

fbdl_summary <- summarise_draws(fbdl_pars$post_warmup_draws, 
                                .num_args = list(sigfig = 4, notation = "dec"))


hpb_pars <- read_cmdstan_csv(files = 
  c("allvars/mixed_gaussian_copula_FINAL-202302060556-1-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-2-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-3-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-4-94c541.csv"), 
  variables = c("constant[10]", "exponent[10]", "offset[10]","noise1[10]",
  "noise2[10]"), sampler_diagnostics = "")

hpb_summary <- summarise_draws(hpb_pars$post_warmup_draws, 
                              .num_args = list(sigfig = 4, notation = "dec"))

maxM1_pars <- read_cmdstan_csv(files = 
  c("allvars/mixed_gaussian_copula_FINAL-202302060556-1-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-2-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-3-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-4-94c541.csv"), 
  variables = c("constant_dentition[1]","u_dentition","threshold_dentition"), 
  sampler_diagnostics = "")

m1_summary <- summarise_draws(maxM1_pars$post_warmup_draws, 
                              .num_args = list(sigfig = 4, notation = "dec"))

manC_pars <- read_cmdstan_csv(files = 
  c("allvars/mixed_gaussian_copula_FINAL-202302060556-1-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-2-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-3-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-4-94c541.csv"), 
  variables = c("constant_dentition[14]"), sampler_diagnostics = "")

manC_summary <- summarise_draws(manC_pars$post_warmup_draws, 
                                .num_args = list(sigfig = 4, notation = "dec"))

TCoss_pars <- read_cmdstan_csv(files = 
  c("allvars/mixed_gaussian_copula_FINAL-202302060556-1-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-2-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-3-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-4-94c541.csv"), 
  variables = c("constant_tarsal", "u_tarsal", "threshold_tarsal"), 
  sampler_diagnostics = "")

TCoss_summary <- summarise_draws(TCoss_pars$post_warmup_draws, 
                                .num_args = list(sigfig = 4, notation = "dec"))

RDEef_pars <- read_cmdstan_csv(files = 
  c("allvars/mixed_gaussian_copula_FINAL-202302060556-1-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-2-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-3-94c541.csv",
  "allvars/mixed_gaussian_copula_FINAL-202302060556-4-94c541.csv"), 
  variables = c("constant_ef[13]", "u_ef", "threshold_ef"), 
  sampler_diagnostics = "")

RDEef_summary <- summarise_draws(RDEef_pars$post_warmup_draws, 
                                 .num_args = list(sigfig = 4, notation = "dec"))

## Prepare variables for analysis

# Continuous first
fdl_mean <- c()
fdl_sd <- c()
fdl_z <- c()
fbdl_mean <- c()
fbdl_sd <- c()
fbdl_z <- c()
hpb_mean <- c()
hpb_sd <- c()
hpb_z <- c()

for(i in 1:nrow(dat)){
  fdl_mean[i] <- unlist(fdl_summary[1,2]*dat$agey[i]^fdl_summary[2,2] + 
                          fdl_summary[3,2])
  fdl_sd[i] <- unlist(fdl_summary[4,2]*(1+fdl_summary[5,2]*dat$agey[i]))
  fdl_z[i] <- (dat$FDL_L[i] - fdl_mean[i]) / fdl_sd[i]
  
  fbdl_mean[i] <- unlist(fbdl_summary[1,2]*dat$agey[i]^fbdl_summary[2,2] + 
                           fbdl_summary[3,2])
  fbdl_sd[i] <- unlist(fbdl_summary[4,2]*(1+fbdl_summary[5,2]*dat$agey[i]))
  fbdl_z[i] <- (dat$FBDL_L[i] - fbdl_mean[i]) / fbdl_sd[i]
  
  hpb_mean[i] <- unlist(hpb_summary[1,2]*dat$agey[i]^hpb_summary[2,2] + 
                          hpb_summary[3,2])
  hpb_sd[i] <- unlist(hpb_summary[4,2]*(1+hpb_summary[5,2]*dat$agey[i]))
  hpb_z[i] <- (dat$HPB_L[i] - hpb_mean[i]) / hpb_sd[i]
}

# Ordinal
m1_u <- m1_summary[2:1317,]
m1_thresholds <- m1_summary %>% 
  filter(variable %in% c("threshold_dentition[1,1]", "threshold_dentition[1,2]",
                         "threshold_dentition[1,3]", "threshold_dentition[1,4]",
                         "threshold_dentition[1,5]", "threshold_dentition[1,6]",
                         "threshold_dentition[1,7]", "threshold_dentition[1,8]",
                         "threshold_dentition[1,9]", 
                         "threshold_dentition[1,10]", 
                         "threshold_dentition[1,11]"))

thresholds <- as.numeric(m1_thresholds$mean)
beta <- as.numeric(m1_summary$mean[1])
dat <- allsub2
colnames(dat)[1] <- "agey"

C_u <- m1_summary[17110:18425,]
C_thresholds <- m1_summary %>% filter(variable %in% 
  c("threshold_dentition[14,1]", "threshold_dentition[14,2]",
  "threshold_dentition[14,3]", "threshold_dentition[14,4]",
  "threshold_dentition[14,5]", "threshold_dentition[14,6]",
  "threshold_dentition[14,7]", "threshold_dentition[14,8]",
  "threshold_dentition[14,9]", "threshold_dentition[14,10]",
   "threshold_dentition[14,11]"))

rde_u <- as.numeric(RDEef_summary$mean[15794:17109])
rde_thresholds <- RDEef_summary %>% filter(variable %in% 
  c("threshold_ef[13,1]", "threshold_ef[13,2]",
  "threshold_ef[13,3]", "threshold_ef[13,4]",
  "threshold_ef[13,5]", "threshold_ef[13,6]"))



tc_oss_constant <- as.numeric(TCoss_summary$mean[1])
tc_oss_u <- as.numeric(TCoss_summary$mean[2:1317])
tc_oss_thresh <- as.numeric(TCoss_summary$mean[1318:1324])

m1 <- dat$max_M1_L
c <- dat$man_C_L
tc_oss <- dat$TC_Oss
rde <- dat$RDE_EF_L

m1_z <- c()
c_z <- c()
tc_oss_z <- c()
rde_z <- c()

for (i in 1:length(m1)){
  if(is.na(m1[i])){
    m1_z[i] <- NA
  } else if (m1[i] == 1){
    bound <- pnorm((m1_thresholds$mean[1] - dat$agey[i]*m1_summary$mean[1]))
    m1_z[i] <- qnorm((bound*m1_u$mean[i]))
  } else if(m1[i] == 12){
    bound <- pnorm((m1_thresholds$mean[11] - dat$agey[i]*m1_summary$mean[1]))
    m1_z[i] <- qnorm((bound + (1-bound)*m1_u$mean[i]))
  } else {
    ub <- pnorm((m1_thresholds$mean[m1[i]] - dat$agey[i]*m1_summary$mean[1]))
    lb <- pnorm((m1_thresholds$mean[m1[i] - 1] - 
                   dat$agey[i]*m1_summary$mean[1]))
    m1_z[i] <- qnorm((lb + (ub-lb)*m1_u$mean[i]))
  }
} 

for (i in 1:length(c)){
  if(is.na(c[i])){
    c_z[i] <- NA
  } else if (c[i] == 1){
    bound <- pnorm((C_thresholds$mean[1] - dat$agey[i]*manC_summary$mean[1]))
    c_z[i] <- qnorm((bound*C_u$mean[i]))
  } else if(c[i] == 12){
    bound <- pnorm((C_thresholds$mean[11] - dat$agey[i]*manC_summary$mean[1]))
    c_z[i] <- qnorm((bound + (1-bound)*C_u$mean[i]))
  } else {
    ub <- pnorm((C_thresholds$mean[c[i]] - dat$agey[i]*manC_summary$mean[1]))
    lb <- pnorm((C_thresholds$mean[c[i] - 1] - 
                   dat$agey[i]*manC_summary$mean[1]))
    c_z[i] <- qnorm((lb + (ub-lb)*C_u$mean[i]))
  }
} 

for (i in 1:length(tc_oss)){
  if(is.na(tc_oss[i])){
    tc_oss_z[i] <- NA
  } else if (tc_oss[i] == 1){
    bound <- pnorm((tc_oss_thresh[1] - dat$agey[i]*tc_oss_constant))
    tc_oss_z[i] <- qnorm((bound*tc_oss_u[i]))
  } else if(tc_oss[i] == 8){
    bound <- pnorm((tc_oss_thresh[7] - dat$agey[i]*tc_oss_constant))
    tc_oss_z[i] <- qnorm((bound + (1-bound)*tc_oss_u[i]))
  } else {
    ub <- pnorm((tc_oss_thresh[tc_oss[i]] - dat$agey[i]*tc_oss_constant))
    lb <- pnorm((tc_oss_thresh[tc_oss[i] - 1] - dat$agey[i]*tc_oss_constant))
    tc_oss_z[i] <- qnorm((lb + (ub-lb)*tc_oss_u[i]))
  }
} 

for (i in 1:length(rde)){
  if(is.na(rde[i])){
    rde_z[i] <- NA
  } else if (rde[i] == 1){
    bound <- pnorm((rde_thresholds$mean[1] - dat$agey[i]*RDEef_summary$mean[1]))
    rde_z[i] <- qnorm((bound*rde_u[i]))
  } else if(rde[i] == 7){
    bound <- pnorm((rde_thresholds$mean[6] - dat$agey[i]*RDEef_summary$mean[1]))
    rde_z[i] <- qnorm((bound + (1-bound)*rde_u[i]))
  } else {
    ub <- pnorm((rde_thresholds$mean[rde[i]] - 
                   dat$agey[i]*RDEef_summary$mean[1]))
    lb <- pnorm((rde_thresholds$mean[rde[i] - 1] - 
                   dat$agey[i]*RDEef_summary$mean[1]))
    rde_z[i] <- qnorm((lb + (ub-lb)*rde_u[i]))
  }
} 

####

all_sub <- cbind(dat$agey, fdl_z, fbdl_z, hpb_z, m1_z, c_z, tc_oss_z, rde_z, m1)

allsub2 <- as.data.frame(all_sub) %>% select(V1, fdl_z, tc_oss_z, m1_z, m1) %>% 
  na.omit()

cor_sub2 <- cop_corr[c("FDL", "max_M1","TC_Oss"),c("FDL", "max_M1","TC_Oss")]

c_m1 <- c()

given_vars <- allsub2 %>% select(1:2)
given_vars <- as.matrix(given_vars)

for(i in 1:nrow(allsub2)){
  c_m1[i] <- condMVN(mean = c(0,0,0), sigma = cor_sub2, dependent = 2,
                     given = c(1,3), X.given = given_vars[i,1:2])
}

c_m1 <- unlist(c_m1)

df <- as.data.frame(cbind(allsub2$m1, cats))
df2 <- gather(df, key = "group")
df2$agey <- c(dat$agey, dat$agey)
df3 <- as.data.frame(cbind(allsub2$V1, allsub2$m1, cats))

df4 <- df4 %>% group_by(as.integer(agey), group)
df4 <- rename(df4, age_group = 'as.integer(agey)')
df4$age_group <- as.factor(df4$age_group)

df2  %>% ggplot(aes(x = value, fill = group)) + geom_histogram(color = "white", 
  position = "dodge") + scale_fill_manual(name ="", values = c("black", 
  "tomato"), labels = c("Observed", "Predicted")) + 
  scale_x_continuous(breaks = seq(1,12,1)) + xlab("Maxillary M1 Score") + 
  theme_classic() + theme(legend.position = "bottom")


df2  %>% ggplot(aes(x = as.factor(as.integer(agey)),y = value, fill = group)) + 
  geom_boxplot() + scale_fill_manual(name ="", values = c("black", "tomato"), 
  labels = c("Observed", "Predicted")) + scale_y_continuous(name = "Score", 
  breaks = seq(1,13,1)) + xlab("Age [years]") +theme_classic() + 
  theme(legend.position = "bottom")

df3 <- df3 %>% group_by(as.integer(V1))
colnames(df3) [4] <- "agey"

df3 %>% ggplot() + geom_boxplot(aes(as.integer(agey), cats, group = agey), 
  position = position_jitter()) + geom_boxplot(aes(as.integer(agey), V2, 
  group = agey), col="tomato")

df4 <- df3 %>% select(-V1)
df4 <- df4 %>% group_by(agey) %>% gather(key = "group")

test <- cbind(dat$agey,fdl_z, fbdl_z, hpb_z, fdl_mean, fdl_sd)
cor_sub <- cop_corr[c("FDL","FBDL","HPB"),c("FDL","FBDL","HPB")]

condMVN(mean=c(0,0,0), sigma=cor_sub, dependent= 1, given=c(2,3),
        X.given=c(3.52,0.80))

test %<>% na.omit()

Cmeam <- c()

for(i in 1:nrow(test)){
  Cmeam[i] <- condMVN(mean = c(0,0,0), sigma = cor_sub, dependent = 1,
                      given = c(2,3), X.given = test[i,2:3])
}
cmean <- unlist(Cmeam)

fdl_unscaled <- c()
cmean_unscaled <- c()

for(n in 1:nrow(test)){
  fdl_unscaled[n] <- (test[n,2]*test[n,6]) + test[n,5]
  cmean_unscaled[n] <- (cmean[n]* test[n,6]) + test[n,5]
}

plot(test[,1], fdl_unscaled, pch=16, xlab = "Age", ylab = "FDL")
points(test[,1], cmean_unscaled, col = "red", pch=16, cex=0.5)
legend("topleft", legend = c("Observed", "Predicted"), col = c("black", "red"), 
       pch = 16)

growth_test <- as.data.frame(cbind(test[,1], fdl_unscaled, cmean_unscaled))

growth_test %>% ggplot() + geom_point(aes(V1, fdl_unscaled, col = "Observed")) + 
  geom_point(aes(V1, cmean_unscaled, col = "Predicted")) + labs(x="Age [years]",
  y="Femur Diaphyseal Length [mm]") + scale_color_manual(name = "", 
  values = c("black", "tomato")) + theme_bw() + theme(legend.position="bottom")


growth_test %>% ggplot() + geom_point(aes(fdl_unscaled, cmean_unscaled)) + 
  labs(x="Observed", y="Predicted") +
  geom_abline(aes(slope=1, intercept=0), color="tomato", linewidth=1) +
  theme_bw()



plot(fdl_unscaled, cmean_unscaled, xlab = "Observed", ylab = "Predicted", 
     pch=16)
abline(0,1,lwd=3, col="tomato")

#####################
# MCP Comparison    #
#####################

# Load in 6 variable MCP model & Copula Correlations
mcp_mod <- readRDS("cdep_model_US_sixvar.rds")
cop_corr <- as.matrix(read.csv("allvars_matrix.csv"))
rownames(cop_corr) <- colnames(cop_corr)


# Pull out 6 variable subset from the copula results
namvec <- c("HME_EF","TC_Oss","max_M1","man_I2","FDL","RDL")
cop_corr <- cop_corr[namvec,namvec]
lower_cop <- cop_corr[lower.tri(cop_corr)]

# Edit the MCP model with the copula correlations
mcp_mod$mod_spec$cdep_groups <- c(1,2,3,4,5,6)
mcp_mod$th_y <- c(mcp_mod$th_y[1:54],lower_cop)

# Load in test data and var info file to prep for posterior calculations
test <- read.csv("us_test_samp.csv")
test %<>% select(1:4, 16, 22, 37, 54, 61)

xcalc <- seq(0,25,by=0.01)
th_x <- readRDS("solutionx_US_sixvar.rds")

multivariate_batch_calc(data_dir="aaba_2022-main/aaba_2022-main/models",
                        analysis_name="US_sixvar",
                        test_samp=test,
                        demo_cols=1:3,
                        model_type="cdep",
                        ci_type="hdi",
                        th_x=th_x, xcalc=xcalc,
                        seed=2021, save_file=T)

multipred <- read.csv("cdep_model_US_sixvar_test_predictions.csv")

ggplot(multipred) + 
  geom_errorbar(aes(x=agey, ymin=lower95, ymax=upper95)) + 
  geom_point(aes(x=agey, y=xmode)) + 
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") + 
  labs(x="Known Age [years]", y="Predicted Age [years]") + 
  theme_bw()

multi_summary <- multipred %>% group_by(agey=as.integer(agey)) %>%
  summarize(point=mean(xmean, na.rm = T),lower95=mean(lower95, na.rm=T),
            upper95=mean(upper95, na.rm=T))

ggplot(multi_summary) +
  geom_point(aes(x=agey, y=point)) +
  geom_errorbar(aes(x=agey, ymin=lower95, ymax=upper95)) +
  geom_abline(aes(slope=1, intercept=0), linetype="dashed", color="black") + 
    labs(x="Age [years]", y="Average Point Estimate and Range")

# Accuracy
df_new <- multipred %>% mutate(accurate=ifelse(agey <= upper95 & 
                                          agey >= lower95,
                                        T,F))
n_samp <- df_new %>% drop_na() %>% nrow()  
accuracy <- length(which(df_new$accurate==T))/n_samp
out <- accuracy
out

# rmse
rmse <- RMSE(multipred$xmean, multipred$agey, na.rm=TRUE)
out <- rmse

# see 
see <- sd(multipred$xmean - multipred$agey, na.rm=TRUE)
out <- paste0("+/- ",round(see*1.96,2))

# TMLP
mcp_mod2 <- readRDS("cdep_model_US_sixvar_AABA.rds")
final2 <- data.frame(matrix(nrow=1, ncol=3))
names(final2) <- c("model","N","tmlp")
final2$model <- "mcp"


df <- test %>% select(3, 8, 9, 6, 7, 4, 5)

if(ncol(df)==2) {
  df <- na.omit(df)
} else {
  df <- df[-idx_all_na(df,2:ncol(df)),]
}

val_vec <- c()

for(j in 1:nrow(df)) {
  print(paste0(j,"/",nrow(df)))
  age <- round(df[j,"agey"],2)
  post <- calc_x_posterior(t(df[j,-1]), th_x, mcp_mod2, xcalc = xcalc, 
                           normalize = F)
  val <- post$density[round(post$x,2)==age]
  
  val_vec <- c(val_vec, log(val))
}

if(length(val_vec) != nrow(df)) {
  stop("val_vec is not equal to number of individuals")
} else {
  N <- final2[,"N"] <- length(val_vec)
}

final2[,"tmlp"] <- round(-sum(val_vec, na.rm=T) / N, 4)

cop_results <- read.csv("six_var_cop.csv")
mcp_results <- read.csv("US_sixvar_cdep_model_test_predictions.csv")

cop_results %<>% mutate(cop_resid = agey - xmean, 
                        cop_abs_resid = abs(agey - xmean))

all <- as.data.frame(cbind(cop_results$agey,mcp_results$resid, 
                           mcp_results$abs_resid, cop_results$cop_resid, 
                           cop_results$cop_abs_resid))

colnames(all) <- c("age", "mcp_bias", "mcp_accuracy", "cop_bias", 
                   "cop_accuracy")

all %>% ggplot() + geom_smooth(aes(x = age, y = as.numeric(mcp_bias)), se=F, 
  color = "black", method = "loess") + geom_smooth(aes(x = age, 
  y = as.numeric(cop_bias)), se=F, color = "tomato") + 
  geom_hline(yintercept=0, linetype="dashed") + theme_classic()

mcp_results %>% ggplot() + geom_smooth(aes(known_age, resid), se=F) + 
  geom_hline(yintercept=0, linetype="dashed") + theme_classic()
cop_results %>% ggplot() + geom_smooth(aes(agey, cop_resid), se=F) + 
  geom_hline(yintercept=0, linetype="dashed") + theme_classic()

plot(loess.smooth(mcp_results$known_age, mcp_results$resid), type = "l", 
     ylim = c(-1,6), lwd = 2, xlab = "Age[years]", ylab = "Bias")
lines(loess.smooth(cop_results$agey, cop_results$cop_resid), col = "tomato", 
      lwd=2)
abline(h = 0, lty = "dashed", lwd = 2)
legend("topleft",legend = c("MCP", "Copula"), lty = c(1,1), lwd = c(2,2), 
       col = c("black", "tomato"), bty = "n")

plot(loess.smooth(mcp_results$known_age, mcp_results$abs_resid), type = "l", 
     ylim = c(0,6), lwd = 2, xlab = "Age[years]", ylab = "Accuracy")
lines(loess.smooth(cop_results$agey, cop_results$cop_abs_resid), col = "tomato",
      lwd=2)
abline(h = 0, lty = "dashed", lwd = 2)
legend("topleft",legend = c("MCP", "Copula"), lty = c(1,1), lwd = c(2,2), 
       col = c("black", "tomato"), bty = "n")

##############################END##############################################
