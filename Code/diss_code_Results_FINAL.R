
#################################################################
#####   The following file includes code to recreate        ##### 
#####   all figures and tables in Chapter 6 Results         ##### 
#####   Figures are reproducible with the files in          ##### 
#####   the data folder.                                     #####
#################################################################

####################################BEGIN####################################### 

## set seed for reproducibility

set.seed(1991)

## load in all necessary packages

library(tidyverse)
library(reshape2)
library(ggExtra)
library(GGally)
library(ggsci)
library(ggpubr)

## define necessary functions and colors

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

col <- c(rep("steelblue", 18), rep("tomato", 16), rep("goldenrod", 20))

## Load in all necessary data files - this includes all posterior correlation 
## matrices and data necessary for the tables. I also create the pseudo-
## observations necessary for the copula plots. I also subset the upper triangle
## for each plot.

### Model 1: Full Model

cor_mat_allvars <- as.matrix(read.csv("allvars_matrix.csv"))
rownames(cor_mat_allvars) <- colnames(cor_mat_allvars)

randnums_all <- MASS::mvrnorm(n = 10000, mu = rep(0,54), 
                              Sigma = cor_mat_allvars, empirical = T)
p_out_all <- pnorm(randnums_all)
q_out_all <- as.data.frame(qnorm(p_out_all))

### Model 2: Infancy Model

cor_mat_infancy <- as.matrix(read.csv("infancy_matrix.csv"))
rownames(cor_mat_infancy) <- colnames(cor_mat_infancy)

randnums_infancy <- MASS::mvrnorm(n = 10000, mu = rep(0,54), 
                                  Sigma = cor_mat_infancy, empirical = T)
p_out_infancy <- pnorm(randnums_infancy)
q_out_infancy <- as.data.frame(qnorm(p_out_infancy))

### Model 3: Childhood Model

cor_mat_child <- as.matrix(read.csv("child_matrix.csv"))
rownames(cor_mat_child) <- colnames(cor_mat_child)

randnums_child<- MASS::mvrnorm(n = 10000, mu = rep(0,54), 
                                  Sigma = cor_mat_child, empirical = T)
p_out_child <- pnorm(randnums_child)
q_out_child <- as.data.frame(qnorm(p_out_child))

### Model 4: Juvenile and Adolescence Model

cor_mat_old <- as.matrix(read.csv("old_matrix.csv"))
rownames(cor_mat_old) <- colnames(cor_mat_old)

randnums_old <- MASS::mvrnorm(n = 10000, mu = rep(0,54), 
                                  Sigma = cor_mat_old, empirical = T)
p_out_old <- pnorm(randnums_old)
q_out_old <- as.data.frame(qnorm(p_out_old))

### Model 5: Biological Male Model

cor_mat_male <- as.matrix(read.csv("male_matrix.csv"))
rownames(cor_mat_male) <- colnames(cor_mat_male)

randnums_male <- MASS::mvrnorm(n = 10000, mu = rep(0,54), 
                                  Sigma = cor_mat_male, empirical = T)
p_out_male <- pnorm(randnums_male)
q_out_male <- as.data.frame(qnorm(p_out_male))

### Model 6: Biological Female Model

cor_mat_female <- as.matrix(read.csv("female_matrix.csv"))
rownames(cor_mat_female) <- colnames(cor_mat_female)

randnums_female <- MASS::mvrnorm(n = 10000, mu = rep(0,54), 
                               Sigma = cor_mat_female, empirical = T)
p_out_female <- pnorm(randnums_female)
q_out_female <- as.data.frame(qnorm(p_out_female))

### Other Data
parcoords <- read.csv("par_coords_data.csv")
names <- c("[DD,DD]","[DD,EF]","[DD,LBB]","[DD,LBL]","[DD,OSS]","[EF,EF]",
           "[EF,LBB]","[EF,LBL]","[EF,OSS]","[LBB,LBB]","[LBB,LBL]","[LBB,OSS]",
           "[LBL,LBL]","[LBL,OSS]","[OSS,OSS]")
colnames(parcoords) <- c("Model",names)

sex_comp <- read.csv("sex_long.csv")
lh_comp <- read.csv("lh_long.csv")

## Make a vector of column names to be used for later analyses

lbl <- colnames(cor_mat_allvars[,1:18])
dent <- colnames(cor_mat_allvars[,19:34])
ef <- colnames(cor_mat_allvars[,35:50])
pelvis <- colnames(cor_mat_allvars[,51:52])
oss <- colnames(cor_mat_allvars[,53:54])


####################
##   TABLE 6.1   ##
###################

# Spearman's Rho and associated permuted p-values comparing corr matrices

## Extract the upper triangle of each matrix
male_upper <- melt(get_upper_tri(cormat = cor_mat_male), na.rm = TRUE)
female_upper <- melt(get_upper_tri(cormat = cor_mat_female), na.rm = TRUE)
infancy_upper <- melt(get_upper_tri(cormat = cor_mat_infancy), na.rm = TRUE)
child_upper <- melt(get_upper_tri(cormat = cor_mat_child), na.rm = TRUE)
old_upper <- melt(get_upper_tri(cormat = cor_mat_old), na.rm = TRUE)

## Combine for correlation test
x <- as.data.frame(cbind(infancy_upper$value, child_upper$value, 
                         old_upper$value, male_upper$value, female_upper$value))

colnames(x) <- c("infancy", "childhood", "old", "male", "female")

## Spearman's correlation test
psych::corr.test(x = x, method = "spearman")


## Non-parametric permutation test for p-value


### Male x Female

rhos <- c()
n_iter <- 5000
true_rho <- cor(male_upper$value,female_upper$value, method = "spearman")

m_ids <- colnames(cor_mat_male)
m2_v <- female_upper$value
for(iter in 1:n_iter){
  m_ids <- sample(m_ids)
  r <- cor(melt(get_upper_tri(cor_mat_male[m_ids,m_ids]), na.rm = TRUE)$value, 
           m2_v, method = "spearman")
  rhos <- c(rhos,r)
}

perm_p_MF <- (sum(abs(true_rho) <= abs(rhos)) + 1)/(n_iter + 1) #two tailed test
perm_p_MF

### Infancy x Old

rhos <- c()
n_iter <- 5000
true_rho <- cor(infancy_upper$value,old_upper$value, method = "spearman")

# matrix permutation
m_ids <- colnames(cor_mat_infancy)
m2_v <- old_upper$value
for(iter in 1:n_iter){
  m_ids <- sample(m_ids)
  r <- cor(melt(get_upper_tri(cor_mat_infancy[m_ids,m_ids]), 
                na.rm = TRUE)$value, m2_v, method = "spearman")
  rhos <- c(rhos,r)
}

perm_p_IO <- (sum(abs(true_rho) <= abs(rhos)) + 1)/(n_iter + 1) #two tailed test
perm_p_IO

### Infancy x Child

rhos <- c()
n_iter <- 5000
true_rho <- cor(infancy_upper$value,child_upper$value, method = "spearman")

# matrix permutation
m_ids <- colnames(cor_mat_infancy)
m2_v <- child_upper$value
for(iter in 1:n_iter){
  m_ids <- sample(m_ids)
  r <- cor(melt(get_upper_tri(cor_mat_infancy[m_ids,m_ids]), 
                na.rm = TRUE)$value, m2_v, method = "spearman")
  rhos <- c(rhos,r)
}

perm_p_IC <- (sum(abs(true_rho) <= abs(rhos)) + 1)/(n_iter + 1) #two tailed test
perm_p_IC

### Child x Old

rhos <- c()
n_iter <- 5000
true_rho <- cor(child_upper$value,old_upper$value, method = "spearman")

# matrix permutation
m_ids <- colnames(cor_mat_child)
m2_v <- old_upper$value
for(iter in 1:n_iter){
  m_ids <- sample(m_ids)
  r <- cor(melt(get_upper_tri(cor_mat_child[m_ids,m_ids]), 
                na.rm = TRUE)$value, m2_v, method = "spearman")
  rhos <- c(rhos,r)
}

perm_p_CO <- (sum(abs(true_rho) <= abs(rhos)) + 1)/(n_iter + 1) #two tailed test
perm_p_CO

#####################
##   Figure 6.1   ##
####################

## Parallel coordinates plot of all models and comparisons

ggparcoord(parcoords,columns = 2:16, groupColumn = 1, showPoints = T, 
  scale = "globalminmax") + geom_line(linewidth = 1) + geom_point(size = 1.5) + 
  scale_color_simpsons() + xlab("Growth Module Pairs") + 
  ylab("Correlation Value (r)") + theme_classic() + 
  theme(axis.text.x=element_text(face = "bold"), 
  axis.text.y = element_text(face = "bold"),
  axis.title = element_text(face = "bold"),
  legend.text = element_text(face = "bold"),
  legend.title = element_text(face = "bold"),
  legend.title.align = 0.5,
  legend.position = "bottom")


#####################
##   Figure 6.2   ##
####################

## Point plot comparing developmental stage data

lh_comp$LH_Stage <- factor(lh_comp$LH_Stage, ordered = T, 
                  levels = c("Infancy", "Childhood", "Juvenile & Adolescence"))

lh_comp %>% ggplot(aes(x = Value, y = Pairs, color = LH_Stage)) + 
  geom_point(size = 5) + facet_wrap(LH_Stage~., ncol=3) + 
  scale_color_simpsons() + labs(color = "Developmental Stage") + 
  xlab("Correlation (r)") + ylab("Module Pairs") + theme_minimal() + 
  theme(axis.text.x = element_text(face = "bold"),
  axis.text.y=element_text(face = "bold"), 
  axis.title = element_text(face = "bold"),
  legend.text = element_text(face = "bold"),
  legend.title = element_text(face = "bold"),
  legend.title.align = 0.5,
  legend.position = "",
  strip.text.x = element_text(color = "black", face = "bold"))


#####################
##   Figure 6.3   ##
####################

## Point plot comparing biological sex models

sex_comp$Model <- factor(sex_comp$Model, ordered = T, 
                         levels = c("Male", "Female", "Full Model"))


sex_comp %>% ggplot(aes(x = Value, y = Module.Pairs, color = Model)) + 
  geom_point(size = 5) + facet_wrap(Model~., ncol=3) +
  scale_color_manual(values = c("#D2AF81FF", "#FD7446FF", "#D5E4A2FF")) + 
  labs(color = "Model") + xlab("Correlation (r)") + ylab("Module Pairs") + 
  theme_minimal() + theme(axis.text.x=element_text(face = "bold"),
  axis.text.y=element_text(face = "bold"),
  axis.title = element_text(face = "bold"),
  legend.text = element_text(face = "bold"),
  legend.title = element_text(face = "bold"),
  legend.title.align = 0.5,
  legend.position = "",
  strip.text.x = element_text(
  color = "black", face = "bold"))

#####################
##   Table 6.2   ##
####################

sub <- melt(cor_mat_allvars) # wide to long format
comb <- paste0("[",sub$Var1,",",sub$Var2,"]") # new matrix of pairs
combos <- cbind(sub,comb)
combos <- combos %>% mutate(full1 = ifelse(Var1 %in% lbl, "LBL", 
      ifelse(Var1 %in% dent, "DD", ifelse(Var1 %in% pelvis, "PEL", 
      ifelse(Var1 %in% oss, "OSS", "EF")))), full2 = ifelse(Var2 %in% lbl, 
      "LBL", ifelse(Var2 %in% dent, "DD", ifelse(Var2 %in% pelvis, "PEL", 
      ifelse(Var2 %in% oss, "OSS", "EF")))),
      full_combo = paste0("[",full1, ",", full2, "]"))

combos <- combos %>% mutate(full1_sub = ifelse(grepl("DL", Var1), "LBL", 
      ifelse(grepl("MSB|PB|DB", Var1), "LBB", ifelse(grepl("max|man", Var1), 
      "DD", ifelse(grepl("_EF", Var1), "EF", "OSS")))),
      full2_sub = ifelse(grepl("DL", Var2), "LBL", ifelse(grepl("MSB|PB|DB", 
      Var2), "LBB", ifelse(grepl("max|man", Var2), "DD", ifelse(grepl("_EF", 
      Var2), "EF", "OSS")))), full_combo2 = paste0("[",full1_sub, ",", 
      full2_sub, "]"))

sub2 <- combos %>% filter(!value == 1.00000)
summ <- sub2 %>% group_by(full_combo2) %>% summarise(mean = mean(value))
summ_allvars <- summ %>% distinct(mean, .keep_all = TRUE)


#####################
##   Figure 6.4   ##
####################

## Extract upper triangle of the matrix
upper_tri <- get_upper_tri(cor_mat_allvars)

## Correlation Plot
ggplot(melt(upper_tri, na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() + 
  labs(x="",y="",fill="Corr") +
  theme(axis.text.x=element_text(size=7, angle=45, vjust=1, hjust=1, 
  margin=margin(-3,0,0,0), colour = col, face = "bold"),
  axis.text.y=element_text(size=7, margin=margin(0,-3,0,0), colour = col, 
  face = "bold"), panel.grid.major=element_blank(), 
  panel.grid = element_line(color="black"))

#####################
##   Figure 6.5   ##
####################

## Marginal relationships between strongest and weakest pairs
p1 <- q_out_all %>% ggplot() + geom_point(aes(UDL, RDL)) + 
  xlab("Ulna Diaphyseal Length") + ylab("Radius Diaphyseal Length") + 
  theme_classic()

p1 <- ggMarginal(p1, type = "histogram")


p2 <- q_out_all %>% ggplot() + geom_point(aes(UPE_EF, RDL)) + 
  xlab("Proximal Ulna Development Score") + ylab("Radius Diaphyseal Length") + 
  theme_classic()
p2 <- ggMarginal(p2, type = "histogram")

allvars_margin <- ggarrange(p1, p2, nrow= 2)
allvars_margin

#####################
##   Figure 6.6   ##
####################

ggplot(melt(upper_tri[1:18,1:18], na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color="black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") +
  theme(axis.text.x=element_text(size=10, angle=45, vjust=1, hjust=1, 
  margin=margin(-3,0,0,0), colour = "steelblue", face = "bold"),
  axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "steelblue", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), 
  color = "black", size = 2.5, fontface = "bold")


#####################
##   Figure 6.7   ##
####################

ggplot(melt(upper_tri[19:34,19:34], na.rm = TRUE), 
  aes(Var1, Var2, fill=value)) + geom_tile(height=0.8, width=0.8, 
  color="black") + scale_fill_gradient2(low="blue", mid="white", high="red", 
  midpoint = 0, limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "tomato", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "tomato", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", 
  size = 2.5, fontface = "bold")


#####################
##   Figure 6.8   ##
####################

ggplot(melt(upper_tri[35:54,35:54], na.rm = TRUE), 
  aes(Var1, Var2, fill=value)) + geom_tile(height=0.8, width=0.8, 
  color = "black") + scale_fill_gradient2(low="blue", mid="white", high="red", 
  midpoint = 0, limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "goldenrod", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "goldenrod", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", size = 2, 
  fontface = "bold")


#####################
##   Figure 6.9   ##
####################

## Extract upper triangle of the matrix
upper_tri <- get_upper_tri(cor_mat_infancy)

ggplot(melt(upper_tri, na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() + 
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=7, angle=45,
  vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = col, face = "bold"),
  axis.text.y=element_text(size=7, margin=margin(0,-3,0,0), colour = col, 
  face = "bold"), panel.grid.major=element_blank(), 
  panel.grid = element_line(color="black"))


#####################
##   Table 6.3   ##
####################

sub <- melt(cor_mat_infancy) # wide to long format
comb <- paste0("[",sub$Var1,",",sub$Var2,"]") # new matrix of pairs
combos <- cbind(sub,comb)
combos <- combos %>% mutate(full1 = ifelse(Var1 %in% lbl, "LBL", 
  ifelse(Var1 %in% dent, "DD", ifelse(Var1 %in% pelvis, "PEL", 
  ifelse(Var1 %in% oss, "OSS", "EF")))), full2 = ifelse(Var2 %in% lbl, "LBL", 
  ifelse(Var2 %in% dent, "DD", ifelse(Var2 %in% pelvis, "PEL", 
  ifelse(Var2 %in% oss, "OSS", "EF")))), full_combo = paste0("[",full1, ",", 
  full2, "]"))

combos <- combos %>% mutate(full1_sub = ifelse(grepl("DL", Var1), "LBL", 
  ifelse(grepl("MSB|PB|DB", Var1), "LBB", ifelse(grepl("max|man", Var1), "DD", 
  ifelse(grepl("_EF", Var1), "EF", "OSS")))), full2_sub = ifelse(grepl("DL", 
  Var2), "LBL", ifelse(grepl("MSB|PB|DB", Var2), "LBB", ifelse(grepl("max|man", 
  Var2), "DD", ifelse(grepl("_EF", Var2), "EF", "OSS")))), 
  full_combo2 = paste0("[",full1_sub, ",", full2_sub, "]"))

sub2 <- combos %>% filter(!value == 1.00000)
summ <- sub2 %>% group_by(full_combo2) %>% summarise(mean = mean(value))
summ_infancy <- summ %>% distinct(mean, .keep_all = TRUE)


#####################
##   Figure 6.10   ##
####################

## Marginal relationships between strongest and weakest pairs
p1 <- q_out_infancy %>% ggplot() + geom_point(aes(UDL, RDL)) + 
  xlab("Ulna Diaphyseal Length") + ylab("Radius Diaphyseal Length") + 
  
  theme_classic()

p1 <- ggMarginal(p1, type = "histogram")

p2 <- q_out_infancy %>% ggplot() + geom_point(aes(RDL, TDE_EF)) + 
  xlab("Radius Diaphyseal Length") + ylab("Distal Tibia Development") + 
  theme_classic()

p2 <- ggMarginal(p2, type = "histogram")

infancy_margin <- ggarrange(p1, p2, nrow= 2)
infancy_margin

######################
##   Figure 6.11   ##
####################

ggplot(melt(upper_tri[1:18,1:18], na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color="black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "steelblue", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "steelblue", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", 
  size = 2.5, fontface = "bold")


#####################
##   Figure 6.12   ##
####################

ggplot(melt(upper_tri[19:34,19:34], na.rm = TRUE), aes(Var1, Var2, 
  fill=value)) + geom_tile(height=0.8, width=0.8, color="black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "tomato", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "tomato", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", 
  size = 2.5, fontface = "bold")


#####################
##   Figure 6.13   ##
####################

ggplot(melt(upper_tri[35:54,35:54], na.rm = TRUE), aes(Var1, Var2, 
  fill=value)) + geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() +
  coord_equal() + labs(x="",y="",fill="Corr") +
  theme(axis.text.x=element_text(size=10, angle=45, vjust=1, hjust=1, 
  margin=margin(-3,0,0,0), colour = "goldenrod", face = "bold"),
  axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "goldenrod", face = "bold"),
  panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", size = 2,
  fontface = "bold")


#####################
##   Figure 6.14   ##
####################

## Extract upper triangle of the matrix
upper_tri <- get_upper_tri(cor_mat_child)

ggplot(melt(upper_tri, na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") +
  theme(axis.text.x=element_text(size=7, angle=45, vjust=1, hjust=1, 
  margin=margin(-3,0,0,0), colour = col, face = "bold"),
  axis.text.y=element_text(size=7, margin=margin(0,-3,0,0), colour = col, 
  face = "bold"), panel.grid.major=element_blank())


#####################
##   Table 6.4   ##
####################

sub <- melt(cor_mat_child) # wide to long format
comb <- paste0("[",sub$Var1,",",sub$Var2,"]") # newmatrix of pairs
combos <- cbind(sub,comb)
combos <- combos %>% mutate(full1 = ifelse(Var1 %in% lbl, "LBL", ifelse(Var1 %in% dent, "DD", ifelse(Var1 %in% pelvis, "PEL", ifelse(Var1 %in% oss, "OSS", "EF")))),
                            full2 = ifelse(Var2 %in% lbl, "LBL", ifelse(Var2 %in% dent, "DD", ifelse(Var2 %in% pelvis, "PEL", ifelse(Var2 %in% oss, "OSS", "EF")))),
                            full_combo = paste0("[",full1, ",", full2, "]"))

combos <- combos %>% mutate(full1_sub = ifelse(grepl("DL", Var1), "LBL", ifelse(grepl("MSB|PB|DB", Var1), "LBB", ifelse(grepl("max|man", Var1), "DD", ifelse(grepl("_EF", Var1), "EF", "OSS")))),
                            full2_sub = ifelse(grepl("DL", Var2), "LBL", ifelse(grepl("MSB|PB|DB", Var2), "LBB", ifelse(grepl("max|man", Var2), "DD", ifelse(grepl("_EF", Var2), "EF", "OSS")))),
                            full_combo2 = paste0("[",full1_sub, ",", full2_sub, "]"))

sub2 <- combos %>% filter(!value == 1.00000)
summ <- sub2 %>% group_by(full_combo2) %>% summarise(mean = mean(value)) %>% print(n = 30)
summ_child <- summ %>% distinct(mean, .keep_all = TRUE)

#####################
##   Figure 6.15   ##
####################

p1 <- q_out_child %>% ggplot() + geom_point(aes(UDL, RDL)) + 
  xlab("Ulna Diaphyseal Length") + ylab("Radius Diaphyseal Length") + 
  theme_classic()

p1 <- ggMarginal(p1, type = "histogram")

p2 <- q_out_child %>% ggplot() + geom_point(aes(TDL, man_PM2)) + 
  xlab("Tibia Diaphyseal Length") + 
  ylab("Mandibular 2nd Premolar Development") + theme_classic()

p2 <- ggMarginal(p2, type = "histogram")

ggarrange(p1, p2, nrow= 2)

#####################
##   Figure 6.16   ##
####################

ggplot(melt(upper_tri[1:18,1:18], na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "steelblue", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "steelblue", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", 
  size = 2.5, fontface = "bold")

#####################
##   Figure 6.17   ##
####################

ggplot(melt(upper_tri[19:34,19:34], na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "tomato", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "tomato", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", 
  size = 2.5, fontface = "bold")


#####################
##   Figure 6.18   ##
####################

ggplot(melt(upper_tri[35:54,35:54], na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "goldenrod", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "goldenrod", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", size = 2,
  fontface = "bold")

#####################
##   Figure 6.19   ##
####################

# Extract upper triangle of the matrix
upper_tri <- get_upper_tri(cor_mat_old)

ggplot(melt(upper_tri, na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=7, angle=45,
  vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = col, face = "bold"),
  axis.text.y=element_text(size=7, margin=margin(0,-3,0,0), colour = col, 
  face = "bold"), panel.grid.major=element_blank())

#####################
##   Table 6.5   ##
####################

sub <- melt(cor_mat_old) # wide to long format
comb <- paste0("[",sub$Var1,",",sub$Var2,"]") # newmatrix of pairs
combos <- cbind(sub,comb)
combos <- combos %>% mutate(full1 = ifelse(Var1 %in% lbl, "LBL", 
  ifelse(Var1 %in% dent, "DD", ifelse(Var1 %in% pelvis, "PEL", 
  ifelse(Var1 %in% oss, "OSS", "EF")))), full2 = ifelse(Var2 %in% lbl, "LBL", 
  ifelse(Var2 %in% dent, "DD", ifelse(Var2 %in% pelvis, "PEL", 
  ifelse(Var2 %in% oss, "OSS", "EF")))), full_combo = paste0("[",full1, ",", 
  full2, "]"))

combos <- combos %>% mutate(full1_sub = ifelse(grepl("DL", Var1), "LBL", 
  ifelse(grepl("MSB|PB|DB", Var1), "LBB", ifelse(grepl("max|man", Var1), "DD", 
  ifelse(grepl("_EF", Var1), "EF", "OSS")))), full2_sub = ifelse(grepl("DL", 
  Var2), "LBL", ifelse(grepl("MSB|PB|DB", Var2), "LBB", ifelse(grepl("max|man", 
  Var2), "DD", ifelse(grepl("_EF", Var2), "EF", "OSS")))), 
  full_combo2 = paste0("[",full1_sub, ",", full2_sub, "]"))

sub2 <- combos %>% filter(!value == 1.00000)
summ <- sub2 %>% group_by(full_combo2) %>% summarise(mean = mean(value))
summ_old <- summ %>% distinct(mean, .keep_all = TRUE)

#####################
##   Figure 6.20   ##
####################

## Marginal relationships between strongest and weakest pairs
p1 <- q_out_old %>% ggplot() + geom_point(aes(TPE_EF, FDE_EF)) + 
  xlab("Proximal Tibia Epiphysis Development") + 
  ylab("Distal Femur Epiphysis Development") + theme_classic()

p1 <- ggMarginal(p1, type = "histogram")


p2 <- q_out_old %>% ggplot() + geom_point(aes(FH_EF, max_I2)) + 
  xlab("Femoral Head Development") + ylab("Maxillary 2nd Incisor Development") + 
  theme_classic()

p2 <- ggMarginal(p2, type = "histogram")


ggarrange(p1, p2, nrow= 2)


#####################
##   Figure 6.21   ##
####################

ggplot(melt(upper_tri[1:18,1:18], na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "steelblue", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "steelblue", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", 
  size = 2.5, fontface = "bold")

#####################
##   Figure 6.22   ##
####################

ggplot(melt(upper_tri[19:34,19:34], na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "tomato", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "tomato", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", 
  size = 2.5, fontface = "bold")


#####################
##   Figure 6.23   ##
####################

ggplot(melt(upper_tri[35:54,35:54], na.rm = TRUE), aes(Var1, Var2, 
  fill=value)) + geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "goldenrod", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "goldenrod", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", size = 2,
  fontface = "bold")


#####################
##   Figure 6.24   ##
####################

# Extract upper triangle of the matrix
upper_tri <- get_upper_tri(cor_mat_male)

ggplot(melt(upper_tri, na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() + 
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=7, angle=45,
  vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = col, face = "bold"),
  axis.text.y=element_text(size=7, margin=margin(0,-3,0,0), colour = col, 
  face = "bold"), panel.grid.major=element_blank())


#####################
##   Table 6.6   ##
####################

sub <- melt(cor_mat_male) # wide to long format
comb <- paste0("[",sub$Var1,",",sub$Var2,"]") # newmatrix of pairs
combos <- cbind(sub,comb)
combos <- combos %>% mutate(full1 = ifelse(Var1 %in% lbl, "LBL", 
  ifelse(Var1 %in% dent, "DD", ifelse(Var1 %in% pelvis, "PEL", 
  ifelse(Var1 %in% oss, "OSS", "EF")))), full2 = ifelse(Var2 %in% lbl, "LBL", 
  ifelse(Var2 %in% dent, "DD", ifelse(Var2 %in% pelvis, "PEL", 
  ifelse(Var2 %in% oss, "OSS", "EF")))), full_combo = paste0("[",full1, ",", 
  full2, "]"))

combos <- combos %>% mutate(full1_sub = ifelse(grepl("DL", Var1), "LBL", 
  ifelse(grepl("MSB|PB|DB", Var1), "LBB", ifelse(grepl("max|man", Var1), "DD", 
  ifelse(grepl("_EF", Var1), "EF", "OSS")))), full2_sub = ifelse(grepl("DL", 
  Var2), "LBL", ifelse(grepl("MSB|PB|DB", Var2), "LBB", ifelse(grepl("max|man", 
  Var2), "DD", ifelse(grepl("_EF", Var2), "EF", "OSS")))), 
  full_combo2 = paste0("[",full1_sub, ",", full2_sub, "]"))

sub2 <- combos %>% filter(!value == 1.00000)
summ <- sub2 %>% group_by(full_combo2) %>% summarise(mean = mean(value))
summ_male <- summ %>% distinct(mean, .keep_all = TRUE)

#####################
##   Figure 6.26   ##
####################

p1 <- q_out_male %>% ggplot() + geom_point(aes(RDL, UDL)) + 
  xlab("Radius Diaphyseal Length") + ylab("Ulna Diaphyseal Length") + 
  theme_classic()

p1 <- ggMarginal(p1, type = "histogram")

p2 <- q_out_male %>% ggplot() + geom_point(aes(UDL, UPE_EF)) + 
  xlab("Ulna Diaphyseal Length") + ylab("Proximal Ulna Development") + 
  theme_classic()

p2 <- ggMarginal(p2, type = "histogram")


ggarrange(p1, p2, nrow= 2)

#####################
##   Figure 6.28   ##
####################

ggplot(melt(upper_tri[1:18,1:18], na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "steelblue", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "steelblue", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", 
  size = 2.5, fontface = "bold")


#####################
##   Figure 6.30   ##
####################

ggplot(melt(upper_tri[19:34,19:34], na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "tomato", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "tomato", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", 
  size = 2.5, fontface = "bold")


#####################
##   Figure 6.32   ##
####################

ggplot(melt(upper_tri[35:54,35:54], na.rm = TRUE), aes(Var1, Var2, 
  fill=value)) + geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "goldenrod", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "goldenrod", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", size = 2,
  fontface = "bold")


#####################
##   Figure 6.25   ##
####################

upper_tri <- get_upper_tri(cor_mat_female)

ggplot(melt(upper_tri, na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=7, angle=45,
  vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = col, face = "bold"),
  axis.text.y=element_text(size=7, margin=margin(0,-3,0,0), colour = col, 
  face = "bold"), panel.grid.major=element_blank())

#####################
##   Table 6.6   ##
####################

sub <- melt(cor_mat_female) # wide to long format
comb <- paste0("[",sub$Var1,",",sub$Var2,"]") # newmatrix of pairs
combos <- cbind(sub,comb)
combos <- combos %>% mutate(full1 = ifelse(Var1 %in% lbl, "LBL", 
  ifelse(Var1 %in% dent, "DD", ifelse(Var1 %in% pelvis, "PEL", 
  ifelse(Var1 %in% oss, "OSS", "EF")))), full2 = ifelse(Var2 %in% lbl, "LBL", 
  ifelse(Var2 %in% dent, "DD", ifelse(Var2 %in% pelvis, "PEL", 
  ifelse(Var2 %in% oss, "OSS", "EF")))),
  full_combo = paste0("[",full1, ",", full2, "]"))

combos <- combos %>% mutate(full1_sub = ifelse(grepl("DL", Var1), "LBL", 
  ifelse(grepl("MSB|PB|DB", Var1), "LBB", ifelse(grepl("max|man", Var1), "DD", 
  ifelse(grepl("_EF", Var1), "EF", "OSS")))), full2_sub = ifelse(grepl("DL", 
  Var2), "LBL", ifelse(grepl("MSB|PB|DB", Var2), "LBB", ifelse(grepl("max|man", 
  Var2), "DD", ifelse(grepl("_EF", Var2), "EF", "OSS")))), 
  full_combo2 = paste0("[",full1_sub, ",", full2_sub, "]"))

sub2 <- combos %>% filter(!value == 1.00000)
summ <- sub2 %>% group_by(full_combo2) %>% summarise(mean = mean(value))
summ_female <- summ %>% distinct(mean, .keep_all = TRUE)

#####################
##   Figure 6.27   ##
####################

# Marginal relationships between strongest and weakest pairs
p1 <- q_out_female %>% ggplot() + geom_point(aes(UDL, RDL)) + 
  xlab("Ulna Diaphyseal Length") + ylab("Radius Diaphyseal Length") + 
  theme_classic()

p1 <- ggMarginal(p1, type = "histogram")

p2 <- q_out_female %>% ggplot() + geom_point(aes(UDL, UPE_EF)) + 
  xlab("Ulna Diaphyseal Length") + ylab("Proximal Ulna Development") + 
  theme_classic()

p2 <- ggMarginal(p2, type = "histogram")


ggarrange(p1, p2, nrow= 2)

#####################
##   Figure 6.29   ##
####################

ggplot(melt(upper_tri[1:18,1:18], na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "steelblue", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "steelblue", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", 
  size = 2.5, fontface = "bold")

#####################
##   Figure 6.31   ##
####################

ggplot(melt(upper_tri[19:34,19:34], na.rm = TRUE), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "tomato", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "tomato", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", 
  size = 2.5, fontface = "bold")

#####################
##   Figure 6.33   ##
####################

ggplot(melt(upper_tri[35:54,35:54], na.rm = TRUE), aes(Var1, Var2, 
  fill=value)) + geom_tile(height=0.8, width=0.8, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
  limits = c(-1,1)) + theme_minimal() + coord_equal() +
  labs(x="",y="",fill="Corr") + theme(axis.text.x=element_text(size=10, 
  angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), colour = "goldenrod", 
  face = "bold"), axis.text.y=element_text(size=10, margin=margin(0,-3,0,0), 
  colour = "goldenrod", face = "bold"), panel.grid.major=element_blank()) + 
  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", size = 2,
  fontface = "bold")

##################################END###########################################
