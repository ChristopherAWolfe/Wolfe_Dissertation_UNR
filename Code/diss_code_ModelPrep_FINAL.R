################################################################################
### The following code prepares the data  for modeling in cmdstanr.           ##
### It takes as input the analysis type, sex, and developmental stage.        ##
### Analysis can be full, sex, or stage. Sex can be male or female,           ##
## stage can be infancy, childhood, and old (juv/adol). 'dat' is a            ## 
## character string of the file path to the requiste data file.               ##
################################################################################

library(tidyverse)

prep <- function(analysis, dat, sex = NA, stage = NA){
  
  
  
  dat <- read.csv(dat)
  
  if(analysis == "Full"){
    
    N <- nrow(dat)
    J <- ncol(dat[,48:65])
    K_dentition <- ncol(dat[,12:27])
    K_ef <- ncol(dat[,28:43])
    K_pelvis <- ncol(dat[,46:47])
    K_all <- K_dentition + K_ef + K_pelvis + 1 + 1
    continuous <- dat %>% dplyr::select(48:65) 
    y_continuous <- continuous[!is.na(continuous)]
    present <- length(y_continuous)
    missing <- sum(is.na(continuous))
    index_pres = which(!is.na(continuous), arr.ind = TRUE)
    index_miss = which(is.na(continuous), arr.ind = TRUE)
    dentition <- dat %>% dplyr::select(12:27)
    dentition[is.na(dentition)] <- 99
    y_dentition <- as.matrix(dentition)
    pelvis <- dat %>% dplyr::select(46:47)
    pelvis[is.na(pelvis)] <- 99
    y_pelvis <- as.matrix(pelvis)
    carpal <- dat %>% dplyr::select(44)
    carpal[is.na(carpal)] <- 99
    y_carpal <- c(carpal$CC_Oss)
    tarsal <- dat %>% dplyr::select(45)
    tarsal[is.na(tarsal)] <- 99
    y_tarsal <- c(tarsal$TC_Oss)
    ef <- dat %>% dplyr::select(28:43)
    ef[is.na(ef)] <- 99
    y_ef <- as.matrix(ef)
    X <- dat$agey
    
    stan_dat <- list(N = N,
                     J = J,
                     K_dentition = K_dentition,
                     K_ef = K_ef,
                     K_pelvis = K_pelvis,
                     K_all = K_all,
                     present = present,
                     y_continuous = y_continuous,
                     missing = missing,
                     index_pres = index_pres,
                     index_miss = index_miss,
                     y_dentition = y_dentition,
                     y_ef = y_ef,
                     y_pelvis = y_pelvis,
                     y_carpal = y_carpal,
                     y_tarsal = y_tarsal,
                     X = X)
    
  } else if(analysis == "Sex"){
    
    male <- dat %>% filter(SEX == "M") %>% droplevels()
    female <- dat %>% filter(SEX == "F") %>% droplevels()
    
    if(sex == "male"){
      
      dat <- male
      N <- nrow(dat)
      J <- ncol(dat[,48:65])
      K_dentition <- ncol(dat[,12:27])
      K_ef <- ncol(dat[,28:43])
      K_pelvis <- ncol(dat[,46:47])
      K_all <- K_dentition + K_ef + K_pelvis + 1 + 1
      continuous <- dat %>% dplyr::select(48:65)
      y_continuous <- continuous[!is.na(continuous)]
      present <- length(y_continuous)
      missing <- sum(is.na(continuous))
      index_pres = which(!is.na(continuous), arr.ind = TRUE)
      index_miss = which(is.na(continuous), arr.ind = TRUE)
      dentition <- dat %>% dplyr::select(12:27)
      dentition[is.na(dentition)] <- 99
      y_dentition <- as.matrix(dentition)
      pelvis <- dat %>% dplyr::select(46:47)
      pelvis[is.na(pelvis)] <- 99
      y_pelvis <- as.matrix(pelvis)
      carpal <- dat %>% dplyr::select(44)
      carpal[is.na(carpal)] <- 99
      y_carpal <- c(carpal$CC_Oss)
      tarsal <- dat %>% dplyr::select(45)
      tarsal[is.na(tarsal)] <- 99
      y_tarsal <- c(tarsal$TC_Oss)
      ef <- dat %>% dplyr::select(28:43)
      ef[is.na(ef)] <- 99
      y_ef <- as.matrix(ef)
      X <- dat$agey
      
      stan_dat <- list(N = N,
                       J = J,
                       K_dentition = K_dentition,
                       K_ef = K_ef,
                       K_pelvis = K_pelvis,
                       K_all = K_all,
                       present = present,
                       y_continuous = y_continuous,
                       missing = missing,
                       index_pres = index_pres,
                       index_miss = index_miss,
                       y_dentition = y_dentition,
                       y_ef = y_ef,
                       y_pelvis = y_pelvis,
                       y_carpal = y_carpal,
                       y_tarsal = y_tarsal,
                       X = X)
      
    } else {
      
      dat <- female
      N <- nrow(dat)
      J <- ncol(dat[,48:65])
      K_dentition <- ncol(dat[,12:27])
      K_ef <- ncol(dat[,28:43])
      K_pelvis <- ncol(dat[,46:47])
      K_all <- K_dentition + K_ef + K_pelvis + 1 + 1
      continuous <- dat %>% dplyr::select(48:65)
      y_continuous <- continuous[!is.na(continuous)]
      present <- length(y_continuous)
      missing <- sum(is.na(continuous))
      index_pres = which(!is.na(continuous), arr.ind = TRUE)
      index_miss = which(is.na(continuous), arr.ind = TRUE)
      dentition <- dat %>% dplyr::select(12:27)
      dentition[is.na(dentition)] <- 99
      y_dentition <- as.matrix(dentition)
      pelvis <- dat %>% dplyr::select(46:47)
      pelvis[is.na(pelvis)] <- 99
      y_pelvis <- as.matrix(pelvis)
      carpal <- dat %>% dplyr::select(44)
      carpal[is.na(carpal)] <- 99
      y_carpal <- c(carpal$CC_Oss)
      tarsal <- dat %>% dplyr::select(45)
      tarsal[is.na(tarsal)] <- 99
      y_tarsal <- c(tarsal$TC_Oss)
      ef <- dat %>% dplyr::select(28:43)
      ef[is.na(ef)] <- 99
      y_ef <- as.matrix(ef)
      X <- dat$agey
      
      stan_dat <- list(N = N,
                       J = J,
                       K_dentition = K_dentition,
                       K_ef = K_ef,
                       K_pelvis = K_pelvis,
                       K_all = K_all,
                       present = present,
                       y_continuous = y_continuous,
                       missing = missing,
                       index_pres = index_pres,
                       index_miss = index_miss,
                       y_dentition = y_dentition,
                       y_ef = y_ef,
                       y_pelvis = y_pelvis,
                       y_carpal = y_carpal,
                       y_tarsal = y_tarsal,
                       X = X)
      
    }
    
  } else {
    
    infancy <- dat %>% filter(agey >= 0 & agey <=2.999)
    child <- dat %>% filter(agey >= 3.0 & agey <=6.999)
    juv_adol <- dat %>% filter(agey >= 7.0)
    
    if(stage == "infancy"){
      
      dat <- infancy
      N <- nrow(dat)
      J <- ncol(dat[,48:65])
      K_dentition <- ncol(dat[,12:27])
      K_ef <- ncol(dat[,28:43])
      K_pelvis <- ncol(dat[,46:47])
      K_all <- K_dentition + K_ef + K_pelvis + 1 + 1
      continuous <- dat %>% dplyr::select(48:65)
      y_continuous <- continuous[!is.na(continuous)]
      present <- length(y_continuous)
      missing <- sum(is.na(continuous))
      index_pres = which(!is.na(continuous), arr.ind = TRUE)
      index_miss = which(is.na(continuous), arr.ind = TRUE)
      dentition <- dat %>% dplyr::select(12:27)
      dentition[is.na(dentition)] <- 99
      y_dentition <- as.matrix(dentition)
      pelvis <- dat %>% dplyr::select(46:47)
      pelvis[is.na(pelvis)] <- 99
      y_pelvis <- as.matrix(pelvis)
      carpal <- dat %>% dplyr::select(44)
      carpal[is.na(carpal)] <- 99
      y_carpal <- c(carpal$CC_Oss)
      tarsal <- dat %>% dplyr::select(45)
      tarsal[is.na(tarsal)] <- 99
      y_tarsal <- c(tarsal$TC_Oss)
      ef <- dat %>% dplyr::select(28:43)
      ef[is.na(ef)] <- 99
      y_ef <- as.matrix(ef)
      X <- dat$agey
      
      stan_dat <- list(N = N,
                       J = J,
                       K_dentition = K_dentition,
                       K_ef = K_ef,
                       K_pelvis = K_pelvis,
                       K_all = K_all,
                       present = present,
                       y_continuous = y_continuous,
                       missing = missing,
                       index_pres = index_pres,
                       index_miss = index_miss,
                       y_dentition = y_dentition,
                       y_ef = y_ef,
                       y_pelvis = y_pelvis,
                       y_carpal = y_carpal,
                       y_tarsal = y_tarsal,
                       X = X)
      
    } else if(stage == "child"){
      
      dat <- child
      N <- nrow(dat)
      J <- ncol(dat[,48:65])
      K_dentition <- ncol(dat[,12:27])
      K_ef <- ncol(dat[,28:43])
      K_pelvis <- ncol(dat[,46:47])
      K_all <- K_dentition + K_ef + K_pelvis + 1 + 1
      continuous <- dat %>% dplyr::select(48:65)
      y_continuous <- continuous[!is.na(continuous)]
      present <- length(y_continuous)
      missing <- sum(is.na(continuous))
      index_pres = which(!is.na(continuous), arr.ind = TRUE)
      index_miss = which(is.na(continuous), arr.ind = TRUE)
      dentition <- dat %>% dplyr::select(12:27)
      dentition[is.na(dentition)] <- 99
      y_dentition <- as.matrix(dentition)
      pelvis <- dat %>% dplyr::select(46:47)
      pelvis[is.na(pelvis)] <- 99
      y_pelvis <- as.matrix(pelvis)
      carpal <- dat %>% dplyr::select(44)
      carpal[is.na(carpal)] <- 99
      y_carpal <- c(carpal$CC_Oss)
      tarsal <- dat %>% dplyr::select(45)
      tarsal[is.na(tarsal)] <- 99
      y_tarsal <- c(tarsal$TC_Oss)
      ef <- dat %>% dplyr::select(28:43)
      ef[is.na(ef)] <- 99
      y_ef <- as.matrix(ef)
      X <- dat$agey
      
      stan_dat <- list(N = N,
                       J = J,
                       K_dentition = K_dentition,
                       K_ef = K_ef,
                       K_pelvis = K_pelvis,
                       K_all = K_all,
                       present = present,
                       y_continuous = y_continuous,
                       missing = missing,
                       index_pres = index_pres,
                       index_miss = index_miss,
                       y_dentition = y_dentition,
                       y_ef = y_ef,
                       y_pelvis = y_pelvis,
                       y_carpal = y_carpal,
                       y_tarsal = y_tarsal,
                       X = X)
      
    } else {
      
      dat <- juv_adol
      N <- nrow(dat)
      J <- ncol(dat[,48:65])
      K_dentition <- ncol(dat[,12:27])
      K_ef <- ncol(dat[,28:43])
      K_pelvis <- ncol(dat[,46:47])
      K_all <- K_dentition + K_ef + K_pelvis + 1 + 1
      continuous <- dat %>% dplyr::select(48:65)
      y_continuous <- continuous[!is.na(continuous)]
      present <- length(y_continuous)
      missing <- sum(is.na(continuous))
      index_pres = which(!is.na(continuous), arr.ind = TRUE)
      index_miss = which(is.na(continuous), arr.ind = TRUE)
      dentition <- dat %>% dplyr::select(12:27)
      dentition[is.na(dentition)] <- 99
      y_dentition <- as.matrix(dentition)
      pelvis <- dat %>% dplyr::select(46:47)
      pelvis[is.na(pelvis)] <- 99
      y_pelvis <- as.matrix(pelvis)
      carpal <- dat %>% dplyr::select(44)
      carpal[is.na(carpal)] <- 99
      y_carpal <- c(carpal$CC_Oss)
      tarsal <- dat %>% dplyr::select(45)
      tarsal[is.na(tarsal)] <- 99
      y_tarsal <- c(tarsal$TC_Oss)
      ef <- dat %>% dplyr::select(28:43)
      ef[is.na(ef)] <- 99
      y_ef <- as.matrix(ef)
      X <- dat$agey
      
      stan_dat <- list(N = N,
                       J = J,
                       K_dentition = K_dentition,
                       K_ef = K_ef,
                       K_pelvis = K_pelvis,
                       K_all = K_all,
                       present = present,
                       y_continuous = y_continuous,
                       missing = missing,
                       index_pres = index_pres,
                       index_miss = index_miss,
                       y_dentition = y_dentition,
                       y_ef = y_ef,
                       y_pelvis = y_pelvis,
                       y_carpal = y_carpal,
                       y_tarsal = y_tarsal,
                       X = X)
      
    }
  }
}

##########################END###################################################
