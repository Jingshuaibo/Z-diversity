#Z diversity calculation functions***

#1.Order_Matrix Function----
#reorder the abundance of each species across sites in descending order

Order_Matrix <- function(mat_unorder){
  mat_tp <- mat_unorder
  for (i in 1:ncol(mat_tp)) {
    mat_tp[i] <- mat_tp[order(mat_tp[[i]], decreasing=T), i]
  }
  return(mat_tp)
}

#2. qF_Cal Function----
#calculate qF components at each order for the overall community (qF_sum) and each species

qF_Cal <- function(mat, species_name=NULL){
  nSp <- ncol(mat)
  nUnit <- nrow(mat)
  mat_tp <- mat
  if(!is.null(species_name)){ 
    sp_name <- species_name }
  else{
    sp_name <- colnames(mat) }
  
  for (i in 1:nSp) {
    for (j in 1:(nUnit-1)) {
      mat_tp[j,i] <- mat_tp[j,i]-mat_tp[j+1,i]
    }}
  mat_tp <- cbind(mat_tp, 'qF_sum'=rowSums(mat_tp))
  
  qF_list <- vector(mode = 'list', length = nSp+1)
  names(qF_list) <- c(sp_name, 'qF_sum')
  for(i in 1:(nSp+1)){
    qF_list[[i]] <- mat_tp[[i]]
    names(qF_list[[i]]) <- paste0('qF_', 1:nUnit)
  }
  
  return(qF_list)
}

#3. CP_Matrix Function----
#generate binomial coefficients for Z diversity calculation

CP_Matrix <- function(dem){
  mat_tp <- data.frame(matrix(ncol = dem, nrow = dem))
  for (i in 1:dem) {
    for (j in 1:dem) {
      if(i >= j){
        mat_tp[i,j] <- choose(i, j)/choose(dem, j) }
      else{
        mat_tp[i,j] <- 0}
    }}
  return(mat_tp)
}

#4.Zdiv_qF_Cal Function----
#calculate Z diversity components at each order for the overall community (Zdiv_sum) and each species

Zdiv_qF_Cal <- function(F_list){
  nSp <- length(F_list)-1
  nUnit <- length(F_list[[1]])
  sp_name <- names(F_list)[1:nSp]
  mat_cp <- CP_Matrix(nUnit)
  
  Zdiv_list <- vector(mode = 'list', length = nSp+1)
  names(Zdiv_list) <- c(sp_name, 'Zdiv_sum')
  for (i in 1:(nSp+1)) {
    Zdiv_list[[i]] <- rep(NA, nUnit)
    names(Zdiv_list[[i]]) <- paste0('Zdiv_', 1:nUnit)
    for (j in 1:nUnit) {
      Zdiv_list[[i]][j] <- sum(F_list[[i]] * mat_cp[[j]])
    }}
  
  return(Zdiv_list)
}

#5.Zr_Cal Function----
#calculate Z ratios at each order for the overall community (Zr_sum) and each species
#meanZr: the average of Z ratios across all orders
#Zd (Z distance): 1 - Z ratio, increases with the decrease of Z ratio, 
#so represents the compositional difference (ecological distance) among sites 

Zr_Cal <- function(Zdiv_list){
  nSp <- length(Zdiv_list)-1
  nUnit <- length(Zdiv_list[[1]])
  sp_name <- names(Zdiv_list)[1:nSp]
  
  Zr_list <- vector(mode = 'list', length = nSp+1)
  names(Zr_list) <- c(sp_name, 'Zr_sum')
  for (i in 1:(nSp+1)) {
    Zr_list[[i]] <- rep(NA, nUnit+1)
    for (j in 1:(nUnit-1)) {
      if(Zdiv_list[[i]][j] == 0){
        Zr_list[[i]][j] <- 0 }
      else{
        Zr_list[[i]][j] <- Zdiv_list[[i]][j+1]/Zdiv_list[[i]][j] }
    }
    Zr_list[[i]][nUnit] <- mean(Zr_list[[i]][1:(nUnit-1)])
    Zr_list[[i]][nUnit+1] <- 1-Zr_list[[i]][nUnit]
    names(Zr_list[[i]]) <- c(paste0('Zr_',2:nUnit,'/',1:(nUnit-1)), 'meanZr', 'Zd')
  }
  
  return(Zr_list)
}

#1.6 Dcl_Model Function----
#model fitting of the Z diversity decline forms
#include 6 models for selection: EX(exponential), PL(power-law), LN(linear), LG(logarithmic), SS(logistic model), GP(Gompertz)
#species_aic_summary & species_rsq_summary: the AIC and R2 of all models and the optimal model for each species
#zdivsum_summary: the AIC and R2 of all models and the optimal model for the whole community

Dcl_Model <- function(Z_list){
  nSp <- length(Z_list)-1
  nUnits <- length(Z_list[[1]])
  sp_name <- names(Z_list)[1:nSp]
  
  #model fitting for each species
  species_model_list <- vector(mode = 'list', length = nSp)
  names(species_model_list) <- c(sp_name)
  for (i in 1:nSp) {
    data_tp <- data.frame('zdiv'=Z_list[[i]], 'x_order'=1:length(Z_list[[i]]))
    data_tp$ss_z <- log10(2*max(data_tp$zdiv)/data_tp$zdiv-1)
    data_tp$sas_z <- log10(log10(2*max(data_tp$zdiv))-log10(data_tp$zdiv))
    
    exp_model <- lm(log10(zdiv) ~ x_order, na.action=stats::na.omit, data=filter(data_tp, zdiv > 0))
    pl_model <- lm(log10(zdiv) ~ log10(x_order), na.action=stats::na.omit, data=filter(data_tp, zdiv > 0))
    lnr_model <- lm(zdiv ~ x_order, na.action=stats::na.omit, data=filter(data_tp, zdiv > 0))
    log_model <- lm(zdiv ~ log10(x_order), na.action=stats::na.omit, data=filter(data_tp, zdiv > 0))
    ss_model <- lm(ss_z ~ x_order, na.action=stats::na.omit, data=filter(data_tp, ss_z != Inf, zdiv > 0))
    sas_model <- lm(sas_z ~ log10(x_order), na.action=stats::na.omit, data=filter(data_tp, sas_z != Inf, zdiv > 0))
    species_model_list[[i]] <- list(exp_model, pl_model, lnr_model, log_model, ss_model, sas_model)
  }
  
  #calculate AIC and R2
  species_aic_summary <- data.frame('Species'=sp_name, 'EX'=NA, 'PL'=NA, 'LN'=NA, 'LG'=NA, 'SS'=NA, 'GP'=NA, 'Model'=NA)
  species_rsq_summary <- data.frame('Species'=sp_name, 'EX'=NA, 'PL'=NA, 'LN'=NA, 'LG'=NA, 'SS'=NA, 'GP'=NA, 'Model'=NA)
  model_name <- c('EX', 'PL', 'LN', 'LG', 'SS', 'GP')
  for (i in 1:nSp) {
    for (j in 1:6) {
      species_aic_summary[i, j+1] <- AIC(species_model_list[[i]][[j]])
      species_rsq_summary[i, j+1] <- summary(species_model_list[[i]][[j]])[["adj.r.squared"]]
      species_rsq_summary[i, j+1] <- ifelse(is.nan(species_rsq_summary[i, j+1]), 0, species_rsq_summary[i, j+1])
    }
    
    if(sum(species_aic_summary[i,2:7] == -Inf) <= 1){
      species_aic_summary$Model[i] <- model_name[which.min(species_aic_summary[i,2:7])]}
    else{species_aic_summary$Model[i] <- 'None'}
    
    max_R2 <- max(species_rsq_summary[i,2:7])
    if(is.na(max_R2)){
      species_rsq_summary$Model[i] <- 'None'}
    else{
      if(sum(species_rsq_summary[i,2:7]==max_R2) > 1){
        species_rsq_summary$Model[i] <- 'None'}
      else{species_rsq_summary$Model[i] <- model_name[which.max(species_rsq_summary[i,2:7])]}
    }
  }
  species_models <- list(species_model_list, species_aic_summary, species_rsq_summary)
  names(species_models) <- c('species_model_list', 'species_aic_summary', 'species_rsq_summary')
  
  #model fitting for the overall community
  data_tp <- data.frame('zdiv'=Z_list[[nSp+1]], 'x_order'=1:length(Z_list[[nSp+1]]))
  data_tp$ss_z <- log10(2*max(data_tp$zdiv)/data_tp$zdiv-1)
  data_tp$sas_z <- log10(log10(2*max(data_tp$zdiv))-log10(data_tp$zdiv))
  
  exp_model <- lm(log10(zdiv) ~ x_order, na.action=stats::na.omit, data=filter(data_tp, zdiv > 0))
  pl_model <- lm(log10(zdiv) ~ log10(x_order), na.action=stats::na.omit, data=filter(data_tp, zdiv > 0))
  lnr_model <- lm(zdiv ~ x_order, na.action=stats::na.omit, data=filter(data_tp, zdiv > 0))
  log_model <- lm(zdiv ~ log10(x_order), na.action=stats::na.omit, data=filter(data_tp, zdiv > 0))
  ss_model <- lm(ss_z ~ x_order, na.action=stats::na.omit, data=filter(data_tp, ss_z != Inf, zdiv > 0))
  sas_model <- lm(sas_z ~ log10(x_order), na.action=stats::na.omit, data=filter(data_tp, sas_z != Inf, zdiv > 0))
  species_model_list[[i]] <- list(exp_model, pl_model, lnr_model, log_model, ss_model, sas_model)
  
  zdivsum_model_list <- list(exp_model, pl_model, lnr_model, log_model, ss_model, sas_model)
  zdivsum_summary <- data.frame(Index=c('AIC','R2'), 'EX'=NA, 'PL'=NA, 'LN'=NA, 'LG'=NA, 'SS'=NA, 'GP'=NA, 'Model'=NA)
  for (i in 1:6) {
    zdivsum_summary[1, i+1] <- AIC(zdivsum_model_list[[i]])
    zdivsum_summary[2, i+1] <- summary(zdivsum_model_list[[i]])[["adj.r.squared"]]
    zdivsum_summary[2, j+1] <- ifelse(is.nan(zdivsum_summary[2, j+1]), 0, zdivsum_summary[2, j+1])
  }
  
  if(sum(zdivsum_summary[1, 2:7] == -Inf) <= 1){
    zdivsum_summary$Model[1] <- model_name[which.min(zdivsum_summary[1,2:7])]}
  else{zdivsum_summary$Model[1] <- 'None'}
  
  max_R2 <- max(zdivsum_summary[2,2:7])
  if(is.na(max_R2)){
    zdivsum_summary$Model[2] <- 'None'}
  else{
    if(sum(zdivsum_summary[2,2:7] == max_R2) > 1){
      zdivsum_summary$Model[2] <- 'None'}
    else{zdivsum_summary$Model[2] <- model_name[which.max(zdivsum_summary[2, 2:7])]}
  }
  
  zdivsum_models <- list(zdivsum_model_list, zdivsum_summary)
  names(zdivsum_models) <- c('zdivsum_model_list', 'zdivsum_summary')
  
  model <- list(species_models, zdivsum_models)
  names(model) <- c('species_models', 'zdivsum_models')
  
  return(model)
}

#7.Zdiv_Cal Function----
#execute Z diversity analysis using the above 6 functions
#input: an abundance matrix with species as columns and sites as rows
#output: results as a list data. It includes:
#data_order: the abundance matrix after reordering the abundances of each species across sites in descending order
#data_qF: qF components for the community and each species
#data_zdiv: Z diversity components for the community and each species
#data_zr: Z ratios for the community and each species
#data_md: model fitting results for the community and each species

Zdiv_Cal <- function(abun_data){
  if(!require("dplyr")){
    print("Please install necessary R packages: dplyr")}
  else{}
  
  data_order <- Order_Matrix(abun_data)
  data_qF <- qF_Cal(data_order)
  data_zdiv <- Zdiv_qF_Cal(data_qF)
  data_zr <- Zr_Cal(data_zdiv)
  data_md <- Dcl_Model(data_zdiv)
  data_summary <- list(abundance=abun_data,
                       abun_order=data_order,
                       qF=data_qF,
                       Zdiv=data_zdiv,
                       Zratio=data_zr,
                       DeclineModel=data_md)
  
  return(data_summary)
}