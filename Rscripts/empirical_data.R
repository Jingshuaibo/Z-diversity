library(readxl)
library(writexl)
library(zetadiv)
library(dplyr)
library(ggplot2)
library(lmtest)
library(ggpubr)

#0.Preparation----
{#0.1 raw data----
#Functions for Z diversity calculation
load("Zdiv_function.RData")

#empirical species distribution data for birds (bird) and butterflies (butf)
bird_data <- read_excel("BirdData.xlsx")
butf_data <- read_excel("ButterflyData.xlsx")
bird_name <- read_excel("Bird_Info.xlsx")
butf_name <- read_excel("Butterfly_Info.xlsx")

#data process
#birds
bird_abun <- bird_data[c(1:3, 9:35)]
bird_abun$treatment[bird_abun$treatment==1] <- 'S1'
bird_abun$treatment[bird_abun$treatment==2] <- 'G5'
bird_abun$treatment[bird_abun$treatment==3] <- 'B16'
bird_abun$treatment[bird_abun$treatment==4] <- 'P32'

bird_abun$soil[bird_abun$soil==1] <- 'SL'
bird_abun$soil[bird_abun$soil==2] <- 'L'
bird_abun$soil[bird_abun$soil==3] <- 'CL'

bird_pres <- bird_abun
cbird <- c(4:30)
bird_pres[cbird] <- ifelse(bird_pres[cbird] > 0, 1, 0)

bird_sp <- colnames(bird_abun)[cbird]
nbird <- length(bird_sp)

#butterfly
butf_abun <- butf_data[c(1:3, 9:39)]
butf_abun$treatment[butf_abun$treatment==1] <- 'S1'
butf_abun$treatment[butf_abun$treatment==2] <- 'G5'
butf_abun$treatment[butf_abun$treatment==3] <- 'B16'
butf_abun$treatment[butf_abun$treatment==4] <- 'P32'

butf_abun$soil[butf_abun$soil==1] <- 'SL'
butf_abun$soil[butf_abun$soil==2] <- 'L'
butf_abun$soil[butf_abun$soil==3] <- 'CL'

butf_abun <- cbind(butf_abun[1:3], butf_abun[order(colnames(butf_abun)[-(1:3)])+3])

butf_pres <- butf_abun
cbutf <- c(4:34)
butf_pres[cbutf] <- ifelse(butf_pres[cbutf] > 0, 1, 0)

butf_sp <- colnames(butf_abun)[cbutf]
nbutf <- length(butf_sp)

#0.2 set variables----
soil <- c('SL', 'L', 'CL')
treat <- c('S1', 'G5', 'B16', 'P32')
site <- c('SL', 'L', 'CL', 'S1', 'G5', 'B16', 'P32', 'ALL')

col_soil <- c('#A32900','#FF6100','#FFAC4D')
col_treat <- c('#E5AB02','#A5761C','#66A61E','#4C6C43')
col_sp <- c('#308192','#B55489')

index_sp <- c('bird','butf')
index_plot <- c('soil', 'treat')
index_metric <- c('zdiv', 'zr')
index_group <- list(c(1:3), c(4:7))
index_color <- list(col_soil, col_treat)
index_function <- list(ZdivLine_Painter, ZrLine_Painter)

#0.3 functions----
#DataTemp: generate species-site matrix with specific plots
DataTemp <- function(dataset, site_i){
  if(site_i %in% c('SL','L','CL')) {
    data_temp <- filter(dataset, soil==site_i)[-c(1:3)]}
  else if(site_i %in% c('S1', 'G5', 'B16', 'P32')) {
    data_temp <- filter(dataset, treatment==site_i)[-c(1:3)]}
  else {data_temp <- dataset[-c(1:3)]}
  return(data_temp)
}

#ZdivLine_Painter: plot zeta diversity as line graph
ZdivLine_Painter <- function(zdiv_data, x_item, y_item, col_gp, ttl=NULL, lg_ttl=NULL){
  if(lg_ttl == 'soil'){
    x_bk <- seq(0, 16, 4)}
  else if(lg_ttl == 'treat'){
    x_bk <- seq(0, 12, 3)}
  else{ x_bk <- seq(0, 50, 10)}
  
  p_tp <- ggplot(data=zdiv_data, aes(x=get(x_item), y=get(y_item), color=Group, linetype=Model))+ 
    geom_point(size=1.5) + 
    geom_line(size=1.1) + 
    scale_color_manual(values=col_gp)+
    scale_linetype_manual(values = c('solid', 'dashed'))+
    scale_x_continuous(breaks = x_bk)+
    labs(x="Order", y="Zeta diversity", color=lg_ttl, title=ttl)+ 
    theme(axis.title = element_text(size=16, face='bold'),
          axis.text = element_text(size=15, face='bold'),
          plot.title=element_text(face="bold", size=15,  hjust=0.5),
          plot.subtitle=element_text(face="bold", size=13,  hjust=0.5),
          legend.title=element_text(face='bold', size=13),
          panel.border=element_rect(fill=NA, linewidth=1.2),
          panel.background = element_rect(color='white', fill=NA),
          panel.grid=element_line(colour='grey', linetype='dashed'))
  return(p_tp)
}

#ZdivLine_Painter: plot zeta ratios as line graph
ZrLine_Painter <- function(zr_data, x_item, y_item, col_gp, ttl=NULL, lg_ttl=NULL){
  if(lg_ttl == 'soil'){
    x_bk <- seq(0, 16, 4)}
  else if(lg_ttl == 'treat'){
    x_bk <- seq(0, 12, 3)}
  else{ x_bk <- seq(0, 50, 10)}
  
  p_tp <- ggplot(data=zr_data, aes(x=get(x_item), y=get(y_item), color=Group, linetype=Model))+ 
    geom_point(size=1.5) + 
    geom_line(size=1.1) + 
    scale_color_manual(values=col_gp)+
    scale_linetype_manual(values = c('solid', 'dashed'))+
    scale_x_continuous(breaks = x_bk)+
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2))+
    labs(x="Order", y="Zeta ratio", color=lg_ttl, title=ttl)+ 
    theme(axis.title = element_text(size=16, face='bold'),
          axis.text = element_text(size=15, face='bold'),
          plot.title=element_text(face="bold", size=15,  hjust=0.5),
          plot.subtitle=element_text(face="bold", size=13,  hjust=0.5),
          legend.title=element_text(face='bold', size=13),
          panel.border=element_rect(fill=NA, linewidth=1.2),
          panel.background = element_rect(color='white', fill=NA),
          panel.grid=element_line(colour='grey', linetype='dashed'))
  return(p_tp)
}

}

#1.Model fitting for the community----
{#1.1 Bird----
zeta_bird <- vector('list', length = 8)
names(zeta_bird) <- site
for (i in 1:8) {
  pres_tp <- DataTemp(bird_pres, site[i])
  pres_tp <- pres_tp[which(colSums(pres_tp) > 0)]
  zeta_bird[[i]] <- Zeta.decline.ex(pres_tp, orders = 1:nrow(pres_tp), plot = F)
}

qZeta_bird <- vector('list', length = 8)
names(qZeta_bird ) <- site
for (i in 1:8) {
  abun_tp <- DataTemp(bird_abun, site[i])
  abun_tp <- abun_tp[which(colSums(abun_tp) > 0)]
  qZeta_bird[[i]] <- Zdiv_Cal(abun_tp)
}

#summary
zmod_bird_smr <- data.frame(Species='Bird', Group=site, nSite=NA, qnSite=NA, 'QuaEX'=NA, 'QuaPL'=NA, 'QanEX'=NA, 
                            'QanPL'=NA, 'QanLNR'=NA, 'QanLOG'=NA, 'QanSS'=NA, 'QanSAS'=NA, 
                            'QuaMod'=NA, 'QanMod'=NA, 'QuaMod_exp'=NA, 'QanMod_exp'=NA)
for (i in 1:8) {
  zeta_tp <- zeta_bird[[i]]
  zmod_bird_smr$nSite[i] <- nrow(zeta_tp[["zeta.exp"]][["model"]])
  zmod_bird_smr$QuaEX[i] <- zeta_tp[["aic"]][["AIC"]][1]
  zmod_bird_smr$QuaPL[i] <- zeta_tp[["aic"]][["AIC"]][2]
  if(zeta_tp[["aic"]][["AIC"]][1] == -Inf & zeta_tp[["aic"]][["AIC"]][1] == -Inf){
    zmod_bird_smr$QuaMod[i] <- 'None'}
  else{zmod_bird_smr$QuaMod[i] <- ifelse(zeta_tp[["aic"]][["AIC"]][1] < zeta_tp[["aic"]][["AIC"]][2], 'EX', 'PL')}
  
  qzeta_tp <- qZeta_bird[[i]]
  zmod_bird_smr$qnSite[i] <- nrow(qzeta_tp[["DeclineModel"]][["zdivsum_models"]][["zdivsum_model_list"]][[1]][["model"]])
  zmod_bird_smr[i, 7:12] <- qzeta_tp[["DeclineModel"]][["zdivsum_models"]][["zdivsum_summary"]][1, 2:7]
  zmod_bird_smr$QanMod[i] <- qzeta_tp[["DeclineModel"]][["zdivsum_models"]][["zdivsum_summary"]][["Model"]][1]
}
zmod_bird_smr$QuaMod_exp <- ifelse(zmod_bird_smr$QuaMod=='EX', 'EX', 'NE')
zmod_bird_smr$QanMod_exp <- ifelse(zmod_bird_smr$QanMod=='EX', 'EX', 'NE')

write_xlsx(zmod_bird_smr, 'zmod_bird_smr.xlsx')

#1.2 Butterfly----
zeta_butf <- vector('list', length = 8)
names(zeta_butf) <- site
for (i in 1:8) {
  pres_tp <- DataTemp(butf_pres, site[i])
  pres_tp <- pres_tp[which(colSums(pres_tp) > 0)]
  zeta_butf[[i]] <- Zeta.decline.ex(pres_tp, orders = 1:nrow(pres_tp))
}

qZeta_butf <- vector('list', length = 8)
names(qZeta_butf ) <- site
for (i in 1:8) {
  abun_tp <- DataTemp(butf_abun, site[i])
  abun_tp <- abun_tp[which(colSums(abun_tp) > 0)]
  qZeta_butf[[i]] <- Zdiv_Cal(abun_tp)
}

#summary
zmod_butf_smr <- data.frame(Species='butf', Group=site, nSite=NA, qnSite=NA, 'QuaEX'=NA, 'QuaPL'=NA, 'QanEX'=NA, 
                            'QanPL'=NA, 'QanLNR'=NA, 'QanLOG'=NA, 'QanSS'=NA, 'QanSAS'=NA, 
                            'QuaMod'=NA, 'QanMod'=NA, 'QuaMod_exp'=NA, 'QanMod_exp'=NA)
for (i in 1:8) {
  zeta_tp <- zeta_butf[[i]]
  zmod_butf_smr$nSite[i] <- nrow(zeta_tp[["zeta.exp"]][["model"]])
  zmod_butf_smr$QuaEX[i] <- zeta_tp[["aic"]][["AIC"]][1]
  zmod_butf_smr$QuaPL[i] <- zeta_tp[["aic"]][["AIC"]][2]
  if(zeta_tp[["aic"]][["AIC"]][1] == -Inf & zeta_tp[["aic"]][["AIC"]][1] == -Inf){
    zmod_butf_smr$QuaMod[i] <- 'None'}
  else{zmod_butf_smr$QuaMod[i] <- ifelse(zeta_tp[["aic"]][["AIC"]][1] < zeta_tp[["aic"]][["AIC"]][2], 'EX', 'PL')}
  
  qzeta_tp <- qZeta_butf[[i]]
  zmod_butf_smr$qnSite[i] <- nrow(qzeta_tp[["DeclineModel"]][["zdivsum_models"]][["zdivsum_model_list"]][[1]][["model"]])
  zmod_butf_smr[i, 7:12] <- qzeta_tp[["DeclineModel"]][["zdivsum_models"]][["zdivsum_summary"]][1, 2:7]
  zmod_butf_smr$QanMod[i] <- qzeta_tp[["DeclineModel"]][["zdivsum_models"]][["zdivsum_summary"]][["Model"]][1]
}
zmod_butf_smr$QuaMod_exp <- ifelse(zmod_butf_smr$QuaMod=='EX', 'EX', 'NE')
zmod_butf_smr$QanMod_exp <- ifelse(zmod_butf_smr$QanMod=='EX', 'EX', 'NE')

write_xlsx(zmod_butf_smr, 'zmod_butf_smr.xlsx')

#1.3 model summary----
site_list <- vector(mode = 'list', length = 7)
names(site_list) <- site[-8]

#for qualitative data
zeta_bird_list <- list(zdiv=site_list, zr=site_list)
for (i in 1:2) {
  for (j in 1:7) {
    index_tp <- ifelse(i==1, "zeta.val", "ratio")
    data_tp <- zeta_bird[[j]][[index_tp]]
    zeta_bird_list[[i]][[j]] <- data.frame(Group=site[j], Order=c(1:length(data_tp)), 
                                           Model=zmod_bird_smr$QuaMod_exp[j], Value=data_tp)
    zeta_bird_list[[i]][[j]]$Value[is.nan(zeta_bird_list[[i]][[j]]$Value)] <- 0
}}
zeta_butf_list <- list(zdiv=site_list, zr=site_list)
for (i in 1:2) {
  for (j in 1:7) {
    index_tp <- ifelse(i==1, "zeta.val", "ratio")
    data_tp <- zeta_butf[[j]][[index_tp]]
    zeta_butf_list[[i]][[j]] <- data.frame(Group=site[j], Order=c(1:length(data_tp)), 
                                           Model=zmod_butf_smr$QuaMod_exp[j], Value=data_tp)
    zeta_butf_list[[i]][[j]]$Value[is.nan(zeta_butf_list[[i]][[j]]$Value)] <- 0
}}

#for quantitative data
qzeta_bird_list <- list(zdiv=site_list, zr=site_list)
for (i in 1:2) {
  for (j in 1:7) {
    if(i == 1){
      data_tp <- qZeta_bird[[j]][["Zdiv"]][["Zdiv_sum"]]}
    else{data_tp <- qZeta_bird[[j]][["Zratio"]][["Zr_sum"]][1:(length(qZeta_bird[[j]][["Zdiv"]][["Zdiv_sum"]])-1)]}
    qzeta_bird_list[[i]][[j]] <- data.frame(Group=site[j], Order=c(1:length(data_tp)), 
                                           Model=zmod_bird_smr$QanMod_exp[j], Value=data_tp)
}}
qzeta_butf_list <- list(zdiv=site_list, zr=site_list)
for (i in 1:2) {
  for (j in 1:7) {
    if(i == 1){
      data_tp <- qZeta_butf[[j]][["Zdiv"]][["Zdiv_sum"]]}
    else{data_tp <- qZeta_butf[[j]][["Zratio"]][["Zr_sum"]][1:(length(qZeta_butf[[j]][["Zdiv"]][["Zdiv_sum"]])-1)]}
    qzeta_butf_list[[i]][[j]] <- data.frame(Group=site[j], Order=c(1:length(data_tp)), 
                                            Model=zmod_butf_smr$QanMod_exp[j], Value=data_tp)
}}

#1.4 plots----
group_list <- vector(mode = 'list', length = 2)
names(group_list) <- index_plot 

for (metric in c('zeta', 'qzeta')) {
  for (sp in index_sp) {
    datalist_tp <- get(paste0(metric,'_',sp,'_list'))
    plist_tp <- list(zdiv=group_list, zr=group_list)
    for (i in 1:2){
      for (j in 1:2) {
        pname_tp <- paste0(sp,'_',metric,'_',index_metric[i],'_',index_plot[j])
        data_tp <-  datalist_tp[[ index_metric[i] ]][[ index_group[[j]][1] ]]
        for (k in index_group[[j]][-1]) {
          data_tp <- rbind(data_tp, datalist_tp[[ index_metric[i] ]][[k]]) }
        data_tp$Group <- factor(data_tp$Group, levels = site)
        data_tp$Model <- factor(data_tp$Model, levels = c('NE','EX'))
        plist_tp[[i]][[j]] <- index_function[[i]](data_tp, x_item='Order', y_item='Value', col_gp=index_color[[j]], 
                                                  lg_ttl=index_plot[j], ttl=pname_tp)
        ggsave(paste0(pname_tp,'.png'), plist_tp[[i]][[j]], height=1105, width=1100, units="px")
      }}
    assign(paste0('p_',sp,'_',metric,'_list'), plist_tp)
    }}
}

#2.Model fitting for the single species----
{#2.1 model summary----
#bird
bird_spmd <- list(SL=NA, L=NA, CL=NA, S1=NA, G5=NA, B16=NA, P32=NA, ALL=NA)
for (i in 1:8) {
  data_tp <- data.frame(bird_name, nSite=NA, Model=NA, AIC=NA, Rsq=NA)
  md_tp <- qZeta_bird[[i]][["DeclineModel"]][["species_models"]]
  
  for (j in 1:nrow(data_tp)){
    spname_tp <- data_tp$ShortName[j]
    if(spname_tp %in% md_tp[["species_aic_summary"]]$Species){
      data_tp$Model[j] <- filter(md_tp[["species_aic_summary"]], Species==spname_tp)$Model }
    else {data_tp$Model[j] <- 'NULL'}
    if(data_tp$Model[j] %in% c('None', 'NULL')){ 
      data_tp$AIC[j] <- 'NULL'
      data_tp$Rsq[j] <- 'NULL'
      data_tp$nSite[j] <- 0 }
    else{
      data_tp$AIC[j] <- filter(md_tp[["species_aic_summary"]], Species==spname_tp)[data_tp$Model[j]]
      data_tp$Rsq[j] <- filter(md_tp[["species_rsq_summary"]], Species==spname_tp)[data_tp$Model[j]]
      data_tp$nSite[j] <- nrow(md_tp[["species_model_list"]][[spname_tp]][[1]][["model"]]) }
  }
  bird_spmd[[i]] <- data_tp
}

bird_spmd_smr <- data.frame(bird_name, SL=NA, L=NA, CL=NA, S1=NA, G5=NA, B16=NA, P32=NA, ALL=NA)
for(i in 1:8){
  md_tp <- bird_spmd[[i]]$Model
  md_tp[bird_spmd[[i]]$nSite <= 3] <- 'NULL'
  md_tp[md_tp=='None'] <- 'None'
  bird_spmd_smr[i+3] <- md_tp
}

qzeta_bird_sp_list <- vector(mode='list', length=nbird)
names(qzeta_bird_sp_list) <- bird_name$ShortName
for (k in 1:nbird) {
  zlist_tp <- list(zdiv=site_list, zr=site_list)
  sp_tp <- bird_name$ShortName[k]
  for (i in 1:2) {
    for (j in 1:7) {
      norder_tp <- ifelse(i==1, length(qZeta_bird[[j]][["Zdiv"]][[1]]), length(qZeta_bird[[j]][["Zdiv"]][[1]])-1)
      zlist_tp[[i]][[j]] <- data.frame(Group=site[j], Order=c(1:norder_tp), Value=0, Model='NULL')
      if(sp_tp %in% names(qZeta_bird[[j]][["Zdiv"]])){
        if(i == 1){
          zlist_tp[[i]][[j]]$Value <- qZeta_bird[[j]][["Zdiv"]][[sp_tp]]}
        else{zlist_tp[[i]][[j]]$Value <- qZeta_bird[[j]][["Zratio"]][[sp_tp]][1:norder_tp] }
        
        md_tp <- filter(bird_spmd_smr, ShortName==sp_tp)[[site[j]]]
        md_tp <- ifelse(md_tp %in% c('PL','LN','LG','SS','GP'), 'NE', md_tp)
        zlist_tp[[i]][[j]]$Model <- rep(md_tp, norder_tp) }
      else{ }
    }}
  qzeta_bird_sp_list[[k]] <- zlist_tp
}

#butf
butf_spmd <- list(SL=NA, L=NA, CL=NA, S1=NA, G5=NA, B16=NA, P32=NA, ALL=NA)
for (i in 1:8) {
  data_tp <- data.frame(butf_name, nSite=NA, Model=NA, AIC=NA, Rsq=NA)
  md_tp <- qZeta_butf[[i]][["DeclineModel"]][["species_models"]]
  
  for (j in 1:nrow(data_tp)){
    spname_tp <- data_tp$ShortName[j]
    if(spname_tp %in% md_tp[["species_aic_summary"]]$Species){
      data_tp$Model[j] <- filter(md_tp[["species_aic_summary"]], Species==spname_tp)$Model }
    else {data_tp$Model[j] <- 'NULL'}
    if(data_tp$Model[j] %in% c('None', 'NULL')){ 
      data_tp$AIC[j] <- 'NULL'
      data_tp$Rsq[j] <- 'NULL'
      data_tp$nSite[j] <- 0 }
    else{
      data_tp$AIC[j] <- filter(md_tp[["species_aic_summary"]], Species==spname_tp)[data_tp$Model[j]]
      data_tp$Rsq[j] <- filter(md_tp[["species_rsq_summary"]], Species==spname_tp)[data_tp$Model[j]]
      data_tp$nSite[j] <- nrow(md_tp[["species_model_list"]][[spname_tp]][[1]][["model"]]) }
  }
  butf_spmd[[i]] <- data_tp
}

butf_spmd_smr <- data.frame(butf_name, SL=NA, L=NA, CL=NA, S1=NA, G5=NA, B16=NA, P32=NA, ALL=NA)
for(i in 1:8){
  md_tp <- butf_spmd[[i]]$Model
  md_tp[butf_spmd[[i]]$nSite <= 3] <- 'NULL'
  md_tp[md_tp=='None'] <- 'None'
  butf_spmd_smr[i+3] <- md_tp
}

qzeta_butf_sp_list <- vector(mode='list', length=nbutf)
names(qzeta_butf_sp_list) <- butf_name$ShortName
for (k in 1:nbutf) {
  zlist_tp <- list(zdiv=site_list, zr=site_list)
  sp_tp <- butf_name$ShortName[k]
  for (i in 1:2) {
    for (j in 1:7) {
      norder_tp <- ifelse(i==1, length(qZeta_butf[[j]][["Zdiv"]][[1]]), length(qZeta_butf[[j]][["Zdiv"]][[1]])-1)
      zlist_tp[[i]][[j]] <- data.frame(Group=site[j], Order=c(1:norder_tp), Value=0, Model='NULL')
      if(sp_tp %in% names(qZeta_butf[[j]][["Zdiv"]])){
        if(i == 1){
          zlist_tp[[i]][[j]]$Value <- qZeta_butf[[j]][["Zdiv"]][[sp_tp]]}
        else{zlist_tp[[i]][[j]]$Value <- qZeta_butf[[j]][["Zratio"]][[sp_tp]][1:norder_tp] }
        
        md_tp <- filter(butf_spmd_smr, ShortName==sp_tp)[[site[j]]]
        md_tp <- ifelse(md_tp %in% c('PL','LN','LG','SS','GP'), 'NE', md_tp)
        zlist_tp[[i]][[j]]$Model <- rep(md_tp, norder_tp) }
      else{ }
    }}
  qzeta_butf_sp_list[[k]] <- zlist_tp
}

#summary
spmd_summary <- rbind(bird_spmd_smr[-11], butf_spmd_smr[-11])
spmd_summary$Tax <- c(rep('bird', nbird), rep('butterfly', nbutf))
spmd_summary <- spmd_summary[c(11,1:10)]
for (i in 1:nrow(spmd_summary)) {
  if(all(spmd_summary[i,5:11] == rep('NULL',7))){
    spmd_summary$tag[i] <- 0 }
  else{spmd_summary$tag[i] <- 1}
}
spmd_summary <- filter(spmd_summary, tag==1)[-12]

write_xlsx(bird_spmd_smr, 'bird_spmd_smr.xlsx')
write_xlsx(butf_spmd_smr, 'butf_spmd_smr.xlsx')
write_xlsx(spmd_summary, 'spmd_summary.xlsx')
  
#2.2 plots for 3 typical species----
sp_repre <- c('chgr', 'cipl', 'dapl')
p_sp_qzeta_list <- list(chgr=list(zdiv=group_list, zr=group_list), 
                        cipl=list(zdiv=group_list, zr=group_list), 
                        dapl=list(zdiv=group_list, zr=group_list))
for (i in 1:3) {
  sp_tp <- ifelse(i==3, 'butf', 'bird')
  datalist_tp <- get(paste0('qzeta_', sp_tp, '_sp_list'))
  
  for (j in 1:2){
    for (k in 1:2) {
      data_tp <-  datalist_tp[[ sp_repre[i] ]][[j]][[ index_group[[k]][1] ]]
      pname_tp <- paste0(sp_tp,'_',sp_repre[i],'_',index_metric[j],'_',index_plot[k])
      for (l in index_group[[k]][-1]) {
        data_tp <- rbind(data_tp, datalist_tp[[ sp_repre[i] ]][[j]][[l]]) }
      data_tp$Model[which(data_tp$Model=='NULL')] <- 'EX'
      data_tp$Group <- factor(data_tp$Group, levels = site)
      data_tp$Model <- factor(data_tp$Model, levels = c('NE','EX'))
      
      p_sp_qzeta_list[[i]][[j]][[k]] <- index_function[[j]](data_tp, x_item='Order', y_item='Value', col_gp=index_color[[k]], 
                                                            lg_ttl=index_plot[k], ttl=pname_tp)
      ggsave(paste0(pname_tp,'.png'), p_sp_qzeta_list[[i]][[j]][[k]], height=1105, width=1100, units="px")
    }}}

}