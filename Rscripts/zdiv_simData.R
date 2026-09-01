library(readxl)
library(writexl)
library(ggplot2)
library(dplyr)
library(zetadiv)
library(pheatmap)
library(ggforce)
library(ggpubr)
library(ggplotify)
library(rsq)
library(lmtest)

#0.Preparation----
{#read data----
##simulated abundance data 
load("abun_list.RData") 
##simulated occurrence data
load("occur_list.RData")  
##quantitative beta diversity for each abundance matrix
qbeta_summary <- read_excel("qbeta_summary.xlsx")
##beta diversity for each occurrence matrix
beta_summary <- read_excel("beta_summary.xlsx") 
##Z diversity calculation functions
load("Zdiv_function.RData")

#Sign_TS Function----
#This function can convert p-value to corresponding significance symbols ('*') 
Sign_TS <- function(p){
  if(!is.numeric(p) | is.nan(p)){
    s <- 'NULL'}
  else{if(p >= 0.05) {s <- 'ns'}
    else if(p < 0.05 && p >= 0.01) {s <- '*'}
    else if(p < 0.01 && p >= 0.001) {s <- '**'}
    else if(p < 0.001 && p >= 0.0001) {s <- '***'}
    else {s <- '****'}}
  return(s)
}

#Function for plot heat maps----
Heatmap_Painter <- function(table_data, title=NULL, ant_r=NULL, ant_c=NULL, ant_color=NULL, 
                            col_bk=NULL, lg_bk=NULL, lg_lab=NULL, gap_r=NULL, gap_c=NULL, 
                            cell_color=colorRampPalette(c('white','red'))(101), dpnum=T, 
                            show_r=T, show_c=T, show_legend=T){
  p_heatmap <- pheatmap(table_data, main = title, fontsize = 11, 
                        show_rownames = show_r, show_colnames = show_c, angle_col=0, row_names_side='left',
                        display_numbers = dpnum, number_format = "%.2f", number_color = "black", fontsize_number=12, 
                        color = cell_color, breaks = col_bk, 
                        cellheight = 40, cellwidth = 40, border_color = "black", 
                        scale = "none", cluster_rows = F, cluster_cols = F, 
                        gaps_row = gap_r, gaps_col = gap_c, #设置间隔
                        annotation_col = ant_c, annotation_row = ant_r, #添加注释
                        annotation_colors = ant_color, #设置注释颜色
                        legend = show_legend, legend_breaks = lg_bk, legend_labels = lg_lab)
  return(p_heatmap)
}

#Function for plot scatter diagrams----
Point_painter <- function(data, x_item, y_item, group=NULL){
  rho <- filter(sprm_test, Beta==x_item, Zeta==y_item)[c(3,5)]
  p_value <- filter(sprm_test, Beta==x_item, Zeta==y_item)$Spearman_p
  p_tp <- ggplot(data=div_smr)+
    geom_abline(slope=1, intercept=0, linetype='dashed')+
    geom_point(mapping=aes(x=get(x_item), y=get(y_item)), alpha=0.5, size=1.7)+
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(limits=c(0,1)) +
    #annotate("text", x=0.12, y=1, color='black', size=3.5, family = "sans",
    #        label = paste0('rho = ', round(rho[[1]],2), ' (', rho[[2]], ')')) +
    labs(x=x_item, y=y_item, colour=group) + 
    theme(panel.border = element_rect(fill=NA, linewidth=1),
          panel.background = element_rect(color='white',fill=NA)) +
    theme(axis.title = element_text(size=12, face='bold'),
          axis.text = element_text(size=10, color='black'),
          axis.ticks.length=unit(-0.15, "cm"))
  return(p_tp)
}

#Winner_TS----
##This function convert the winner (i.e., the metric with a higher zeta ratio) to corresponding value
##If winner is: 1) quantitative zeta ratio -> 1; 
##              2) qualitative zeta ratio -> -1;
##              3) None (equal value) -> 0
Winner_TS <- function(winner){
  s <- ifelse(winner=='None', 0, ifelse(winner=='qZr', 1, -1))
  return(s)
}

}

#1.Model_fitting----
model_name <- c('EX', 'PL', 'LN', 'LG', 'SS', 'GP')

#qualitative zeta diversity
zmod_pres_list <- data_list
for (x in 1:6) {
  for (y in 1:6) {
    data_tp <- data_list[[x]][[y]]
    nsite_tp <- nrow(data_tp)
    zmod_pres_list[[x]][[y]] <- Zeta.decline.ex(data_tp, orders=1:nsite_tp, plot=F)
  }}

zmod_pres_list2 <- data_list
for (x in 1:6) {
  for (y in 1:6) {
    zmod_pres_list2[[x]][[y]] <- Zdiv_Cal(data_list[[x]][[y]])
  }}

#quantitative zeta diversity
zmod_abun_list <- abun_list
for (x in 1:3) { for (y in 1:3) { 
    for (i in 1:6) { for (j in 1:6) {
      zmod_abun_list[[x]][[y]][[i]][[j]] <- Zdiv_Cal(abun_list[[x]][[y]][[i]][[j]]) 
    }}}}

#2.Model_analysis----
{#2.1 Model number statistics----
#qualitative models
pres_aic_mod <- data.frame(Nest=rep(1:6, each=6), Turn=rep(1:6, 6), nSite=NA, EX=NA, PL=NA, Model=NA, Rsq=NA)
for (i in 1:nrow(pres_aic_mod)) {
  x <- pres_aic_mod$Nest[i]
  y <- pres_aic_mod$Turn[i]
  zeta_tp <- zmod_pres_list[[x]][[y]]
  
  pres_aic_mod$nSite[i] <- nrow(zeta_tp[["zeta.exp"]][["model"]])
  pres_aic_mod$EX[i] <- zeta_tp[["aic"]][["AIC"]][1]
  pres_aic_mod$PL[i] <- zeta_tp[["aic"]][["AIC"]][2]

  if(pres_aic_mod$EX[i] < pres_aic_mod$PL[i]){
    pres_aic_mod$Model[i] <- 'EX'
    pres_aic_mod$Rsq[i] <- summary(zeta_tp[["zeta.exp"]])[["adj.r.squared"]] }
  else if(pres_aic_mod$EX[i] > pres_aic_mod$PL[i]){
    pres_aic_mod$Model[i] <- 'PL'
    pres_aic_mod$Rsq[i] <- summary(zeta_tp[["zeta.pl"]])[["adj.r.squared"]] }
  else{pres_aic_mod$Model[i] <- 'None'
    pres_aic_mod$Rsq[i] <- NaN }
}

pres_aic_mod2 <- data.frame(Nest=rep(1:6, each=6), Turn=rep(1:6, 6), nSite=NA, 
                            EX=NA, PL=NA, LN=NA, LG=NA, SS=NA, GP=NA, Model=NA, Rsq=NA)
for (i in 1:nrow(pres_aic_mod)) {
  x <- pres_aic_mod$Nest[i]
  y <- pres_aic_mod$Turn[i]
  zeta_tp <- zmod_pres_list2[[x]][[y]][["DeclineModel"]][["zdivsum_models"]]
  
  pres_aic_mod2$nSite[i] <- nrow(zeta_tp[["zdivsum_model_list"]][[1]][["model"]])
  pres_aic_mod2[i, 4:10] <- zeta_tp[["zdivsum_summary"]][1, 2:8]
  if(pres_aic_mod2$Model[i] == 'None'){
    pres_aic_mod2$Rsq[i] <- NaN }
  else{
    pres_aic_mod2$Rsq[i] <- zeta_tp[["zdivsum_summary"]][which(colnames(zeta_tp[["zdivsum_summary"]])== pres_aic_mod2$Model[i])][2,1]}
}

#quantitative models
abun_aic_mod <- data.frame(qbeta_summary[1:4], nSite=NA, EX=NA, PL=NA, LN=NA, LG=NA, SS=NA, GP=NA, Model=NA, Rsq=NA)
for (k in 1:nrow(abun_aic_mod)) {
  x <- abun_aic_mod$Grad_row[k]
  y <- abun_aic_mod$Grad_col[k]
  i <- abun_aic_mod$Nest[k]
  j <- abun_aic_mod$Turn[k]
  zeta_tp <- zmod_abun_list[[x]][[y]][[i]][[j]][["DeclineModel"]][["zdivsum_models"]]
  
  abun_aic_mod$nSite[k] <- nrow(zeta_tp[["zdivsum_model_list"]][[1]][["model"]])
  abun_aic_mod[k, 6:12] <- zeta_tp[["zdivsum_summary"]][1, 2:8]
  if(abun_aic_mod$Model[k] == 'None'){
    abun_aic_mod$Rsq[k] <- NaN }
  else{
    abun_aic_mod$Rsq[k] <- zeta_tp[["zdivsum_summary"]][which(colnames(zeta_tp[["zdivsum_summary"]])==abun_aic_mod$Model[k])][2,1]}
}

#count model numbers
#for qualitative models
mdcount_pres <- data.frame('Model'=c('EX', 'PL', 'None'), 'Count'=NA, 'Pcent'=NA)
for (i in 1:3) {
  mdcount_pres$Count[i] <- nrow(filter(pres_aic_mod, Model==mdcount_pres$Model[i]))
}
mdcount_pres$Pcent <- mdcount_pres$Count/36
mdcount_pres$Model <- factor(mdcount_pres$Model, levels = c('EX', 'PL', 'None'))

mdcount_pres2 <- data.frame('Model'=c(model_name, 'None'), 'Count'=NA, 'Pcent'=NA)
for (i in 1:7) {
  mdcount_pres2$Count[i] <- nrow(filter(pres_aic_mod2, Model==mdcount_pres2$Model[i]))
}
mdcount_pres2$Pcent <- mdcount_pres2$Count/36
mdcount_pres2$Model <- factor(mdcount_pres2$Model, levels = c(model_name, 'None'))

#for quantitative models
mdcount_abun <- data.frame('Model'=c(model_name, 'None'), 'Count'=NA, 'Pcent'=NA)
for (i in 1:7) {
  mdcount_abun$Count[i] <- nrow(filter(abun_aic_mod, Model==mdcount_abun$Model[i]))
}
mdcount_abun$Pcent <- mdcount_abun$Count/324
mdcount_abun$Model <- factor(mdcount_abun$Model, levels = c(model_name, 'None'))

#summary
mdcount_smr <- data.frame('Model'=c('EX', 'NE', 'None'), Pres=NA, Abun=NA)
for (i in c(1,3)) {
  md_tp <- mdcount_smr$Model[i]
  mdcount_smr$Pres[i] <- filter(mdcount_pres, Model==md_tp)$Pcent %>% sum()
  mdcount_smr$Abun[i] <- filter(mdcount_abun, Model==md_tp)$Pcent %>% sum()
}
mdcount_smr[2,2:3] <- 1-colSums(mdcount_smr[-2,2:3])

#R2 summary 
abun_model_rsq_count <- data.frame(Rsq=c(0.99, 0.95, 0.9, 0.85, 0.8, 0.75), Count=NA, Pcent=NA)
for (i in 1:6) {
  rsq_tp <- abun_aic_mod$Rsq[!is.nan(abun_aic_mod$Rsq)]
  abun_model_rsq_count$Count[i] <- sum(rsq_tp > abun_model_rsq_count$Rsq[i])
}
abun_model_rsq_count$Pcent <- abun_model_rsq_count$Count/324

write_xlsx(pres_aic_mod, 'pres_model_summary.xlsx')
write_xlsx(abun_aic_mod, 'abun_model_summary.xlsx')
write_xlsx(mdcount_pres, 'pres_model_count.xlsx')
write_xlsx(mdcount_abun, 'abun_model_count.xlsx')
write_xlsx(abun_model_rsq_count, 'abun_model_rsq_count.xlsx')

#2.2 plot fan diagrams----
pres_mdcolor <- c('#23BAC5','#00A14E','grey')
abun_mdcolor <- c('#23BAC5','#00A14E','#B092B6','#e38d26','#FFC000')

p_pie_pres <- ggplot(mdcount_pres, aes(x=2, y=Pcent, fill=Model)) +
                geom_bar(stat="identity", color="black") + 
                coord_polar("y", start = 0)+
                scale_fill_manual(values = pres_mdcolor)+
                xlim(1.545, 2.5)+
                theme_void()
             
p_pie_abun <- ggplot(filter(mdcount_abun, Pcent>0), aes(x=2, y=Pcent, fill=Model)) +
                geom_bar(stat="identity", color="black") + 
                coord_polar("y", start = 0)+
                scale_fill_manual(values = abun_mdcolor)+
                xlim(1.545, 2.5)+
                theme_void()

p_pie_pres2 <- ggplot(filter(mdcount_smr[1:2], Pres>0), aes(x=2, y=Pres, fill=Model)) +
  geom_bar(stat="identity", color="black") + 
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c('#23BAC5','#00A14E','grey'))+
  xlim(1.545, 2.5)+
  theme_void()

p_pie_abun2 <- ggplot(filter(mdcount_smr[c(1,3)], Abun>0), aes(x=2, y=Abun, fill=Model)) +
  geom_bar(stat="identity", color="black") + 
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c('#23BAC5','#00A14E','grey'))+
  xlim(1.545, 2.5)+
  theme_void()

ggsave('p_pie_pres.png', p_pie_pres, height=2000, width=2000, units='px')
ggsave('p_pie_abun.png', p_pie_abun, height=2000, width=2000, units='px')
ggsave('p_pie_pres2.png', p_pie_pres2, height=2000, width=2000, units='px')
ggsave('p_pie_abun2.png', p_pie_abun2, height=2000, width=2000, units='px')

#2.3 model comparison----
model_comp <- data.frame(abun_aic_mod[1:4], OccurMod=rep(pres_aic_mod$Model, 9), AbunMod=abun_aic_mod$Model, OccurMod_exp=NA, AbunMod_exp=NA, Result=NA)
model_comp$OccurMod_exp <- ifelse(model_comp$OccurMod == 'EX', 'E', 'NE')
model_comp[which(model_comp$OccurMod=='None'), 'OccurMod_exp'] <- 'NA'
model_comp$AbunMod_exp <- ifelse(model_comp$AbunMod == 'EX', 'E', 'NE')
model_comp[which(model_comp$AbunMod=='None'), 'AbunMod'] <- 'NE'

for (k in 1:nrow(model_comp)) {
  md_tp <- c(model_comp$OccurMod_exp[k], model_comp$AbunMod_exp[k])
  if(all(md_tp == c('E', 'E'))){
    model_comp$Result[k] <- 1}
  else if(all(md_tp == c('NE', 'NE'))){
    model_comp$Result[k] <- 2}
  else if(all(md_tp==c('NE', 'E'))){
    model_comp$Result[k] <- 3}
  else if(all(md_tp==c('E', 'NE'))){
    model_comp$Result[k] <- 4}
  else{model_comp$Result[k] <- 0}
}

#summary
mdcount_comp <- data.frame(OccurMod=c('None','E','NE','NE','E'), AbunMod=c('*','E','NE','E','NE'), Count=NA, Pcent=NA)
for (i in 1:5) {
  mdcount_comp$Count[i] <- sum(model_comp$Result == (i-1))
}
mdcount_comp$Pcent <- mdcount_comp$Count/324

write_xlsx(model_comp, 'model_comp.xlsx')
write_xlsx(mdcount_comp, 'mdcount_comp.xlsx')

#2.4 plot heat maps----
#generate 36*36 tables for heat map
table_rowname <- table_colname <- rep(NA, 18)
for (x in 1:3) {
  for (y in 1:6) {
    table_rowname[6*(x-1)+y] <- paste0('G',x,'N',y)
    table_colname[6*(x-1)+y] <- paste0('G',x,'T',y)
  }}

ht_table <- data.frame(matrix(ncol=18, nrow=18))
rownames(ht_table) <- table_rowname
colnames(ht_table) <- table_colname

#heat map for model comparison
model_ant_table <- model_table <- ht_table
for (k in 1:nrow(model_comp)) {
  i <- (model_comp$Grad_row[k]-1)*6 + model_comp$Nest[k]
  j <- (model_comp$Grad_col[k]-1)*6 + model_comp$Turn[k]
  
  res_tp <- model_comp$Result[k]
  if(res_tp %in% c(1,2)) {model_table[i,j] <-1}
  else if(res_tp %in% c(3,4)) {model_table[i,j] <-2}
  else {model_table[i,j] <- 0}
  
  presmd_tp <- ifelse(model_comp$OccurMod[k]=='None', 'NA', model_comp$OccurMod[k])
  model_ant_table[i,j] <- paste0(presmd_tp, '/', abunmd_tp)
}

model_heatmap <- pheatmap(model_table, main = "Occurence/Abundance Model", fontsize = 11, 
                          show_rownames = T, show_colnames = T, angle_col=0, row_names_side='left',
                          display_numbers = model_ant_table, number_color = "black", fontsize_number=10.5, 
                          color=c('#BFBFBF', '#70AD47', '#F4B183'), breaks=seq(-0.5,2.5,1), 
                          cellheight = 42, cellwidth = 42, border_color = "black", 
                          scale = "none", cluster_rows = F, cluster_cols = F, 
                          gaps_row = seq(6,18,6), gaps_col = seq(6,18,6), 
                          legend = F, legend_breaks = c(0,1,2), 
                          legend_labels = c('NA/*','Same', 'Different'))
ggsave('model_comp_heatmap.png', model_heatmap, width = 4000, height = 3800, units = 'px')
}

#3.Zeta ratio comparison----
{#3.1 Zeta ratio comparison----
#for occurrence data
oZr_summary <- data.frame(Nest=rep(1:6, each=6), Turn=rep(1:6, 6), oZr=NA, oZr1=NA, oZr4=NA, oZr7=NA)
for(i in 1:36){
  x <- oZr_summary$Nest[i]
  y <- oZr_summary$Turn[i]
  zr_tp <- zmod_pres_list[[x]][[y]][["ratio"]]
  zr_tp[is.nan(zr_tp)] <- 0
  
  oZr_summary$oZr[i] <- mean(zr_tp)
  oZr_summary$oZr1[i] <- mean(zr_tp[1:3])
  oZr_summary$oZr4[i] <- mean(zr_tp[4:6])
  oZr_summary$oZr7[i] <- mean(zr_tp[7:9])
}

#for abundance data
qZr_summary <- data.frame(qbeta_summary[1:4], qZr=NA, qZr1=NA, qZr4=NA, qZr7=NA)
for (k in 1:nrow(qZr_summary)) {
  x <- qZr_summary$Grad_row[k]
  y <- qZr_summary$Grad_col[k]
  i <- qZr_summary$Nest[k]
  j <- qZr_summary$Turn[k]
  zr_tp <- zmod_abun_list[[x]][[y]][[i]][[j]][["Zratio"]][["Zr_sum"]]
  
  qZr_summary$qZr[k] <- zr_tp['meanZr']
  qZr_summary$qZr1[k] <- mean(zr_tp[1:3])
  qZr_summary$qZr4[k] <- mean(zr_tp[4:6])
  qZr_summary$qZr7[k] <- mean(zr_tp[7:9])
}

#结果汇总
zr_summary <- data.frame(qZr_summary[1:4], oZr=NA, oZr1=NA, oZr4=NA, oZr7=NA, qZr_summary[5:8])
for (i in 3:6){
  zr_summary[i+2] <- rep(oZr_summary[[i]], 9)
}
zr_summary$Grad_row <- factor(zr_summary$Grad_row, levels = c(1,2,3))
zr_summary$Grad_col <- factor(zr_summary$Grad_col, levels = c(1,2,3))

write_xlsx(oZr_summary, 'oZr_summary.xlsx')
write_xlsx(qZr_summary, 'qZr_summary.xlsx')
write_xlsx(zr_summary, 'zr_summary.xlsx')

#3.2 normality test----
#1vs36 (fix nestedness and turnover)
zr_noltest1 <- data.frame(Zr=rep(c('qZr','qZr1','qZr4','qZr7'), each=36), 
                         Nest=rep(1:6, each=6), Turn=rep(1:6), Shapiro_p=NA, Sign=NA)
for (i in 1:nrow(zr_noltest1)) {
  i_zd <- zr_noltest1$Zr[i]
  i_nest <- zr_noltest1$Nest[i]
  i_turn <- zr_noltest1$Turn[i]
  data_tp <- filter(zr_summary, Nest==i_nest, Turn==i_turn)[[i_zd]]
  
  if(var(data_tp) != 0){
    zr_noltest1$Shapiro_p[i] <- shapiro.test(data_tp)$p.value
    zr_noltest1$Sign[i] <- Sign_TS(zr_noltest1$Shapiro_p[i])}
  else{
    zr_noltest1$Shapiro_p[i] <- NaN
    zr_noltest1$Sign[i] <- 'NULL'}
}

#9vs9 (fix inter-site and inter-species abundance gradient)
zr_noltest2 <- data.frame(Zr=rep(c('qZr','qZr1','qZr4','qZr7'), each=9), 
                          Grad_row=rep(1:3, each=3), Grad_col=rep(1:3), Shapiro_p=NA, Sign=NA)
for (i in 1:nrow(zr_noltest2)) {
  i_zr <- zr_noltest2$Zr[i]
  i_grar <- zr_noltest2$Grad_row[i]
  i_grac <- zr_noltest2$Grad_col[i]
  data_tp <- filter(zr_summary, Grad_row==i_grar, Grad_col==i_grac)[[i_zr]]
  if(var(data_tp) != 0){
    zr_noltest2$Shapiro_p[i] <- shapiro.test(data_tp)$p.value
    zr_noltest2$Sign[i] <- Sign_TS(zr_noltest2$Shapiro_p[i])}
  else{
    zr_noltest2$Shapiro_p[i] <- NaN
    zr_noltest2$Sign[i] <- 'NULL'}
}

write_xlsx(zr_noltest1, 'zr_noltest1.xlsx')
write_xlsx(zr_noltest2, 'zr_noltest2.xlsx')

#3.3 difference test----
#1vs36 (fix nestedness and turnover)
zr_wilcox1 <- data.frame(oZr=rep(c('oZr', 'oZr1', 'oZr4', 'oZr7'), each=36), 
                        qZr=rep(c('qZr', 'qZr1', 'qZr4', 'qZr7'), each=36), 
                        Nest=rep(1:6, each=6), Turn=rep(1:6), 
                        Wilcox_W=NA, Wilcox_p=NA, Sign=NA, oZr_med=NA, qZr_med=NA, delta_med=NA, Winner=NA)
for(i in 1:nrow(zr_wilcox1)){
  i_ozr <- zr_wilcox1$oZr[i]
  i_qzr <- zr_wilcox1$qZr[i]
  i_nest <- zr_wilcox1$Nest[i]
  i_turn <- zr_wilcox1$Turn[i]
  
  data_tp <- filter(zr_summary, Nest==i_nest, Turn==i_turn)[[i_qzr]]
  zr_wilcox1$oZr_med[i] <- filter(oZr_summary, Nest==i_nest, Turn==i_turn)[[i_ozr]]
  zr_wilcox1$qZr_med[i] <- median(data_tp)
  zr_wilcox1$delta_med[i] <- zr_wilcox1$qZr_med[i]-zr_wilcox1$oZr_med[i]
  
  test_tp <- wilcox.test(data_tp, mu=zr_wilcox1$oZr_med[i], alternative='two.sided', exact=F, correct=F)
  zr_wilcox1$Wilcox_W[i] <- test_tp$statistic[["V"]]
  zr_wilcox1$Wilcox_p[i] <- test_tp$p.value
  zr_wilcox1$Sign[i] <- Sign_TS(zr_wilcox1$Wilcox_p[i])
  if(is.nan(zr_wilcox1$Wilcox_p[i]) | zr_wilcox1$Wilcox_p[i] > 0.05){
    zr_wilcox1$Winner[i] <- 'None'}
  else {if(zr_wilcox1$oZr_med[i] > zr_wilcox1$qZr_med[i]){
        zr_wilcox1$Winner[i] <- 'oZr'}
       else{zr_wilcox1$Winner[i] <- 'qZr'}}
}

#9vs9 (fix inter-site and inter-species abundance gradient)
zr_wilcox2 <- data.frame(oZr=rep(c('oZr', 'oZr1', 'oZr4', 'oZr7'), each=9), 
                         qZr=rep(c('qZr', 'qZr1', 'qZr4', 'qZr7'), each=9), 
                         Grad_row=rep(1:3, each=3), Grad_col=rep(1:3), 
                         Wilcox_W=NA, Wilcox_p=NA, Sign=NA, oZr_med=NA, qZr_med=NA, delta_med=NA,Winner=NA)
for(i in 1:nrow(zr_wilcox2)){
  i_ozr <- zr_wilcox2$oZr[i]
  i_qzr <- zr_wilcox2$qZr[i]
  i_grar <- zr_wilcox2$Grad_row[i]
  i_grac <- zr_wilcox2$Grad_col[i]
  
  data_tp1 <- filter(zr_summary, Grad_row==i_grar, Grad_col==i_grac)[[i_ozr]]
  data_tp2 <- filter(zr_summary, Grad_row==i_grar, Grad_col==i_grac)[[i_qzr]]
  zr_wilcox2$oZr_med[i] <- median(data_tp1)
  zr_wilcox2$qZr_med[i] <- median(data_tp2)
  zr_wilcox2$delta_med[i] <- zr_wilcox2$qZr_med[i]-zr_wilcox2$oZr_med[i]
  
  test_tp <- wilcox.test(data_tp1, data_tp2, alternative='two.sided', paired=T,  exact=F, correct=F)
  zr_wilcox2$Wilcox_W[i] <- test_tp$statistic[["V"]]
  zr_wilcox2$Wilcox_p[i] <- test_tp$p.value
  zr_wilcox2$Sign[i] <- Sign_TS(zr_wilcox2$Wilcox_p[i])
  if(is.nan(zr_wilcox2$Wilcox_p[i]) | zr_wilcox2$Wilcox_p[i] > 0.05){
    zr_wilcox2$Winner[i] <- 'None'}
  else {if(zr_wilcox2$oZr_med[i] > zr_wilcox2$qZr_med[i]){
    zr_wilcox2$Winner[i] <- 'oZr'}
    else{zr_wilcox2$Winner[i] <- 'qZr'}}
}

write_xlsx(zr_wilcox1, 'zr_wilcox1.xlsx')
write_xlsx(zr_wilcox2, 'zr_wilcox2.xlsx')

#3.4 generate tables for heat maps----
ht_table2 <- data.frame(matrix(nrow=6, ncol=6))
rownames(ht_table2) <- c('Nest1','Nest2','Nest3','Nest4','Nest5','Nest6')
colnames(ht_table2) <- c('Turn1','Turn2','Turn3','Turn4','Turn5','Turn6')

ht_table3 <- data.frame(matrix(nrow=3, ncol=3))
colnames(ht_table3) <- rownames(ht_table3) <- c('Grad1','Grad2','Grad3')

zr_name <- c('qZr', 'qZr1', 'qZr4', 'qZr7')

#1vs36 (fix nestedness and turnover)
zr_winner_list1 <- zr_table_list1 <- as.vector(rep(NA, 4), mode = 'list')
names(zr_winner_list1) <- names(zr_table_list1) <- zr_name
for (k in 1:4) {
  zr_table_list1[[k]] <- ht_table2
  zr_winner_list1[[k]] <- ht_table2
  for (i in 1:6) {
    for (j in 1:6) {
      pvalue_tp <- filter(zr_wilcox1, qZr==zr_name[k], Nest==i, Turn==j)$delta_med
      sign_tp <- filter(zr_wilcox1, qZr==zr_name[k], Nest==i, Turn==j)$Sign
      zr_table_list1[[k]][i,j] <- pvalue_tp
      zr_winner_list1[[k]][i,j] <- paste0(round(pvalue_tp, 3), '\n', sign_tp)
}}}

#9vs9 (fix abundance gradient)
zr_winner_list2 <- zr_table_list2 <- as.vector(rep(NA, 4), mode = 'list')
names(zr_winner_list2) <- names(zr_table_list2) <- zr_name
for (k in 1:4) {
  zr_table_list2[[k]] <- ht_table3
  zr_winner_list2[[k]] <- ht_table3
  for (i in 1:3) {
    for (j in 1:3) {
      pvalue_tp <- filter(zr_wilcox2, qZr==zr_name[k], Grad_row==i, Grad_col==j)$delta_med
      sign_tp <- filter(zr_wilcox2, qZr==zr_name[k], Grad_row==i, Grad_col==j)$Sign
      zr_table_list2[[k]][i,j] <- pvalue_tp
      zr_winner_list2[[k]][i,j] <- paste0(round(pvalue_tp, 3), '\n', sign_tp)
}}}

#3.5 plot heat maps----
heatmap_color <- c('#507aaf','#fbf9fa','#BE5C37')
p_zr_list2 <- p_zr_list1 <- as.vector(rep(NA, 4), mode = 'list')
names(p_zr_list2) <- names(p_zr_list1) <- zr_name

for (i in 1:4) {
  p_zr_list1[[i]] <- Heatmap_Painter(zr_table_list1[[i]], title=zr_name[i], col_bk=seq(-0.15, 0.15, 0.002), 
                                     cell_color=colorRampPalette(heatmap_color)(151), dpnum=zr_winner_list1[[i]])
  p_zr_list2[[i]] <- Heatmap_Painter(zr_table_list2[[i]], title=zr_name[i], col_bk=seq(-0.15, 0.15, 0.002), 
                                     cell_color=colorRampPalette(heatmap_color)(151), dpnum=zr_winner_list2[[i]])
  ggsave(paste0(zr_name[i],'_1.png'), p_zr_list1[[i]], width = 1500, height = 1400, units = 'px')
  ggsave(paste0(zr_name[i],'_2.png'), p_zr_list2[[i]], width = 1000, height = 800, units = 'px')
}

}