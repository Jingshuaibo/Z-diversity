library(betapart)
library(dplyr)
library(readxl)
library(writexl)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(ggplotify)

#0.Prepare----
#读取数据
load("abun_list.RData")
load("occur_list.RData")

#function for heatmap
Heatmap_Painter <- function(table_data, title=NULL, ant_r=NULL, ant_c=NULL, ant_color=NULL, col_bk=NULL, lg_bk=NULL, 
                            gap_r=NULL, gap_c=NULL, cell_color=NULL, show_legend=T, show_r=T, show_c=T){
  p_heatmap <- pheatmap(table_data, main = title, fontsize = 11, 
                        show_rownames = show_r, show_colnames = show_c, angle_col=0, 
                        display_numbers = T, number_format = "%.2f", number_color = "black", fontsize_number=12, 
                        color = colorRampPalette(cell_color)(201), breaks = col_bk, 
                        cellheight = 40, cellwidth = 40, border_color = "black", 
                        scale = "none", cluster_rows = F, cluster_cols = F, 
                        gaps_row = gap_r, gaps_col = gap_c, 
                        annotation_col = ant_c, annotation_row = ant_r, 
                        annotation_colors = ant_color, 
                        legend = show_legend, legend_breaks = lg_bk, legend_labels = lg_bk)
  return(p_heatmap)
}

#1.Occurrence_data----
{#1.1 beta diversity calculation----
beta_summary <- data.frame('Nest'=rep(1:6, each=6), 'Turn'=rep(1:6, 6), 
                           'SNE'=NA, 'SIM'=NA, 'SOR'=NA, 'SNE_ratio'=NA, 'SIM_ratio'=NA, 'TNR'=NA)
for (x in 1:6) {
  for (y in 1:6) {
    beta_tp <- beta.multi(data_list[[x]][[y]], index.family = "sorensen")
    beta_summary[which(beta_summary$Nest==x & beta_summary$Turn==y), 3:5] <- beta_tp[c(2,1,3)]
  }
}
beta_summary$SNE_ratio <- beta_summary$SNE/beta_summary$SOR
beta_summary$SIM_ratio <- beta_summary$SIM/beta_summary$SOR
beta_summary$TNR <- beta_summary$SIM_ratio-beta_summary$SNE_ratio
beta_summary[is.na(beta_summary)] <- 0

write_xlsx(beta_summary, 'beta_summary.xlsx')

#1.2 generate tables for heatmap----
ht_table <- data.frame(matrix(ncol=6, nrow=6))
rownames(ht_table) <- paste0('Nest', c(1:6))
colnames(ht_table) <- paste0('Turn', c(1:6))
beta_name <- c('SOR', 'SIM', 'SNE', 'TNR')

beta_table_list <- as.vector(rep(NA, 4), mode = 'list')
names(beta_table_list) <- beta_name
for (k in 1:4) {
  beta_table_list[[k]] <- ht_table
  for (i in 1:6) {
    beta_table_list[[k]][i] <- filter(beta_summary, Turn==i)[beta_name[k]]
}}

#1.3 plot heatmap----
ht_color <- c('#507aaf','#fbf9fa','#bd3d3f')

p_beta_list <- as.vector(rep(NA, 4), mode = 'list')
names(p_beta_list) <- beta_name
for (i in 1:4) {
  p_beta_list[[i]] <- Heatmap_Painter(beta_table_list[[i]], title=beta_name[[i]], 
                                      col_bk=seq(-1,1,0.01), lg_bk=seq(-1,1,0.5), cell_color=ht_color)
  ggsave(paste0(beta_name[i],'.png'), p_beta_list[[i]], width = 1500, height = 1400, units = 'px')
}

}

