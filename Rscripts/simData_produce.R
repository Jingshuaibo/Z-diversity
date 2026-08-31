library(betapart)
library(dplyr)
library(readxl)
library(writexl)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

#0.Import_Data----
data_list1 <- vector(mode = 'list', length = 6)
data_list <- data_list1
for (i in 1:6) {
  data_list[[i]] <- data_list1
}
#Nest=1, Turn=1-6
data_list[[1]][[1]] <- read_excel("sim_occur_data/data11.xlsx")
data_list[[1]][[2]] <- read_excel("sim_occur_data/data12.xlsx")
data_list[[1]][[3]] <- read_excel("sim_occur_data/data13.xlsx")
data_list[[1]][[4]] <- read_excel("sim_occur_data/data14.xlsx")
data_list[[1]][[5]] <- read_excel("sim_occur_data/data15.xlsx")
data_list[[1]][[6]] <- read_excel("sim_occur_data/data16.xlsx")
#Nest=2, Turn=1-6
data_list[[2]][[1]] <- read_excel("sim_occur_data/data21.xlsx")
data_list[[2]][[2]] <- read_excel("sim_occur_data/data22.xlsx")
data_list[[2]][[3]] <- read_excel("sim_occur_data/data23.xlsx")
data_list[[2]][[4]] <- read_excel("sim_occur_data/data24.xlsx")
data_list[[2]][[5]] <- read_excel("sim_occur_data/data25.xlsx")
data_list[[2]][[6]] <- read_excel("sim_occur_data/data26.xlsx")
#Nest=3, Turn=1-6
data_list[[3]][[1]] <- read_excel("sim_occur_data/data31.xlsx")
data_list[[3]][[2]] <- read_excel("sim_occur_data/data32.xlsx")
data_list[[3]][[3]] <- read_excel("sim_occur_data/data33.xlsx")
data_list[[3]][[4]] <- read_excel("sim_occur_data/data34.xlsx")
data_list[[3]][[5]] <- read_excel("sim_occur_data/data35.xlsx")
data_list[[3]][[6]] <- read_excel("sim_occur_data/data36.xlsx")
#Nest=4, Turn=1-6
data_list[[4]][[1]] <- read_excel("sim_occur_data/data41.xlsx")
data_list[[4]][[2]] <- read_excel("sim_occur_data/data42.xlsx")
data_list[[4]][[3]] <- read_excel("sim_occur_data/data43.xlsx")
data_list[[4]][[4]] <- read_excel("sim_occur_data/data44.xlsx")
data_list[[4]][[5]] <- read_excel("sim_occur_data/data45.xlsx")
data_list[[4]][[6]] <- read_excel("sim_occur_data/data46.xlsx")
#Nest=5, Turn=1-6
data_list[[5]][[1]] <- read_excel("sim_occur_data/data51.xlsx")
data_list[[5]][[2]] <- read_excel("sim_occur_data/data52.xlsx")
data_list[[5]][[3]] <- read_excel("sim_occur_data/data53.xlsx")
data_list[[5]][[4]] <- read_excel("sim_occur_data/data54.xlsx")
data_list[[5]][[5]] <- read_excel("sim_occur_data/data55.xlsx")
data_list[[5]][[6]] <- read_excel("sim_occur_data/data56.xlsx")
#Nest=6, Turn=1-6
data_list[[6]][[1]] <- read_excel("sim_occur_data/data61.xlsx")
data_list[[6]][[2]] <- read_excel("sim_occur_data/data62.xlsx")
data_list[[6]][[3]] <- read_excel("sim_occur_data/data63.xlsx")
data_list[[6]][[4]] <- read_excel("sim_occur_data/data64.xlsx")
data_list[[6]][[5]] <- read_excel("sim_occur_data/data65.xlsx")
data_list[[6]][[6]] <- read_excel("sim_occur_data/data66.xlsx")

save(data_list, file='occur_list.RData')

#1.Occurrence_data----
occur_table <- data.frame(matrix(nrow = 60, ncol = 60))
colnames(occur_table) <- paste0('T', rep(1:6, each=10), 'Sp', rep(1:10, 6))
rownames(occur_table) <- paste0('N', rep(1:6, each=10), 'S', rep(1:10, 6))
for (i in 1:6) {
  for (j in 1:6) {
    i_begin <- (i-1)*10+1
    j_begin <- (j-1)*10+1
    occur_table[i_begin:(i_begin+9), j_begin:(j_begin+9)] <- data_list[[i]][[j]]
  }
}

#plot occurrence matrices
p_occurdata <- pheatmap(occur_table, main = 'Occurence data', fontsize = 11, 
                        show_rownames = T, show_colnames = T, angle_col=90,   
                        display_numbers = F, number_format = "%d", number_color = "black", fontsize_number=10, 
                        color = c('white', 'red'), 
                        cellheight = 25, cellwidth = 25, border_color = "black", 
                        scale = "none", cluster_rows = F, cluster_cols = F, 
                        gaps_row = seq(10,60,10), gaps_col = seq(10,60,10), 
                        legend = T, legend_breaks = c(0,1), legend_labels = c('Absense', 'Occurence'))
ggsave('p_occurdata.png', p_occurdata, width = 7000, height = 7000, units = 'px')

#3.Abundance_data----
#generate simulated abundance matrices
grad <- data.frame(grad1=rep(100, 10), 
                   grad2=c(150, 140, 130, 120, 110, 90, 80, 70, 60, 50),
                   grad3=c(190, 170, 150, 130, 110, 90, 70, 50, 30, 10))

abun_list1 <- vector(mode = 'list', length = 3)
abun_list <- abun_list1

for (i in 1:3) {
  abun_list[[i]] <- abun_list1
}
for (x in 1:3) {
  for(y in 1:3){
    grad_row <- grad[x]
    grad_col <- grad[y]
    data_tp <- data_list
    for (i in 1:6) {
      for (j in 1:6) {
        data_tp[[i]][[j]] <- sweep(data_tp[[i]][[j]], MARGIN=1, grad[[x]], '*')
        data_tp[[i]][[j]] <- sweep(data_tp[[i]][[j]], MARGIN=2, grad[[y]], '*')
      }
      abun_list[[x]][[y]] <- data_tp
    }}}

save(abun_list, file='abun_list2.RData')