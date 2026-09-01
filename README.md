# Z-diversity
R scripts and raw data for Z diversity

##1.data  
The raw data for Z diversity test, including:  

##sim_occur_data  
  36 simulated presence-absence matrices named as "dataxy.xlsx" (x: nestedness gradient; y: turnover gradient), with speices as columns and site as rows   
  
2)abun(occur)_list.RData  
  simulated abundance and presence-absence species-site matrices stored in the form of RData file for use in R  

3)empirical data  
Bird(Butterfly)Data.xlsx: the species-site abundance matrices of birds and butterflies  
Bird(Butterfly)_info.xlsx: the family, species name and short name of each species  
  
2.Rscripts  
  The R scripts for Z diversity test, including:  
1)Zdiv_Function.R: functions for the calculation of Z diversity and Z ratios and model fitting of Z decline forms.  
*Zdiv_cal(x): The function for Z diversity analysis. The input x is an abundance matrix with species as columns and sites as rows, and the output is the results as a list data, including:
  data_order: the abundance matrix after reordering the abundances of each species across sites in descending order
#data_qF: qF components for the community and each species
#data_zdiv: Z diversity components for the community and each species
#data_zr: Z ratios for the community and each species
#data_md: model fitting results for the community and each species 

2)
 
