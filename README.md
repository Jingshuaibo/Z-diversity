# Z-diversity
R scripts and raw data for Z diversity

## 1. Data  
The raw data for Z diversity test, including:  
### 1.1 sim_occur_data  
  36 simulated presence-absence matrices named as **"dataxy.xlsx"** (**x**: nestedness gradient; **y**: turnover gradient), with speices as columns and site as rows.
### 1.2 abun(occur)_list.RData  
  Simulated presence-absence/abundance species-site matrices stored in the form of RData file for use in R.
### 1.3 empirical data  
Empirical data for Z diversity test  
**Bird(Butterfly)Data.xlsx**: the species-site abundance matrices of birds and butterflies  
**Bird(Butterfly)_info.xlsx**: the family, species name and short name of each species  
Data source: Myers, M. C., Mason, J. T., Hoksch, B. J., Cambardella, C. A. & Pfrimmer, J. D. (2015). Birds and butterflies respond to soil-induced habitat heterogeneity in experimental plantings of tallgrass prairie species managed as agroenergy crops in Iowa, USA. J. Appl. Ecol. 52, 1176-1187.  

## 2. Rscripts  
  The R scripts for Z diversity test, including:  
### 2.1 Zdiv_Function.R  
Functions for the calculation of Z diversity and Z ratios and model fitting of Z decline forms.  
#### Zdiv_cal(x)  
The function for Z diversity analysis. The input x is an abundance matrix with species as columns and sites as rows, and the output is the results as a list data, including the qF components, Z diversity components, Z ratios and model fitting results for the community and each species.  
**Note**: this program needs R package *dplyr* to run
### 2.2 simData_produce.R
The generation and process for the simulated presence-absence/abundance matrices.
### 2.3 beta_cal.R
The calculation of beta diversity metrices for the simulated presence-absence/abundance matrices
### 2.4 zdiv_simData.R
Z diversity analysis on the simulated presence-absence/abundance matrices
### 2.5 empirical_data.R  
Z diversity analysis on the empirical matrices  

**For any questions, contact me via email: 1747225066@qq.com**

 
