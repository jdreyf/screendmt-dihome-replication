# screendmt-dihome-replication
Data, analysis code, and results from Dreyfuss et al. (2024).

## Download
1. Click on the green icon "Clone or download" then "Download ZIP"
2. Unzip

## Prerequisites
1. Download and Install R: https://cran.r-project.org/
2. Download and Install RStudio Desktop: https://rstudio.com/products/rstudio/download/
3. On Windows, you should have [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

## Lipidomics analysis
1. Open the main document here, `analyze_screendmt_dihome_replication.Rmd`, in R Studio
2. Install and load necessary R packages, including ours such as `DirectionalMaxPTest`.
	+ If asked "Do you want to install from sources the package which needs compilation?" say "no".
3. You can execute the R code in this file in blocks, or you can execute the entire file, including text, to produce an HTML file of the same name using the R Studio button "Knit". This will load our data from the `data` folder. It then applies our and other packages to  our data to reproduce the analyses from the paper.
4. You can compare the results to those at https://github.com/jdreyf/screendmt-dihome-replication/tree/main/results.

## Simulations
1. Open `simulations_screendmt.Rmd` in R Studio
2. Install and load necessary R packages, including ours such as `DirectionalMaxPTest`.
	+ If asked "Do you want to install from sources the package which needs compilation?" say "no".
3. You can execute the R code in this file in blocks, or you can execute the entire file, including text, to produce an HTML file of the same name using the R Studio button "Knit". However, this is very slow.
4. You can see the simulation results for FDR by opening `pwr_fdr_arr_b1000_m1000_mu3.5_logpi1.RDS` or for FWER by opening `pwr_fwer_arr_b1000_m1000_mu3.5_logpi1.RDS` in R studio. These are multidimensional arrays that store all simulation results. These reults are plotting in the simulations_screendmt.Rmd chunk `ggplot` in the `Plots` section.
