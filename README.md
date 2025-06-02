# Introduction
This repository houses the R code used for simulation studies within the project titled “An augmented pseudo-likelihood approach for population-level modeling of HPV infection and cervical precancer natural history using current-status data.” Inside, there are two folders. One folder labeled utility contains all utility functions; the other folder labeled sim contains the running scripts designed for executing simulations.

# utility
This folder contains:
	- dt_gen_fct.R and dt_csd_gen.R contain utility functions for data generation;and
	- commonf.cpp, commonf.R, est_rec20_g_j.R, est.R, est_rec20.R contain functions for estimation and inference.

# sim
This folder contains:
	- sim_g_j.R to execute an example simulation run with a current-status data on complete state observation (Data type 1:  D_1);
	- sim_g_S.R to execute an example simulation run with a current-status data on HPV infection with unknown cervical precancer status (Data type 2: D_2);
	- sim.R to execute an example simulation run with a current-status data on HPV infection with missing sexual debut age (Data type 3:  D_3); and
	- sim_g_S_boot.R to execute an example simulation run with a bootstrapped D_2.

Each script is designed to execute an example simulation study. By modifying the data configurations, it can replicate all the simulation results outlined in the manuscript. The script employs a for-loop to conduct repeated simulations. Within each simulation iteration, three key steps are executed: 1) data generation; 2) estimation and inference; and 3) saving the results.
