# NMR_sims

This repository contains code for executing Bayesian modeling simulations using JAGS and related tools.

File Structure

    run_sims: Main script that controls execution. It defines which scenarios will be run in the simulations.

	simulate_data/ contains R scripts that generate the simulated networks.
	
	run_models/ contain R scripts that execute netmetareg and netmeta frequentist NMA and NMR, Bayesian models in JAGS and transforms the results into tabular format 
	
	JAGS/: Contains JAGS model scripts used for Bayesian analysis. These are executed via the Run_JAGS_Models_*.R files.
	
	functions/: Contains NetMeta-related functions used across the modeling process.
  


Usage

Run run_sims to start the full modeling and simulation pipeline