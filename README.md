# NMR_sims

This repository contains code for executing Bayesian modeling simulations using JAGS and related tools.

File Structure

    run_sims: Main script that controls execution. It defines which scenarios will be run in the simulations.

    NMR_simulated_ex_*_generic.R/: It defines which simulation syntax to run. It calls the modeling syntax. It combines results from Bayesian and Frequentist models in a tabular format. If run_sims is not called, it contains its own defaults
	
	sim_*_.R: contains R scripts that generate the simulated networks.
	
	Run_Frequentist_Models_*.R: contain R scripts that execute netmetareg and netmeta frequentist NMA and NMR and transforms the results into tabular format
	
	Run_JAGS_Models_*.R: contain R scripts that execute models in JAGS and transforms the results into tabular format
	
	JAGS/: Contains JAGS model scripts used for Bayesian analysis. These are executed via the Run_JAGS_Models_*.R files.
	
	functions/: Contains NetMeta-related functions used across the modeling process.
  


Usage

Run run_sims to start the full modeling and simulation pipeline for multiple scenarios

Run NMR_simulated_ex_*_generic.R if interested in a specific scenario
