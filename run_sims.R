# Simulation Script. Calls Generation, Model fit, and creates summary statistics.
# In practice segments may need to be called separately to accommodate memory
### Setup ############################
options(save.workspace = "no")
num_cores <- 2
gendir <- getwd()
seedval <- 545
N.sim <- 1000
NT <- 7   # number of treatments
NS <- 4   # number of studies per comparison (NS*NT for multi-arm)

# Custom library path
.libPaths(c("/h/kwarteng/mypkgs", .libPaths()))
# 
### Scenario Definitions ############################
consistency_options <- c(TRUE, FALSE)
assumption_options <- c("common", "independent","NMA")#
arm_options <- c("mixarm", "multiarm", "twoarm")
tau_options <- c("none", "moderate","moderate3") # none moderate moderate3
density_options <- c("low", "moderate_dense","moderate","moderate_sparse","full")# "moderate_dense","moderate","moderate_sparse","full"
event_rate_options <- c("moderate")
outcome_options <- c("binary")
interaction_direction_options <- c("positive")

# Libraries and Functions ############################
library(parallel)
library(dplyr)
library(dtplyr)
library(tidyr)
library(jagsUI)
library(netmeta)
library(R.utils)
library(mcstatsim)
source("netmetareg.R") # custom functions for UMIE

# Functions ############################
# Convert tau text label to numeric value
get_tau <- function(label) {
  dplyr::case_when(
    label == "none"     ~ 0,
    # label == "moderate" ~ 0.25,
    # label == "moderate" ~0.24, #doi:10.1093/ije/dys041. Semi Objective Non Pharmacological
    # label == "moderate2" ~0.36, #doi:10.1093/ije/dys041. Subjective Non Pharmacological
    label == "moderate" ~ 0.28, # https://onlinelibrary.wiley.com/doi/epdf/10.1002/sim.6381 # median rlnorm(10000,-2.56,1.74)
    label == "moderate3" ~0.5#, # 3rd quartile rlnorm(10000,-2.56,1.74)
    # label == "high"     ~ 1 # Seide
  )
}

# Convert density label to numeric proportion
get_keep_prop <- function(label) {
  dplyr::case_when(
    label == "full"      ~ 1,
    label == "all"      ~ 1,
    label == "moderate_dense" ~ 0.75,
    label == "moderate" ~ 0.5,
    label == "moderate_sparse" ~ 0.35,
    label == "low"      ~ 0.25
  )
}

# Event rate ranges
get_event_range <- function(rate) {
  minr <- dplyr::case_when(
    rate == "rare"     ~ 0.03,
    rate == "moderate" ~ 0.3,
    rate == "frequent" ~ 0.95
  )
  maxr <- dplyr::case_when(
    rate == "rare"     ~ 0.05,
    rate == "moderate" ~ 0.6,
    rate == "frequent" ~ 0.97
  )
  list(min = minr, max = maxr)
}

# Build directory path for results
build_loopdir <- function(outcome, arms, consistency, assumption, tau, density, interaction_direction, event_rate) {
  paste(
    gendir, paste0("/seed_", seedval), outcome, arms,
    ifelse(consistency, "consistency", "UMIE"),
    assumption, paste0("tau_", tau), paste0("connect_", density),
    paste0(interaction_direction, "_beta"),
    paste0(ifelse(outcome == "continuous", "NA", event_rate), "_events"),
    sep = "/"
  )
}

# Exclusion logic
should_skip <- function(outcome, event_rate, interaction_direction, arms, density, consistency, assumption) {
  if (outcome == "continuous" && (event_rate != "moderate" || interaction_direction != "positive")) return(TRUE)
  if (arms == "multiarm" && (! density %in% c("full","all") || !consistency)) return(TRUE)
  if (assumption == "NMA" && !consistency) return(TRUE)
  FALSE
}
# Check if cluster is active
is_cluster_active <- function() {
  tryCatch({
    response <- parLapply(cl, 1:length(cl), function(x) TRUE)
    all(unlist(response))
  }, error = function(e) FALSE)
}

# Calculate number of mixarm studies if relevant
# 10% comparisons from three-arm studies 
# calculate multiarm based on NT. satisfies 
# .132=nmulti/(nmulti+n2arm)
# n2arm=(3 * choose(NT, 2))-(nmulti*3)

ctmixarm <- function(NT,NS) {
  nmulti <- 0.132*NS * choose(NT, 2) / (1-0.132+(0.132*NS))
  n2arm <- (NS * choose(NT, 2)) - (nmulti * NS)
  list(nmulti = nmulti, n2arm = n2arm)
}


keepsettings<-c(ls(),"keepsettings")
# generate simulated data ##########
source("gen_sims.R")
source("simulate_data/combine_simulated_data.R")
gc()
keepsettings<-c(keepsettings,"all_refdat", "alldata2")
rm(list=ls()[!ls()%in%keepsettings])#reset workspace for memory management
gc()
### if separating data generation and modelling, import generated data ##########
# all_refdat <- readRDS(paste(gendir, paste("seed_", seedval,sep=""),"all_refdat.rds", sep="/"))
# alldata2 <- readRDS(paste(gendir, paste("seed_", seedval,sep=""),"all_data.rds", sep="/"))
# run models ##########
source("run_models.R")
rm(list=ls()[!ls()%in%keepsettings])
gc()

# extract and combine results ###############
source("explore_results/check_sims.R")
rm(list=ls()[!ls()%in%keepsettings])
gc()
# process results for shiny app #############
source("explore_results/shiny_transformations.R")
rm(list=ls()[!ls()%in%keepsettings])
gc()