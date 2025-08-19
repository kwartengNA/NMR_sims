# Model Fitting Script
# Setup ############################
options(save.workspace = "no")
num_cores <- 10
gendir <- getwd()
seedval <- 120
N.sim <- 1000

# Custom library path
.libPaths(c("/h/kwarteng/mypkgs", .libPaths()))

# Scenario Definitions ############################
consistency_options <- c(TRUE, FALSE)
assumption_options <- c("common", "independent","NMA")# 
arm_options <- c("mixarm", "multiarm", "twoarm") 
tau_options <- c("none", "moderate", "high") 
density_options <- c("low", "moderate", "all")
event_rate_options <- c("moderate") 
outcome_options <- c("binary") 
interaction_direction_options <- c("positive")

# Libraries and Functions ############################
library(parallel)
library(dplyr)
library(tidyr)
library(jagsUI)
library(netmeta)
source("netmetareg.R") # custom functions for UMIE

# Functions ############################

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
  if (arms == "multiarm" && (density != "all" || !consistency)) return(TRUE)
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


# Tracking Information for Loop ############################
loopnum <- 0
status_frame <- data.frame()
starttime_loop <- Sys.time()


# Model Fitting Loop ############################
for (true_consistency in consistency_options) {
  for (true_assumption in assumption_options) {
    for (true_arms in arm_options) {
      for (tau_txt in tau_options) {
        for (keep_prop_txt in density_options) {
          for (interaction_direction in interaction_direction_options) {
            for (event_rate in event_rate_options) {
              for (outcome in outcome_options) {
                
                # Check Inputs ############################
                if (should_skip(outcome, event_rate, interaction_direction, true_arms, keep_prop_txt, true_consistency, true_assumption)==TRUE){ next}
                keepls <- c(ls(), 
                            "N.arms","keepls")
                
                loopdir <- build_loopdir(outcome, true_arms, true_consistency, true_assumption, tau_txt, keep_prop_txt, interaction_direction, event_rate)
                
                status_list <- list(
                  true_consistency = true_consistency,
                  true_assumption = true_assumption,
                  true_arms = true_arms,
                  tau_txt = tau_txt,
                  keep_prop_txt = keep_prop_txt,
                  outcome = outcome,
                  interaction_direction = interaction_direction,
                  event_rate = ifelse(outcome == "continuous", NA, event_rate))
                
                
                
                # Fit Models  ############################
                if (outcome == "continuous") {
                  tryCatch({
                    refdat<-readRDS(file = paste(loopdir,"/refdat_continuous.rds",sep=""))
                    data2_df<-readRDS(file = paste(loopdir,"/data2_continuous.rds",sep=""))
                    data2<-split(data2_df, data2_df$sim)
                    es_MD<-refdat$Estimate[refdat$Reference=="A" & refdat$paramtype=="es"]
                    beta_x<-refdat$Estimate[refdat$Reference=="A" & refdat$paramtype=="beta_x"]
                    jagsseedval<-seedval
                    N.sim<-ifelse(exists("N.sim")==FALSE,length(data2),N.sim)
                    gc()
                    
                    source("run_models/Run_JAGS_Models_Cont_all.R")
                    gc()
                    source("run_models/Run_Frequentist_Models_Cont2.R")
                    gc()
                    source("run_models/Output_Models.R") 
                    saveRDS(all_res, file = paste(loopdir,"/res_continuous.rds",sep=""), compress="xz")
                    gc()
                    
                  }, error = function(e) {
                    message("Error in continuous sim: ", e$message)
                    status_list$success <- FALSE
                  })
                } else if (outcome == "binary") {
                  tryCatch({
                    refdat<-readRDS(file = paste(loopdir,"/refdat_binary.rds",sep=""))
                    data2_df<-readRDS(file = paste(loopdir,"/data2_binary.rds",sep=""))
                    data2<-split(data2_df, data2_df$sim)
                    rm(data2_df)
                    gc()
                    es_lOR<-refdat$Estimate[refdat$Reference=="A" & refdat$paramtype=="es"]
                    beta_x<-refdat$Estimate[refdat$Reference=="A" & refdat$paramtype=="beta_x"]
                    jagsseedval<-seedval
                    N.sim<-ifelse(exists("N.sim")==FALSE,length(data2),N.sim)
                    data2<-data2[1:N.sim]
                    
                    source("run_models/Run_JAGS_Models_Bin_all.R")#update
                    gc()
                    source("run_models/Run_Frequentist_Models_Bin3.R")# update
                    gc()
                    source("run_models/Output_Models.R")
                    saveRDS(all_res, file = paste(loopdir,"/res_binary.rds",sep=""), compress="xz")
                    gc()
                    status_list["success"]<-TRUE
                  }, error = function(e) {
                    message("Error in binary sim: ", e$message)
                    status_list$success <- FALSE
                  })
                }
                
                on.exit(stopCluster(cl))
                
                # Update status, print to screen, cleanup ############################
                # assign(paste("dat.studies_merge", keep_prop_txt, true_arms, sep="_"), dat.studies_merge)
                loopnum <- loopnum + 1
                status_list$loopnum <- loopnum
                status_list$duration <- Sys.time() - starttime_loop
                status_frame <- bind_rows(status_frame, status_list)
                
                print(loopnum)
                print(loopdir)
                print(Sys.time() - starttime_loop)
                
                rm(list = ls()[!ls() %in% keepls])
                gc()
              }
            }
          }
        }
      }
    }
  }
}

# Summary Check  ############################
endtime_loop <- Sys.time()
total_time <- endtime_loop - starttime_loop

expected <- expand.grid(
  true_consistency = consistency_options,
  true_assumption = assumption_options,
  true_arms = arm_options,
  tau_txt = tau_options,
  keep_prop_txt = density_options,
  outcome = outcome_options,
  interaction_direction = interaction_direction_options,
  event_rate = c(event_rate_options, NA)) %>%
  filter(
    !(outcome == "continuous" & !is.na(event_rate)),
    !(true_arms == "multiarm" & keep_prop_txt != "all"),
    !(outcome == "binary" & is.na(event_rate))) %>%
  mutate(expected = TRUE)

expected %>%
  full_join(status_frame %>%
              mutate(attempted = TRUE,
                     true_consistency = as.logical(true_consistency))) %>%
  count(true_arms, outcome, expected, attempted, success)

