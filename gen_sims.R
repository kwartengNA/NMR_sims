# Simulation Generating Script
# Setup ############################
options(save.workspace = "no")
num_cores <- 2
gendir <- getwd()
seedval <- 389
N.sim <- 1000

# Custom library path
.libPaths(c("/h/kwarteng/mypkgs", .libPaths()))

# Scenario Definitions ############################
consistency_options <- c(TRUE, FALSE)
assumption_options <- c("common", "independent","NMA")# 
arm_options <- c("mixarm", "multiarm", "twoarm") 
tau_options <- c("none", "moderate", "high") 
density_options <- c("low", "moderate_dense","moderate_sparse", "full")
event_rate_options <- c("moderate") 
outcome_options <- c("binary") 
interaction_direction_options <- c("positive")

# Libraries and Functions ############################
# library(parallel)
library(dplyr)
library(tidyr)
# library(jagsUI)
# library(netmeta)
# source("netmetareg.R") # custom functions for UMIE

# Functions ############################
# Convert tau text label to numeric value
get_tau <- function(label) {
  dplyr::case_when(
    label == "none"     ~ 0,
    label == "moderate" ~ 0.28, #Röver
    label == "high"     ~ 0.5
  )
}

# Convert density label to numeric proportion
get_keep_prop <- function(label) {
  dplyr::case_when(
    label == "full"      ~ 1,
    label == "moderate_dense" ~ 0.75,
    label == "moderate_sparse" ~ 0.5,
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
  if (arms == "multiarm" && (density != "full" || !consistency)) return(TRUE)
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

# Tracking Information for Loop ############################
loopnum <- 0
status_frame <- data.frame()
starttime_loop <- Sys.time()


# Data Generation Loop ############################
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
                            "N.arms","keepls", paste("dat.studies_merge", keep_prop_txt,true_arms, sep="_"))
                
                loopdir <- build_loopdir(outcome, true_arms, true_consistency, true_assumption, tau_txt, keep_prop_txt, interaction_direction, event_rate)
                dir.create(loopdir, recursive = TRUE, showWarnings = FALSE)
                
                status_list <- list(
                  true_consistency = true_consistency,
                  true_assumption = true_assumption,
                  true_arms = true_arms,
                  tau_txt = tau_txt,
                  keep_prop_txt = keep_prop_txt,
                  outcome = outcome,
                  interaction_direction = interaction_direction,
                  event_rate = ifelse(outcome == "continuous", NA, event_rate))
                
                tau <- get_tau(tau_txt)
                keep_prop <- get_keep_prop(keep_prop_txt)
                
                
                # retain structure for identical networks for each density and arm structure
                if(exists(paste("dat.studies_merge", keep_prop_txt, true_arms, sep="_"))){
                  assign("dat.studies_merge",get(paste("dat.studies_merge", keep_prop_txt, true_arms, sep="_")))
                }
                
                # Simulate Data ############################
                if (outcome == "continuous") {
                  tryCatch({
                    source("simulate_data/sim_cont_all.R")
                    saveRDS(refdat, file = file.path(loopdir, "refdat_continuous.rds"), compress = "xz")
                    saveRDS(bind_rows(data2, .id = "sim"), file = file.path(loopdir, "data2_continuous.rds"), compress = "xz")
                    status_list$success <- TRUE
                  }, error = function(e) {
                    message("Error in continuous sim: ", e$message)
                    status_list$success <- FALSE
                  })
                } else if (outcome == "binary") {
                  tryCatch({
                    source("simulate_data/sim_bin_all.R")
                    saveRDS(refdat, file = file.path(loopdir, "refdat_binary.rds"), compress = "xz")
                    saveRDS(bind_rows(data2, .id = "sim"), file = file.path(loopdir, "data2_binary.rds"), compress = "xz")
                    status_list$success <- TRUE
                  }, error = function(e) {
                    message("Error in binary sim: ", e$message)
                    status_list$success <- FALSE
                  })
                }
                
                # on.exit(stopCluster(cl))
                
                # Update status, print to screen, cleanup ############################
                assign(paste("dat.studies_merge", keep_prop_txt, true_arms, sep="_"), dat.studies_merge)
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
    !(true_arms == "multiarm" & keep_prop_txt != "full"),
    !(outcome == "binary" & is.na(event_rate))) %>%
  mutate(expected = TRUE)

expected %>%
  full_join(status_frame %>%
              mutate(attempted = TRUE,
                     true_consistency = as.logical(true_consistency))) %>%
  count(true_arms, outcome, expected, attempted, success)

