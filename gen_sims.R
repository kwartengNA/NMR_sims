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

