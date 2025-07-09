options(save.workspace = "no")
num_cores<-2
# run simulations
# libs and functions ############################
library(parallel)
library(dplyr)
library(ggplot2)
# library(tibble)
library(tidyr)
library(jagsUI)
library(netmeta)
library(knitr)
knitr::opts_chunk$set(echo = TRUE)
source("netmetareg.R")
# sessionInfo()
is_cluster_active <- function(){tryCatch({
  response <- parLapply(cl, 1:length(cl), function(x) TRUE)
  all(unlist(response))  # Check if all nodes responded
}, error = function(e) FALSE)  # If an error occurs, assume inactive)
}
gendir <- getwd() 
seedval <- 120 # seed value
N.sim <-20
# items to catch within loop ###################
loopnum<-0
status_frame<-as.data.frame(NULL)
starttime_loop <- Sys.time()
# scenario definitions  ############################
consistency_options<-c(TRUE,FALSE)
assumption_options<-c("common", "independent", "NMA")
arm_options<-c("multiarm", "mixarm", "twoarm") #"multiarm", "mixarm", "twoarm"
tau_options<-c("none", "moderate", "high") # #0,0.25,0.5 
density_options<-c( "low","moderate","all") #1,0.5,0.25 
event_rate_options<-c("moderate") #"rare", "moderate", "frequent"
outcome_options<-c("continuous", "binary")#, "binary"
interaction_direction_options<-c("positive") #"positive", "negative", "random"

# start loop  ############################
for (true_consistency in consistency_options) {
  for (true_assumption in assumption_options) {
    for (true_arms in arm_options) {
      for (tau_txt in tau_options) {
        for (keep_prop_txt in density_options) {
          for(interaction_direction in interaction_direction_options){
            for(event_rate in event_rate_options){
              for (outcome in outcome_options) {
                
                # exclude irrelevant combinations
                if (outcome == "continuous" &
                    (event_rate !="moderate"|
                    interaction_direction != "positive")) {
                  next
                  }
                if (true_arms=="multiarm" & (keep_prop_txt!="all"|true_consistency!=TRUE)){
                  next
                }
                if (true_assumption=="NMA" & (true_consistency!=TRUE)){
                  next
                }
                

                keepls <- c(ls(), "N.arms","keepls", paste("keep_studies", keep_prop_txt,true_arms, sep="_"))
                loopdir <- paste(gendir, paste("/seed_", seedval,sep=""),outcome, ifelse(true_consistency == FALSE, "UMIE", "consistency"),
                                 true_assumption,true_arms, paste("tau",tau_txt,sep="_"), paste("connect",keep_prop_txt,sep="_"), 
                                 paste(interaction_direction, "beta", sep="_"), paste(ifelse(outcome=="continuous", "NA", event_rate), "events", sep="_"), sep="/")
                dir.create(loopdir,
                           recursive = TRUE)
                status_list<-c(true_consistency=true_consistency ,true_assumption=true_assumption ,true_arms=true_arms ,tau_txt=tau_txt ,
                               keep_prop_txt=keep_prop_txt ,outcome=outcome, interaction_direction=interaction_direction,  event_rate=ifelse(outcome=="continuous", NA, event_rate))
                tau<-case_when(tau_txt=="none"~0,
                               tau_txt=="moderate"~0.25,
                               tau_txt=="high"~0.5)
                keep_prop<-case_when(keep_prop_txt=="all"~1,
                                     keep_prop_txt=="moderate"~0.5,#sparse. doesn't matter for multiarm.
                                     keep_prop_txt=="low"~0.25)
                
                keep_studies<-NULL
                if(exists(paste("keep_studies", keep_prop_txt, true_arms, sep="_"))){
                  assign("keep_studies",get(paste("keep_studies", keep_prop_txt, true_arms, sep="_")))
                }
                
                if (outcome == "continuous") {
                  tryCatch({
                    loopdir<-ifelse(exists("loopdir")==FALSE,getwd(),loopdir)
                    source("simulate_data/sim_cont_all.R")
                    gc()
                    saveRDS(refdat, file = paste(loopdir,"/refdat_continuous.rds",sep=""), compress="xz")
                    saveRDS(bind_rows(data2, .id = "sim"), file = paste(loopdir,"/data2_continuous.rds",sep=""), compress="xz")
                    gc()
                    source("run_models/Run_JAGS_Models_Cont_all.R")#WARNING HARDCODING
                    gc()
                    source("run_models/Run_Frequentist_Models_Cont2.R") #WARNING HARDCODING
                    gc()
                    source("run_models/Output_Models.R") #WARNING HARDCODING
                    saveRDS(all_res, file = paste(loopdir,"/res_continuous.rds",sep=""), compress="xz")
                    gc()
                    status_list["success"]<-TRUE
                    
                  }, error = function(e) 
                    print("error in sim"),
                    status_list["success"]<-FALSE
                  )
                  
                  } else if (outcome == "binary") {
                        keepls <- c(ls(), "keepls")
                        minr<-case_when(event_rate=="rare"~0.03,
                                        event_rate=="moderate"~0.3,
                                        event_rate=="frequent"~1-0.05)
                        maxr<-case_when(event_rate=="rare"~0.05,
                                        event_rate=="moderate"~0.6,
                                        event_rate=="frequent"~1-0.03)
                        tryCatch({
                          # ifelse(rbinom(n=1,size=1,prob=0.5),error,TRUE)
                          loopdir<-ifelse(exists("loopdir")==FALSE,getwd(),loopdir)
                          source("simulate_data/sim_bin_all.R")#create
                          gc()
                          saveRDS(refdat, file = paste(loopdir,"/refdat_binary.rds",sep=""), compress="xz")
                          saveRDS(bind_rows(data2, .id = "sim"), file = paste(loopdir,"/data2_binary.rds",sep=""), compress="xz")
                          gc()
                          source("run_models/Run_JAGS_Models_Bin_all.R")#update
                          gc()
                          source("run_models/Run_Frequentist_Models_Bin2.R")# update
                          gc()
                          source("run_models/Output_Models.R")
                          saveRDS(all_res, file = paste(loopdir,"/res_binary.rds",sep=""), compress="xz")
                          gc()
                          status_list["success"]<-TRUE
                        
                        }, error = function(e)
                          print("error in sim"),
                          status_list["success"]<-FALSE)
                  }
                
                
    
                    
            # if(is_cluster_active()==TRUE) {
            #   stopCluster(cl)
            # }
                on.exit(stopCluster(cl))
            assign(paste("keep_studies", keep_prop_txt, true_arms, sep="_"), keep_studies)
            loopnum<-loopnum+1
            status_list[["loopnum"]]<-loopnum
            status_list[["duration"]]<-Sys.time() - starttime_loop
            status_frame<-status_frame%>%bind_rows(status_list)
            print(loopnum)
            print(loopdir)
            print(Sys.time() - starttime_loop)
            print(utils::object.size(ls()))
            rm(list = ls()[!ls() %in% keepls])
            gc()
                  }
              } # end event rate
            
          } # end interaction direction
        }
      }
    }
  }
  
}

# brief check ##################################################################
endtime_loop <- Sys.time()
total_time <- endtime_loop - starttime_loop

expected<-data.frame(true_consistency = consistency_options) %>%
  cross_join( data.frame(true_assumption = assumption_options) )%>%
  cross_join( data.frame(true_arms = arm_options) )%>%
  cross_join( data.frame(tau_txt = tau_options) )%>%
  cross_join( data.frame(keep_prop_txt = density_options) )%>%
  cross_join( data.frame(outcome = outcome_options) )%>%
  cross_join(data.frame(interaction_direction = interaction_direction_options))%>%
  cross_join( data.frame(event_rate = c(event_rate_options,NA)))%>%
  # filter(!event_rate %in% c("rare", "frequent"))%>%
  # filter(!interaction_direction %in% c("random", "negative"))%>%
  filter(!(outcome=="continuous" & is.na(event_rate)==FALSE) &
           # !(outcome=="continuous" & interaction_direction!="random")&
           !(true_arms=="multiarm" & keep_prop_txt!="all")&
           # !(true_arms=="multiarm" & true_consistency!=TRUE)&
           !(outcome=="binary" & is.na(event_rate)==TRUE))%>%
  mutate(expected=TRUE)

expected%>%
  full_join(status_frame%>%
              mutate(attempted=TRUE,
                     true_consistency=as.logical(true_consistency)))%>%
  # filter(success ==FALSE)%>%
  count(true_arms,outcome,#keep_prop_txt,
        expected, attempted, success)
