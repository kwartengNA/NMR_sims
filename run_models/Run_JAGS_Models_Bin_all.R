# Bayesian Binary outcome

# function to check clusters ########
# Check if all nodes are responsive (defined in run_sims. here in case it needs to run standalone)
# is_cluster_active <- function(){tryCatch({
#   response <- parLapply(cl, 1:length(cl), function(x) TRUE)
#   all(unlist(response))  # Check if all nodes responded
# }, error = function(e) FALSE)  # If an error occurs, assume inactive)
#   }

# function to start clusters ################
startjagscl<-function(){
  if(exists("cl")==FALSE|is_cluster_active()==FALSE){
  cl <- makeCluster(num_cores, outfile="temp.txt")
  clusterExport(cl,c("N.sim","data2","es_lOR","beta_x","refdat",#"all_JAGS_lists", 
                     "model_AD_NMA_BIN_RE",
                     "model_AD_NMR_BIN_RE_const","model_AD_NMR_BIN_RE_const_UME",
                     "model_AD_NMR_BIN_RE_independent","model_AD_NMR_BIN_RE_independent_UME",
                     "build_data_filter",
                     "n.iter","n.burnin","seedval", "jagsseedval"))
  clusterEvalQ(cl, {.libPaths( c( "/h/kwarteng/mypkgs" , .libPaths() ) )}) # needed for cluster
  clusterEvalQ(cl, {library(dplyr)})
  clusterEvalQ(cl, {library(jagsUI)})
  clusterEvalQ(cl, {set.seed(jagsseedval)})# if the JAGS models are run separately, they require the same seed in order to obtain the same results
  }
  return(cl)
  }

# Model files iterations, and burnin ####################
model_AD_NMR_BIN_RE_independent<-"JAGS/model_AD_NMR_BIN_RE_independent.txt"
model_AD_NMR_BIN_RE_const<-"JAGS/model_AD_NMR_BIN_RE_constant.txt"
model_AD_NMA_BIN_RE<-"JAGS/model_AD_NMA_BIN_RE.txt"
model_AD_NMR_BIN_RE_const_UME<-"JAGS/model_AD_NMR_BIN_RE_constant_inconsistency_Idir.txt"
model_AD_NMR_BIN_RE_independent_UME<-"JAGS/model_AD_NMR_BIN_RE_independent_inconsistency_Idir.txt"
n.iter = 30000
n.burnin = 10000
# JAGS data prep ###################

build_data_filter <- function(simid){
  temp_mydata_wide<-
    data2[[simid]]%>%
    rename(treat_num=treat)%>%
    mutate(study=study,
           # Ivalrand=ifelse(Ival!=0 & Ivalrand==0,-1,Ivalrand)
    )%>% # for JAGS it needs to be 1 and -1, unless neither ref nor comp is favored.
    arrange(study,treat_num)%>%
    group_by(study)%>%
    mutate(nstudy=sum(n),
           treat_order=row_number())%>%
    ungroup()%>%
    tidyr::pivot_wider(names_from = treat_order,
                       values_from =c("events","n","treat_num","treatment","r",starts_with("Ival_")))#"Ival_std", "Ival_gen"
  
  data_filter <- list("ns"=length(unique(temp_mydata_wide$study)),
                      "nt"=max(temp_mydata_wide%>%
                                 select(starts_with("treat_num")),na.rm=TRUE),
                      "na"=rowSums(is.na(select(temp_mydata_wide,
                                                starts_with("treat_num")))==FALSE),
                      "x"=temp_mydata_wide$xcov,
                      "mx"=0,#mean(temp_mydata_wide$xcov),
                      ## need a column for each arm
                      "t"=temp_mydata_wide%>%select(starts_with("treat_num")),
                      ## need a column for each arm
                      "r"=temp_mydata_wide%>%select(starts_with("r_")),
                      "n"=temp_mydata_wide%>%select(starts_with("n_")),
                      ## need a column for each arm
                      "I_dir"= temp_mydata_wide%>%select(starts_with("Ival_order_"))
  )
  
  return(data_filter)
}

# Define functions for models/output ####################

## NMA ################
run_NMA_JAGS<-function(simid){ # each individual dataset
  data_filter<-build_data_filter(simid = simid)
  
  tryCatch({
    temp_nma<- jags(
      data = data_filter,
      inits = NULL,
      parameters.to.save = c("lor","sd"),
      n.chains = 2, n.iter = n.iter,
      n.burnin = n.burnin,DIC=F,
      model.file = model_AD_NMA_BIN_RE)
    
    NMA_res<-as.data.frame(temp_nma$summary)%>%
      select(TE_NMR=mean,seTE_NMR=sd,
             lci_NMR=`2.5%`,
             uci_NMR=`97.5%`,
             Rhat)%>%
      tibble::rownames_to_column()%>%
      mutate(paramtype=gsub("\\[.*","",rowname),
             tempname=gsub(".*\\[|\\],*","",rowname))%>%
      tidyr::separate(tempname, into=c("comp_num","ref_num"),sep = ",")%>%
      mutate(Comparator=LETTERS[as.numeric(comp_num)],
             Reference =LETTERS[as.numeric(ref_num)],
             paramtype=case_when(paramtype=="lor"~"es",
                                 paramtype=="betadiff"~"beta_x",
                                 paramtype=="sd"~"tau",
                                 TRUE~paramtype))%>%
      left_join(refdat,by=c("Comparator", "Reference","paramtype"))
    
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMA_res")==FALSE){
    NMA_res<-data.frame(error=TRUE)}
  
  )
  rm(temp_nma)
  gc()
  return(NMA_res)
}



## IND ################
run_ind_JAGS<-function(simid){ # each individual dataset
  data_filter<-build_data_filter(simid = simid)
  tryCatch({
    temp_nmr<- jags(
      data = data_filter,
      inits = NULL,
      parameters.to.save = c("lor", "lor0", "betadiff","sd"),
      n.chains = 2, n.iter = n.iter,
      n.burnin = n.burnin,DIC=F,
      model.file = model_AD_NMR_BIN_RE_independent)
    
    NMR_res<-as.data.frame(temp_nmr$summary)%>%
      select(TE_NMR=mean,seTE_NMR=sd,
             lci_NMR=`2.5%`,
             uci_NMR=`97.5%`,
             Rhat)%>%
      tibble::rownames_to_column()%>%
      mutate(paramtype=gsub("\\[.*","",rowname),
             tempname=gsub(".*\\[|\\],*","",rowname))%>%
      tidyr::separate(tempname, into=c("comp_num","ref_num"),sep = ",")%>%
      mutate(Comparator=LETTERS[as.numeric(comp_num)],
             Reference =LETTERS[as.numeric(ref_num)],
             paramtype=case_when(paramtype=="lor0"~"es",
                                 paramtype=="betadiff"~"beta_x",
                                 paramtype=="sd"~"tau",
                                 TRUE~paramtype))%>%
      left_join(refdat,by=c("Comparator", "Reference","paramtype"))
    
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  rm(temp_nmr)
  gc()
  return(NMR_res)
  
}

## COMMON ########################
run_common_JAGS<-function(simid){ # each individual dataset
  data_filter<-build_data_filter(simid = simid)
  
  tryCatch({
    temp_nmr<- jags(
      data = data_filter,
      inits = NULL,
      parameters.to.save = c("lor", "lor0","betadiff","sd"),
      n.chains = 2, n.iter = n.iter,
      n.burnin = n.burnin,DIC=F,
      model.file = model_AD_NMR_BIN_RE_const)
    
    NMR_res<-as.data.frame(temp_nmr$summary)%>%
      select(TE_NMR=mean,seTE_NMR=sd,
             lci_NMR=`2.5%`,
             uci_NMR=`97.5%`,
             Rhat)%>%
      tibble::rownames_to_column()%>%
      mutate(paramtype=gsub("\\[.*","",rowname),
             tempname=gsub(".*\\[|\\],*","",rowname))%>%
      tidyr::separate(tempname, into=c("comp_num","ref_num"),sep = ",")%>%
      mutate(Comparator=LETTERS[as.numeric(comp_num)],
             Reference =LETTERS[as.numeric(ref_num)],
             paramtype=case_when(paramtype=="lor0"~"es",
                                 paramtype=="betadiff"~"beta_x",
                                 paramtype=="sd"~"tau",
                                 TRUE~paramtype))%>%
      left_join(refdat,by=c("Comparator", "Reference","paramtype"))
    
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  rm(temp_nmr)
  gc()
  return(NMR_res)
  
}

## COMMON UMIE ##############

run_common_UMIE_JAGS<-function(simid){ # each individual dataset
  data_filter<-build_data_filter(simid = simid)
  tryCatch({
    temp_nmr<- jags(
      data = data_filter,
      inits = NULL,
      parameters.to.save = c("LOR","lor0", "beta","tau"),
      n.chains = 2, n.iter = n.iter,
      n.burnin = n.burnin,DIC=F,
      model.file = model_AD_NMR_BIN_RE_const_UME)
    
    NMR_res<-as.data.frame(temp_nmr$summary)%>%
      select(TE_NMR=mean,seTE_NMR=sd,
             lci_NMR=`2.5%`,
             uci_NMR=`97.5%`,
             Rhat)%>%
      tibble::rownames_to_column()%>%
      mutate(paramtype=gsub("\\[.*","",rowname),
             tempname=gsub(".*\\[|\\],*","",rowname))%>%
      tidyr::separate(tempname, into=c("comp_num","ref_num"),sep = ",")%>%
      mutate(Comparator=LETTERS[as.numeric(comp_num)],
             Reference =LETTERS[as.numeric(ref_num)],
             paramtype=case_when(paramtype=="lor0"~"es",
                                 paramtype=="beta"~"beta_x",
                                 TRUE~paramtype))%>%
      left_join(refdat,by=c("Comparator", "Reference","paramtype"))
    
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  rm(temp_nmr)
  gc()
  return(NMR_res)
  
}

## IND UMIE #########################
run_ind_UMIE_JAGS<-function(simid){ # each individual dataset
  data_filter<-build_data_filter(simid = simid)
  tryCatch({
    temp_nmr<- jags(
      data = data_filter,
      inits = NULL,
      parameters.to.save = c("LOR","lor0", "beta","tau"),
      n.chains = 2, n.iter = n.iter,
      n.burnin = n.burnin,DIC=F,
      model.file = model_AD_NMR_BIN_RE_independent_UME)
    
    NMR_res<-as.data.frame(temp_nmr$summary)%>%
      select(TE_NMR=mean,seTE_NMR=sd,
             lci_NMR=`2.5%`,
             uci_NMR=`97.5%`,
             Rhat)%>%
      tibble::rownames_to_column()%>%
      mutate(paramtype=gsub("\\[.*","",rowname),
             tempname=gsub(".*\\[|\\],*","",rowname))%>%
      tidyr::separate(tempname, into=c("comp_num","ref_num"),sep = ",")%>%
      mutate(Comparator=LETTERS[as.numeric(comp_num)],
             Reference =LETTERS[as.numeric(ref_num)],
             paramtype=case_when(paramtype=="lor0"~"es",
                                 paramtype=="beta"~"beta_x",
                                 TRUE~paramtype))%>%
      left_join(refdat,by=c("Comparator", "Reference","paramtype"))
    
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  rm(temp_nmr)
  gc()
  return(NMR_res)
  
}



# run models in clusters ############################
if (!is_cluster_active()) {
  cl <- startjagscl()
}

# if(is_cluster_active()==TRUE) {
#   stopCluster(cl)
# }
# 

ws_snapshot<-ls()
## NMA #################
clusterEvalQ(cl, {set.seed(jagsseedval)})
clusterExport(cl, "run_NMA_JAGS")
NMA_res_JAGS=parLapply(cl,1:N.sim, function(x) run_NMA_JAGS(simid=x))
# NMA_res_JAGS=parLapply(cl,1:N.sim, function(x) run_NMA_JAGS(all_JAGS_lists[[x]]))
NMA_res_JAGS<-bind_rows(NMA_res_JAGS, .id = "simulation")%>%
  rename(es_logOR=Estimate)%>%
  mutate(jagsseedval=jagsseedval,bias=TE_NMR -es_logOR,
         coverage_TE=ifelse(es_logOR>=lci_NMR &es_logOR<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )%>%
  mutate(fileloc=unique(data2[[1]]$fileloc))
# write.csv(NMA_res_JAGS,paste(loopdir,"/NMAjagsres_binary_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".csv",sep=""))
saveRDS(NMA_res_JAGS, file = paste(loopdir,"/NMAjagsres_binary_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".rds",sep=""), compress="xz")


gc()
rm(list=ls()[!ls()%in%c(ws_snapshot,"ws_snapshot")])
if (!is_cluster_active()) {
  cl <- startjagscl()
}

## IND ###############
clusterEvalQ(cl, {set.seed(jagsseedval)})
clusterExport(cl, "run_ind_JAGS")
ind_res_JAGS=parLapply(cl,1:N.sim, function(x) run_ind_JAGS(simid=x))
# ind_res_JAGS=parLapply(cl,1:N.sim, function(x) run_ind_JAGS(all_JAGS_lists[[x]]))
ind_res_JAGS<-bind_rows(ind_res_JAGS, .id = "simulation")%>%
  rename(es_logOR=Estimate)%>%
  mutate(jagsseedval=jagsseedval,bias=TE_NMR -es_logOR,
         coverage_TE=ifelse(es_logOR>=lci_NMR &es_logOR<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )%>%
  mutate(fileloc=unique(data2[[1]]$fileloc))
# write.csv(ind_res_JAGS,paste(loopdir,"/Ijagsres_binary_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".csv",sep=""))
saveRDS(ind_res_JAGS, file = paste(loopdir,"/Ijagsres_binary_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".rds",sep=""), compress="xz")

gc()
rm(list=ls()[!ls()%in%c(ws_snapshot,"ws_snapshot")])
if (!is_cluster_active()) {
  cl <- startjagscl()
}


## IND UMIE ###############

clusterEvalQ(cl, {set.seed(jagsseedval)})
clusterExport(cl, "run_ind_UMIE_JAGS")
ind_UMIE_res_JAGS=parLapply(cl,1:N.sim, function(x) run_ind_UMIE_JAGS(simid=x))
# ind_UMIE_res_JAGS=parLapply(cl,1:N.sim, function(x) run_ind_UMIE_JAGS(all_JAGS_lists[[x]]))
ind_UMIE_res_JAGS<-bind_rows(ind_UMIE_res_JAGS, .id = "simulation")%>%
  rename(es_logOR=Estimate)%>%
  mutate(jagsseedval=jagsseedval,bias=TE_NMR -es_logOR,
         coverage_TE=ifelse(es_logOR>=lci_NMR &es_logOR<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )%>%
  mutate(fileloc=unique(data2[[1]]$fileloc))
# write.csv(ind_UMIE_res_JAGS,paste(loopdir,"/IUjagsres_binary_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".csv",sep=""))
saveRDS(ind_UMIE_res_JAGS, file = paste(loopdir,"/IUjagsres_binary_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".rds",sep=""), compress="xz")


gc()
rm(list=ls()[!ls()%in%c(ws_snapshot,"ws_snapshot")])
if (!is_cluster_active()) {
  cl <- startjagscl()
}


## COMMON ###############

clusterEvalQ(cl, {set.seed(jagsseedval)})
clusterExport(cl, "run_common_JAGS")
common_res_JAGS=parLapply(cl,1:N.sim, function(x) run_common_JAGS(simid=x))
# common_res_JAGS=parLapply(cl,1:N.sim, function(x) run_common_JAGS(all_JAGS_lists[[x]]))
common_res_JAGS<-bind_rows(common_res_JAGS, .id = "simulation")%>%
  rename(es_logOR=Estimate)%>%
  mutate(jagsseedval=jagsseedval,bias=TE_NMR -es_logOR,
         coverage_TE=ifelse(es_logOR>=lci_NMR &es_logOR<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )%>%
  mutate(fileloc=unique(data2[[1]]$fileloc))
# write.csv(common_res_JAGS,paste(loopdir,"/Cjagsres_binary_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".csv",sep=""))
saveRDS(common_res_JAGS, file = paste(loopdir,"/Cjagsres_binary_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".rds",sep=""), compress="xz")


gc()
rm(list=ls()[!ls()%in%c(ws_snapshot,"ws_snapshot")])
if (!is_cluster_active()) {
  cl <- startjagscl()
}

## COMMON UMIE ###############
clusterEvalQ(cl, {set.seed(jagsseedval)})
clusterExport(cl, "run_common_UMIE_JAGS")
common_UMIE_res_JAGS=parLapply(cl,1:N.sim, function(x) run_common_UMIE_JAGS(simid=x))
# common_UMIE_res_JAGS=parLapply(cl,1:N.sim, function(x) run_common_UMIE_JAGS(all_JAGS_lists[[x]]))
common_UMIE_res_JAGS<-bind_rows(common_UMIE_res_JAGS, .id = "simulation")%>%
  rename(es_logOR=Estimate)%>%
  mutate(jagsseedval=jagsseedval,
         bias=TE_NMR -es_logOR,
         coverage_TE=ifelse(es_logOR>=lci_NMR &es_logOR<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )%>%
  mutate(fileloc=unique(data2[[1]]$fileloc))
# write.csv(common_UMIE_res_JAGS,paste(loopdir,"/CUjagsres_binary_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".csv",sep=""))
saveRDS(common_UMIE_res_JAGS, file = paste(loopdir,"/CUjagsres_binary_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".rds",sep=""), compress="xz")

stopCluster(cl)
gc()
rm(list=ls()[!ls()%in%c(ws_snapshot,"ws_snapshot")])
