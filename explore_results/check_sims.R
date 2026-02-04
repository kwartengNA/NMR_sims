# libs and settings (if file  not called by an external function) ######
# library(dtplyr)
# library(magrittr)
# library(dplyr)
# library(parallel)
# library(R.utils)
# library(mcstatsim)
# setwd("C:/Users/kwarteng/Documents/NMR_simulation_project/mobaxterm")
# gendir<-getwd()
# seedval <- 120 
# import results from implemented models and export (in practice this would benefit from subgrouping the runs) ##############
outputfiles<-
  list.files(paste(gendir, paste("seed_", seedval,sep=""), sep="/"),#paste(gendir, paste("seed_", seedval, "_multi",sep=""), sep="/"),
             full.names = TRUE, recursive = TRUE)%>%
  lazy_dt()%>%
  as.data.frame()%>%
  rename(val=1)%>%
  mutate(relative_loc=gsub(gendir,"",val))%>%
  tidyr::separate(col=relative_loc, into = paste("vals", 1:11, sep=""), sep = "/", remove = FALSE)%>%
  rename(seed=vals2, outcome=vals3, true_arms=vals4, consistency_txt=vals5, true_assumption =vals6, tau =vals7, density=vals8, interaction=vals9, events =vals10, output=vals11)%>%
  mutate(output=ifelse(is.na(output) & outcome=="continuous" & grepl("_events", interaction), events, output),
         events=ifelse(outcome=="continuous" & grepl("_events", interaction), "NA", events),
         true_consistency=ifelse(consistency_txt=="consistency", TRUE, FALSE),
         tau_txt=gsub("tau_","",tau), 
         keep_prop_txt=gsub("connect_","",density),
         interaction_direction=gsub("_beta","", interaction),
         event_rate=gsub("_events", "", events),
         event_rate=ifelse(event_rate=="NA", NA, event_rate))%>%
  select(-starts_with("vals"))

my_rds<-outputfiles%>%filter(grepl("\\.rds", output)==TRUE)%>%filter(grepl("^res_binary", output)==TRUE)

all_res <- vector("list", nrow(my_rds))
for(csvi in 1:nrow(my_rds)){
  all_res[[csvi]]<-readRDS(my_rds$val[csvi])%>%
    mutate(fileloc=trimws(gsub("/[^/]*$", "", filename)))}
gc()

readr::write_rds(x=bind_rows(all_res), file=paste(gendir, paste("seed_", seedval,sep=""),"all_res.rds", sep="/"), compress="xz",compression = 9L)
# create summary statistics ###########

outputfiles<-
  list.files(paste(gendir, paste("seed_", seedval,sep=""), sep="/"),#paste(gendir, paste("seed_", seedval, "_multi",sep=""), sep="/"),
             full.names = TRUE, recursive = FALSE)%>%
  as.data.frame()%>%
  rename(val=1)%>%
  tidyr::separate(col=val, into = paste("vals", 1:8, sep=""), sep = "/", remove = FALSE)%>%
  rename(seed=vals7, output=vals8)%>%
  select(-starts_with("vals"))%>%
  filter(grepl("all_res", output))

diag_res_all<-NULL
param_res_all<-NULL

keepls<-unique(c(ls(), "param_res_all", "diag_res_all", "filei", "outputfiles", "keepls"))

for(filei in c(outputfiles$val)){

  all_res <- readRDS(filei)
  
  all_res<-split(all_res, f = all_res$fileloc)
  gc()
  # prepare diag_res and param_res mean (TAU is missing) #######################
  # need to manually combine if i split all_res into files
  # get N.sim
  N.sim<- max(as.numeric(all_res[[1]]$simulation))
  
  # filter into param_res and diagres
  
  
  diag_res<-lapply(all_res, function(x){
    x%>%
      filter(paramtype %in%c("AIC", "AICc", "BIC", "dev", "ll"))%>%
      select(where(~ any(!is.na(.))))
  })
  keep_diag_res <- sapply(diag_res, function(x) nrow(x) > 0)
  diag_res <- diag_res[keep_diag_res] #if the file is a jags file it is no longer relevant
  
  param_res<-lapply(all_res, function(x){
    if("es_logOR" %in% names(x)){
      x<-x%>%
        mutate(es_NMR=coalesce(es_NMR, es_logOR))}
    
    x%>%
      filter(paramtype %in%c("es", "beta_x", "tau")|error==TRUE)%>%
      filter((comp_num>ref_num)|paramtype=="tau"|error==TRUE)%>%
      mutate(comp=coalesce(comp, paste(Comparator, Reference)))
    
  })
  
  other_res<-lapply(all_res, function(x){
    x%>%
      filter(!paramtype %in%c("AIC", "AICc", "BIC", "dev", "ll", "es", "beta_x", "tau"))%>%
      filter(is.na(error))
    
  })
  
  
  rm(all_res)
  gc()
  
  # define attributes and link to fileloc to prevent multiple separate() runs
  attribute_set<-data.frame("fileloc"=names(param_res))%>%# if imported
    mutate(fileloc=gsub("C:/Users/kwarteng/Documents/NMR_simulation_project/mobaxterm/","",fileloc))%>%
    # data.frame("fileloc"=gsub("C:/Users/kwarteng/Documents/NMR_sims/","",outputfiles$val))%>%#if continuation
    tidyr::separate(col=fileloc, into = paste("vals", 1:10, sep=""), sep = "/", remove = FALSE)%>%
    mutate(filenum=as.numeric(as.factor(fileloc)))%>%
    rename(seed=vals2, outcome=vals3,true_arms=vals4, 
           # output, 
           consistency_txt=vals5, true_assumption=vals6, tau=vals7, 
           density=vals8, interaction=vals9, events=vals10)%>%
    mutate(events=ifelse(outcome=="continuous" & grepl("_events", interaction), "NA", events),
           true_consistency=ifelse(consistency_txt=="consistency", TRUE, FALSE),
           tau_txt=gsub("tau_","",tau), 
           keep_prop_txt=gsub("connect_","",density),
           interaction_direction=gsub("_beta","", interaction),
           event_rate=gsub("_events", "", events),
           event_rate=ifelse(event_rate=="NA", NA, event_rate))%>%
    select(-starts_with("vals"))
  
  # get scenario attributes from fileloc(param_res)
  starttime<-Sys.time()
  param_res<-withTimeout({lapply(param_res, function(x){
    x%>%
      left_join(attribute_set, by=c("fileloc"))%>%
      filter((true_assumption=="NMA" & true_consistency==TRUE)|
               (true_assumption!="NMA"))%>%
      mutate(density=factor(density, levels=c("connect_all", "connect_moderate", "connect_low")),true_model=trimws(paste(ifelse(true_assumption=="independent","ind",true_assumption), ifelse(true_consistency==FALSE, "UMIE",""))))
  })}, timeout = 30*length(param_res))
  Sys.time()-starttime
  
  gc()
  
  # get scenario attributes from fileloc (diag_res)
  starttime<-Sys.time()
  diag_res<-withTimeout({lapply(diag_res, function(x){
    x%>%
      left_join(attribute_set, by=c("fileloc"))%>%
      filter((true_assumption=="NMA" & true_consistency==TRUE)|
               (true_assumption!="NMA"))%>%
      mutate(density=factor(density, levels=c("connect_all", "connect_moderate", "connect_low")),true_model=trimws(paste(ifelse(true_assumption=="independent","ind",true_assumption), ifelse(true_consistency==FALSE, "UMIE",""))))
  })}, timeout = 30*length(param_res))
  Sys.time()-starttime
  
  gc()
  
  # assign implied NMA interaction estimates in JAGS (0) and associated metrics
  starttime<-Sys.time()
  param_res<-lapply(param_res, function(x){
    x%>%
      filter(model=="NMA"&src=="RMA_MV" & paramtype=="beta_x")%>%
      mutate(rowname="manualNMAbeta",
             src="JAGS",
             # paramtype=="beta_x", TE_NMR=0, seTE_NMR=0, lci_NMR=0, uci_NMR=0, Rhat=NA,
             es_MD=es_NMR,
             rowname.y=NA,graphname=NA,graphname2=NA, jagsseedval=NA,
             bias=TE_NMR-es_MD,
             coverage_TE=ifelse(es_MD>=lci_NMR &es_MD<=uci_NMR,TRUE,FALSE),
             mse=(bias^2)+(seTE_NMR^2),
             # length=0, filename="manual",
             .comp=as.character(NA), cov_NMR=NA, error=NA, es_NMR=NA, ML=NA, REML=NA,# for RMA.MV only
             bias2=NA# I don't use it for anything
      )%>%
      bind_rows(x)
  })
  Sys.time()-starttime
  
  
  # Rename ind UMIE to ind UMIE min and id converged estimates 
  starttime<-Sys.time()
  param_res<-lapply(param_res, function(x){
    if(sum(grepl("Rhat", names(x))==TRUE)==0){
      x$Rhat<-NA
    }
    
    x%>%
      mutate(model=ifelse(model=="ind\nUMIE", "ind\nUMIE\nmin", model))%>%
      mutate(Rhat_temp=ifelse(is.na(Rhat)==TRUE, -999, Rhat), # Missing Rhats do not mean non-convergence, they are okay
             error_temp=ifelse(is.na(error)==TRUE, -999, as.numeric(error)), # the NAs just mean they didn't fail, they are okay
             seTE_NMR_temp=case_when(is.na(seTE_NMR)==TRUE & paramtype=="tau"~-999, # tau does not have CI, so the NAs are okay
                                     TRUE~seTE_NMR),
             nonest=ifelse(Rhat_temp>=1.1|(seTE_NMR_temp>300)|
                             is.na(TE_NMR)| 
                             error_temp==TRUE| # RMA.MV runs that completely failed
                             (seTE_NMR_temp==0&model!="common"), # implied 0s are meaningful in common, so they are estimations
                           TRUE, FALSE))%>%
      select(-c(Rhat_temp, error_temp, seTE_NMR_temp))
  })
  Sys.time()-starttime
  # param_res_bkup<-param_res
  
  # Proportion Reduction in values################
  
  param_res<-lapply(param_res, function(x){
    x%>%
      group_by(simulation, src, paramtype, comp)%>%
      mutate(NMA_est=ifelse(length(TE_NMR[model=="NMA"])==0,NA,TE_NMR[model=="NMA"]),
             prop_red=(TE_NMR-es_NMR)/es_NMR,
             NMA_prop_red=(TE_NMR-NMA_est)/NMA_est)%>%
      ungroup()
  })
  
  
  ## Start Optional Filter for converged estimates #################
  fail_res<-lapply(param_res, function(x){
    x%>%filter(nonest==TRUE)
  })
  param_res<-lapply(param_res, function(x){
    x%>%filter(nonest==FALSE)
  })
  gc()
  ## End Optional Filter for converged estimates #################
  
  
  # summary statistics for param res #############
  # bkupparam_res<-param_res
  param_res<-bind_rows(param_res)
  param_res<-split(param_res, list(param_res$comp, param_res$src, param_res$paramtype, param_res$model, param_res$filenum))
  param_res<-param_res[lapply(param_res, nrow)>0]
  
  
  starttime<-Sys.time()
  bias<-lapply(param_res, function(x) calc_bias(x$TE_NMR, unique(x$es_NMR)))
  
  # only calculate coverage where profile CI is available
  coverage_param_res<-lapply(param_res, function(x) x%>%filter(!(is.na(lci_NMR) & paramtype=="tau")))
  coverage<-lapply(coverage_param_res, function(x) calc_coverage(upper_bound = x$uci_NMR, lower_bound = x$lci_NMR,true_param = unique(x$es_NMR)))
  rm(coverage_param_res)
  
  mse<-lapply(param_res,function(x)calc_mse(estimates = x$TE_NMR,true_param = unique(x$es_NMR)))
  rmse<-lapply(param_res,function(x)calc_rmse(estimates = x$TE_NMR,true_param = unique(x$es_NMR)))
  variance_info<-lapply(param_res,function(x)calc_variance(estimates = x$TE_NMR))
  
  man_mse<-lapply(param_res,function(x) mean(x$mse)) # this mse is bias^2 +se^2
  names(man_mse)<-paste(names(man_mse),".man_mse",sep="")
  
  CI_length<-lapply(param_res,function(x) mean(x$length)) # this mse is bias^2 +se^2
  names(CI_length)<-paste(names(CI_length),".CI_length",sep="")
  # param_resbkup<-param_res
  param_res<-data.frame("est"=c(unlist(bias),
                                unlist(coverage),
                                unlist(mse),
                                unlist(rmse),
                                unlist(variance_info),
                                unlist(man_mse),
                                unlist(CI_length)))%>%
    mutate(rowname=rownames(.))%>%
    tidyr::separate(col = rowname, into = c("comp",  "src",
                                            "paramtype", "model", "filenum",
                                            "name"), sep="\\.")%>%
    arrange(comp)#
  
  param_res<-param_res%>%
    tidyr::pivot_wider(values_from = est, names_from = name)%>%
    mutate(filenum=as.numeric(filenum))%>%
    left_join(attribute_set%>%distinct(), by=c("filenum"))
  
  param_res<-param_res%>%
    group_by(fileloc,src,model, paramtype)%>%
    mutate(#mean_length=mean(length,na.rm=TRUE), # summarise mean values across runs and comparisons
      mean_bias=mean(bias,na.rm=TRUE),
      mean_man_mse=mean(man_mse,na.rm=TRUE),
      mean_mse=mean(mse,na.rm=TRUE),
      mean_rmse=mean(rmse,na.rm=TRUE),
      mean_variance=mean(variance,na.rm=TRUE),
      mean_CI_length=mean(CI_length,na.rm=TRUE),
      
      mean_bias_mcse=mean(bias_mcse,na.rm=TRUE),
      mean_mse_mcse=mean(mse_mcse,na.rm=TRUE),
      mean_rmse_mcse=mean(rmse_mcse,na.rm=TRUE),
      mean_variance_mcse=mean(variance_mcse,na.rm=TRUE),
      
      prop_coverage=mean(coverage,na.rm=TRUE),
      
      # median_length=median(length,na.rm=TRUE), # summarise median values across runs and comparisons
      median_bias=median(bias,na.rm=TRUE),
      median_man_mse=median(man_mse,na.rm=TRUE),
      median_mse=median(mse,na.rm=TRUE),
      median_rmse=median(rmse,na.rm=TRUE),
      median_variance=median(variance,na.rm=TRUE),
      median_coverage=median(coverage,na.rm=TRUE),
      median_CI_length=median(CI_length,na.rm=TRUE),
      
      median_bias_mcse=median(bias_mcse,na.rm=TRUE),
      median_mse_mcse=median(mse_mcse,na.rm=TRUE),
      median_rmse_mcse=median(rmse_mcse,na.rm=TRUE),
      median_variance_mcse=median(variance_mcse,na.rm=TRUE),
      
      n_groups=n())%>%
    ungroup()
  Sys.time()-starttime
  
  # combine/export param_res and diag_res#################
  diag_res_all<-c(diag_res, diag_res_all)
  param_res_all<-c(list(param_res), param_res_all)
  
  rm(list=ls()[!ls()%in%keepls])
  gc()
  print(filei)
}

diag_res<-bind_rows(diag_res_all)%>%
  mutate(scenarionum=as.numeric(as.factor(fileloc)),
         tau_txt=coalesce(tau_txt, gsub("tau_", "",tau)))

param_res<-bind_rows(param_res_all)%>%
  mutate(scenarionum=as.numeric(as.factor(fileloc)))
# export ###########################
readr::write_rds(x=diag_res, file=paste(gendir, paste("seed_", seedval,sep=""),"diag_res.rds", sep="/"), compress="xz",compression = 9L)
readr::write_rds(x=param_res, file=paste(gendir, paste("seed_", seedval,sep=""),"filter_param_res.rds", sep="/"), compress="xz",compression = 9L)
