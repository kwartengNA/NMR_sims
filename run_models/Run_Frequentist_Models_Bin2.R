# Define functions for models/output using netmetareg ##########

## NMA with try catch ##########################
run_netmeta<-function(x){ # each individual dataset

  
  tryCatch({
    mynma<-netmeta(x, method.tau = "REML", reference.group = "A")
      }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("mynma")==FALSE){
    mynma<-data.frame(error=TRUE)}
  
  )
  
  return(mynma)
  
}



## NMA results ##########################
run_NMA<-function(data_filter){ # each individual dataset
  
  mynet<-data_filter
  
  tryCatch({
    # mynet<-netmeta(tmp_pairwise, method.tau = "REML", reference.group = "A")
    temp_nmr<-netmetareg(mynet, covar=NULL, consistency = TRUE, assumption = "independent")
    temp_summary<-summary.netmetareg_old(object = temp_nmr)
    
    NMR_res<-
      temp_summary$all_combinations%>%
      mutate(error=FALSE,
             lci_d_NMR=TE_d_NMR-(qnorm(1-0.05/2)*seTE_d_NMR),
             uci_d_NMR=TE_d_NMR+(qnorm(1-0.05/2)*seTE_d_NMR),
             lci_beta_NMR=TE_beta_NMR-(qnorm(1-0.05/2)*seTE_beta_NMR),
             uci_beta_NMR=TE_beta_NMR+(qnorm(1-0.05/2)*seTE_beta_NMR))%>%
      left_join(refdat%>%filter(paramtype=="es")%>%select(es_d_NMR=Estimate,Comparator,Reference),
                by=c("treat1"="Comparator", "treat2"="Reference"))%>%
      left_join(refdat%>%filter(paramtype=="beta_x")%>%select(es_beta_NMR=Estimate,Comparator,Reference),
                by=c("treat1"="Comparator", "treat2"="Reference"))%>%
      bind_rows(  as.data.frame(as.data.frame(confint(temp_nmr)))%>%
                    tibble::rownames_to_column(var="paramtype")%>%
                    rename(TE_tau_NMR=estimate,
                           lci_tau_NMR=ci.lb,
                           uci_tau_NMR=ci.ub)%>%
                    left_join(refdat%>%filter(paramtype%in%c("tau"))%>%
                                select(paramtype,es_tau_NMR=Estimate)))%>%
      bind_rows(temp_nmr$fit.stats%>%as.data.frame()%>%tibble::rownames_to_column(var="paramtype"))
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  
  return(NMR_res)
  
}


## ind ##########################
run_ind<-function(data_filter){ # each individual dataset
  
  mynet<-data_filter
  
  tryCatch({
    #mynet<-netmeta(tmp_pairwise, method.tau = "REML", reference.group = "A")
    temp_nmr<-netmetareg(mynet, covar=xcov, consistency = TRUE, assumption = "independent")
    temp_summary<-summary.netmetareg_old(object = temp_nmr)
    
    NMR_res<-
      temp_summary$all_combinations%>%
      mutate(error=FALSE,
             lci_d_NMR=TE_d_NMR-(qnorm(1-0.05/2)*seTE_d_NMR),
             uci_d_NMR=TE_d_NMR+(qnorm(1-0.05/2)*seTE_d_NMR),
             lci_beta_NMR=TE_beta_NMR-(qnorm(1-0.05/2)*seTE_beta_NMR),
             uci_beta_NMR=TE_beta_NMR+(qnorm(1-0.05/2)*seTE_beta_NMR))%>%
      left_join(refdat%>%filter(paramtype=="es")%>%select(es_d_NMR=Estimate,Comparator,Reference),
                by=c("treat1"="Comparator", "treat2"="Reference"))%>%
      left_join(refdat%>%filter(paramtype=="beta_x")%>%select(es_beta_NMR=Estimate,Comparator,Reference),
                by=c("treat1"="Comparator", "treat2"="Reference"))%>%
      bind_rows(  as.data.frame(as.data.frame(confint(temp_nmr)))%>%
                    tibble::rownames_to_column(var="paramtype")%>%
                    rename(TE_tau_NMR=estimate,
                           lci_tau_NMR=ci.lb,
                           uci_tau_NMR=ci.ub)%>%
                    left_join(refdat%>%filter(paramtype%in%c("tau"))%>%
                                select(paramtype,es_tau_NMR=Estimate)))%>%
      bind_rows(temp_nmr$fit.stats%>%as.data.frame()%>%tibble::rownames_to_column(var="paramtype"))
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  
  return(NMR_res)
  
}

## ind umie ##########################
run_ind_UMIE<-function(data_filter){ # each individual dataset
  
  mynet<-data_filter
  
  tryCatch({
    #mynet<-netmeta(tmp_pairwise, method.tau = "REML", reference.group = "A")
    temp_nmr<-netmetareg(mynet, covar=xcov, consistency = FALSE, assumption = "independent")
    temp_summary<-summary.netmetareg_old(object = temp_nmr)
    
    NMR_res<-
      temp_summary$all_combinations%>%
      mutate(error=FALSE,
             lci_d_NMR=TE_d_NMR-(qnorm(1-0.05/2)*seTE_d_NMR),
             uci_d_NMR=TE_d_NMR+(qnorm(1-0.05/2)*seTE_d_NMR),
             lci_beta_NMR=TE_beta_NMR-(qnorm(1-0.05/2)*seTE_beta_NMR),
             uci_beta_NMR=TE_beta_NMR+(qnorm(1-0.05/2)*seTE_beta_NMR))%>%
      left_join(refdat%>%filter(paramtype=="es")%>%select(es_d_NMR=Estimate,Comparator,Reference),
                by=c("treat1"="Comparator", "treat2"="Reference"))%>%
      left_join(refdat%>%filter(paramtype=="beta_x")%>%select(es_beta_NMR=Estimate,Comparator,Reference),
                by=c("treat1"="Comparator", "treat2"="Reference"))%>%
      bind_rows(  as.data.frame(as.data.frame(confint(temp_nmr)))%>%
                    tibble::rownames_to_column(var="paramtype")%>%
                    rename(TE_tau_NMR=estimate,
                           lci_tau_NMR=ci.lb,
                           uci_tau_NMR=ci.ub)%>%
                    left_join(refdat%>%filter(paramtype%in%c("tau"))%>%
                                select(paramtype,es_tau_NMR=Estimate)))%>%
      bind_rows(temp_nmr$fit.stats%>%as.data.frame()%>%tibble::rownames_to_column(var="paramtype"))
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  
  return(NMR_res)
  
}
## common ##########################
run_common<-function(data_filter){ # each individual dataset
  
  mynet<-data_filter
  
  tryCatch({
    #mynet<-netmeta(tmp_pairwise, method.tau = "REML", reference.group = "A")
    temp_nmr<-netmetareg(mynet, covar=xcov, consistency = TRUE, assumption = "common")
    temp_summary<-summary.netmetareg_old(object = temp_nmr)
    
    NMR_res<-
      temp_summary$all_combinations%>%
      mutate(error=FALSE,
             lci_d_NMR=TE_d_NMR-(qnorm(1-0.05/2)*seTE_d_NMR),
             uci_d_NMR=TE_d_NMR+(qnorm(1-0.05/2)*seTE_d_NMR),
             lci_beta_NMR=TE_beta_NMR-(qnorm(1-0.05/2)*seTE_beta_NMR),
             uci_beta_NMR=TE_beta_NMR+(qnorm(1-0.05/2)*seTE_beta_NMR))%>%
      left_join(refdat%>%filter(paramtype=="es")%>%select(es_d_NMR=Estimate,Comparator,Reference),
                by=c("treat1"="Comparator", "treat2"="Reference"))%>%
      left_join(refdat%>%filter(paramtype=="beta_x")%>%select(es_beta_NMR=Estimate,Comparator,Reference),
                by=c("treat1"="Comparator", "treat2"="Reference"))%>%
      bind_rows(  as.data.frame(as.data.frame(confint(temp_nmr)))%>%
                    tibble::rownames_to_column(var="paramtype")%>%
                    rename(TE_tau_NMR=estimate,
                           lci_tau_NMR=ci.lb,
                           uci_tau_NMR=ci.ub)%>%
                    left_join(refdat%>%filter(paramtype%in%c("tau"))%>%
                                select(paramtype,es_tau_NMR=Estimate)))%>%
      bind_rows(temp_nmr$fit.stats%>%as.data.frame()%>%tibble::rownames_to_column(var="paramtype"))
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  
  return(NMR_res)
  
}

## ind umie min ##########################
run_ind_UMIEmax<-function(data_filter){ # each individual dataset
  
  mynet<-data_filter
  
  tryCatch({
    #mynet<-netmeta(tmp_pairwise, method.tau = "REML", reference.group = "A")
    temp_nmr<-netmetareg(mynet, covar=xcov, consistency = FALSE, assumption = "independent", max_interactions=TRUE)
    temp_summary<-summary.netmetareg_old(object = temp_nmr)
    
    NMR_res<-
      temp_summary$all_combinations%>%
      mutate(error=FALSE,
             lci_d_NMR=TE_d_NMR-(qnorm(1-0.05/2)*seTE_d_NMR),
             uci_d_NMR=TE_d_NMR+(qnorm(1-0.05/2)*seTE_d_NMR),
             lci_beta_NMR=TE_beta_NMR-(qnorm(1-0.05/2)*seTE_beta_NMR),
             uci_beta_NMR=TE_beta_NMR+(qnorm(1-0.05/2)*seTE_beta_NMR))%>%
      left_join(refdat%>%filter(paramtype=="es")%>%select(es_d_NMR=Estimate,Comparator,Reference),
                by=c("treat1"="Comparator", "treat2"="Reference"))%>%
      left_join(refdat%>%filter(paramtype=="beta_x")%>%select(es_beta_NMR=Estimate,Comparator,Reference),
                by=c("treat1"="Comparator", "treat2"="Reference"))%>%
      bind_rows(  as.data.frame(as.data.frame(confint(temp_nmr)))%>%
                    tibble::rownames_to_column(var="paramtype")%>%
                    rename(TE_tau_NMR=estimate,
                           lci_tau_NMR=ci.lb,
                           uci_tau_NMR=ci.ub)%>%
                    left_join(refdat%>%filter(paramtype%in%c("tau"))%>%
                                select(paramtype,es_tau_NMR=Estimate)))%>%
      bind_rows(temp_nmr$fit.stats%>%as.data.frame()%>%tibble::rownames_to_column(var="paramtype"))
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  
  return(NMR_res)
  
}


## common umie ##########################
run_common_UMIE<-function(data_filter){ # each individual dataset
  
  mynet<-data_filter
  
  tryCatch({
    #mynet<-netmeta(tmp_pairwise, method.tau = "REML", reference.group = "A")
    temp_nmr<-netmetareg(mynet, covar=xcov, consistency = FALSE, assumption = "common")
    temp_summary<-summary.netmetareg_old(object = temp_nmr)
    
    NMR_res<-
      temp_summary$all_combinations%>%
      mutate(error=FALSE,
             lci_d_NMR=TE_d_NMR-(qnorm(1-0.05/2)*seTE_d_NMR),
             uci_d_NMR=TE_d_NMR+(qnorm(1-0.05/2)*seTE_d_NMR),
             lci_beta_NMR=TE_beta_NMR-(qnorm(1-0.05/2)*seTE_beta_NMR),
             uci_beta_NMR=TE_beta_NMR+(qnorm(1-0.05/2)*seTE_beta_NMR))%>%
      left_join(refdat%>%filter(paramtype=="es")%>%select(es_d_NMR=Estimate,Comparator,Reference),
                by=c("treat1"="Comparator", "treat2"="Reference"))%>%
      left_join(refdat%>%filter(paramtype=="beta_x")%>%select(es_beta_NMR=Estimate,Comparator,Reference),
                by=c("treat1"="Comparator", "treat2"="Reference"))%>%
      bind_rows(  as.data.frame(as.data.frame(confint(temp_nmr)))%>%
                    tibble::rownames_to_column(var="paramtype")%>%
                    rename(TE_tau_NMR=estimate,
                           lci_tau_NMR=ci.lb,
                           uci_tau_NMR=ci.ub)%>%
                    left_join(refdat%>%filter(paramtype%in%c("tau"))%>%
                                select(paramtype,es_tau_NMR=Estimate)))%>%
      bind_rows(temp_nmr$fit.stats%>%as.data.frame()%>%tibble::rownames_to_column(var="paramtype"))
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  
  return(NMR_res)
  
}

# run models in clusters #########################


# clusterExport(cl,"X1")
# clusterExport(cl,"NT")
# num_cores <- floor(detectCores()/2)
cl <- makeCluster(num_cores)
clusterExport(cl,c("N.sim",
                   "es_lOR","beta_x","refdat"))
clusterEvalQ(cl, {.libPaths( c( "/h/kwarteng/mypkgs" , .libPaths() ) )}) # needed for cluster
clusterEvalQ(cl, {library(netmeta)})
clusterEvalQ(cl, {library(dplyr)})
clusterEvalQ(cl, {source("netmetareg.R")})

# while netmetareg is internal
clusterEvalQ(cl, {
  chknumeric <- netmeta:::chknumeric
  chklevel <- netmeta:::chklevel
  chkclass <- netmeta:::chkclass
  updateversion <- netmeta:::updateversion
  chklogical <- netmeta:::chklogical
  setchar <- netmeta:::setchar
  setref <- netmeta:::setref
  catch <- netmeta:::catch
  calcV <- netmeta:::calcV
  runNN <- netmeta:::runNN
  formatN<-netmeta:::formatN
  contrmat <- metafor::contrmat
  bldiag <- metafor::bldiag
  rma.mv <- metafor::rma.mv
  formatCI<-netmeta:::formatCI
  replaceNA<-netmeta:::replaceNA
  formatPT<-netmeta:::formatPT
  # as_df<-netmeta:::as_df
})

clusterExport(cl, "data2")
data2_pairwise<-parLapply(cl,data2, function(x) pairwise(treat=treatment,n = n,event = events,studlab=study,data=x,sm="OR"))

clusterExport(cl,c("data2_pairwise"))
clusterExport(cl, "run_netmeta") #netmeta with trycatch
data2_netmeta<-parLapply(cl,1:N.sim, function(x)run_netmeta(data2_pairwise[[x]]))


clusterExport(cl,c("data2_netmeta"))
clusterExport(cl, "run_NMA")
NMA_res=parLapply(cl,1:N.sim, function(x) run_NMA(data2_netmeta[[x]]))

clusterExport(cl, "run_ind")
ind_res=parLapply(cl,1:N.sim, function(x) run_ind(data2_netmeta[[x]]))

clusterExport(cl, "run_ind_UMIE")
ind_UMIE_res=parLapply(cl,1:N.sim, function(x) run_ind_UMIE(data2_netmeta[[x]]))

clusterExport(cl, "run_common")
common_res=parLapply(cl,1:N.sim, function(x) run_common(data2_netmeta[[x]]))

clusterExport(cl, "run_common_UMIE")
common_UMIE_res=parLapply(cl,1:N.sim, function(x) run_common_UMIE(data2_netmeta[[x]]))

clusterExport(cl, "run_ind_UMIEmax")
ind_UMIE_resmin=parLapply(cl,1:N.sim, function(x) run_ind_UMIEmax(data2_netmeta[[x]]))

# extract bias estimates ################
## NMA #################
NMA_res<-bind_rows(NMA_res, .id = "simulation")%>%
  mutate(across(.cols=all_of(c("TE_beta_NMR", "seTE_beta_NMR", "lci_beta_NMR", "uci_beta_NMR")), ~ ifelse(is.na(.comp)==FALSE, 0, .x)))%>%
  tidyr::pivot_longer(cols = c(contains("_beta_"),
                               contains("_d_"),
                               contains("_tau_"),
                               "ML",
                               "REML"))%>%
  filter(is.na(value)==FALSE|error==TRUE)%>%
  mutate(paramtype=case_when(grepl("_beta",name)~"beta",
                             grepl("_d",name)~"d",
                             TRUE~paramtype),
         name_new=gsub("_beta|_d|_tau","",name))%>%
  select(-name)%>%
  tidyr::pivot_wider(names_from = name_new,values_from = value)%>%
  mutate(bias=TE_NMR -es_NMR,
         coverage_TE=ifelse(es_NMR>=lci_NMR &es_NMR<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )
## IND ###############
ind_res<-bind_rows(ind_res, .id = "simulation")%>%
  tidyr::pivot_longer(cols = c(contains("_beta_"),
                               contains("_d_"),
                               contains("_tau_"),
                               "ML",
                               "REML"))%>%
  filter(is.na(value)==FALSE|error==TRUE)%>%
  mutate(paramtype=case_when(grepl("_beta",name)~"beta",
                             grepl("_d",name)~"d",
                             TRUE~paramtype),
         name_new=gsub("_beta|_d|_tau","",name))%>%
  select(-name)%>%
  tidyr::pivot_wider(names_from = name_new,values_from = value)%>%
  mutate(bias=TE_NMR -es_NMR,
         coverage_TE=ifelse(es_NMR>=lci_NMR &es_NMR<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )
## COMMON ###############
common_res<-bind_rows(common_res, .id = "simulation")%>%
  tidyr::pivot_longer(cols = c(contains("_beta_"),
                               contains("_d_"),
                               contains("_tau_"),
                               "ML",
                               "REML"))%>%
  filter(is.na(value)==FALSE|error==TRUE)%>%
  mutate(paramtype=case_when(grepl("_beta",name)~"beta",
                             grepl("_d",name)~"d",
                             TRUE~paramtype),
         name_new=gsub("_beta|_d|_tau","",name))%>%
  select(-name)%>%
  tidyr::pivot_wider(names_from = name_new,values_from = value)%>%
  mutate(bias=TE_NMR -es_NMR,
         coverage_TE=ifelse(es_NMR>=lci_NMR &es_NMR<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )

## COMMON UMIE ###############
common_UMIE_res<-bind_rows(common_UMIE_res, .id = "simulation")%>%
  tidyr::pivot_longer(cols = c(contains("_beta_"),
                               contains("_d_"),
                               contains("_tau_"),
                               "ML",
                               "REML"))%>%
  filter(is.na(value)==FALSE|error==TRUE)%>%
  mutate(paramtype=case_when(grepl("_beta",name)~"beta",
                             grepl("_d",name)~"d",
                             TRUE~paramtype),
         name_new=gsub("_beta|_d|_tau","",name))%>%
  select(-name)%>%
  tidyr::pivot_wider(names_from = name_new,values_from = value)%>%
  mutate(bias=TE_NMR -es_NMR,
         coverage_TE=ifelse(es_NMR>=lci_NMR &es_NMR<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )
## IND UMIE ###############
ind_UMIE_res<-bind_rows(ind_UMIE_res, .id = "simulation")%>%
  tidyr::pivot_longer(cols = c(contains("_beta_"),
                               contains("_d_"),
                               contains("_tau_"),
                               "ML",
                               "REML"))%>%
  filter(is.na(value)==FALSE|error==TRUE)%>%
  mutate(paramtype=case_when(grepl("_beta",name)~"beta",
                             grepl("_d",name)~"d",
                             TRUE~paramtype),
         name_new=gsub("_beta|_d|_tau","",name))%>%
  select(-name)%>%
  tidyr::pivot_wider(names_from = name_new,values_from = value)%>%
  mutate(bias=TE_NMR -es_NMR,
         coverage_TE=ifelse(es_NMR>=lci_NMR &es_NMR<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )


## IND UMIE v2 ###############
ind_UMIE_resmin<-bind_rows(ind_UMIE_resmin, .id = "simulation")%>%
  tidyr::pivot_longer(cols = c(contains("_beta_"),
                               contains("_d_"),
                               contains("_tau_"),
                               "ML",
                               "REML"))%>%
  filter(is.na(value)==FALSE|error==TRUE)%>%
  mutate(paramtype=case_when(grepl("_beta",name)~"beta",
                             grepl("_d",name)~"d",
                             TRUE~paramtype),
         name_new=gsub("_beta|_d|_tau","",name))%>%
  select(-name)%>%
  tidyr::pivot_wider(names_from = name_new,values_from = value)%>%
  mutate(bias=TE_NMR -es_NMR,
         coverage_TE=ifelse(es_NMR>=lci_NMR &es_NMR<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )

# combine res ######################
all_freq_res_df<-bind_rows("ind\nUMIE\nmax"=ind_UMIE_res,
							"ind\nUMIE\nmin"=ind_UMIE_resmin,
                           "ind"=ind_res,
                           "common"=common_res,
                           "common\nUMIE"=common_UMIE_res,
                           "NMA"=NMA_res,
                           .id = "model")%>%
  mutate(src="RMA_MV",
         comp=gsub("_vs_"," ",.comp),
         comp_num=as.numeric(factor(treat1,levels=LETTERS)),
         ref_num=as.numeric(factor(treat2,levels=LETTERS)),
         paramtype=case_when(paramtype=="d"~"es",
                             paramtype=="beta"~"beta_x",
                             TRUE~paramtype))%>%
  rename(Reference=treat2,
         Comparator=treat1)
