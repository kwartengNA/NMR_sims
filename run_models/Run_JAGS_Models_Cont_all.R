#WARNING HARDCODING
# Bayesian Continuous Outcome
#hardwd<-"C:/Users/kwarteng/Documents/NMR_git"
hardwd<-loopdir

# function to check clusters ########
# Check if all nodes are responsive
is_cluster_active <- function(){tryCatch({
  response <- parLapply(cl, 1:length(cl), function(x) TRUE)
  all(unlist(response))  # Check if all nodes responded
}, error = function(e) FALSE)  # If an error occurs, assume inactive)
  }

# function to start clusters ################
startjagscl<-function(){
  if(exists("cl")==FALSE|is_cluster_active()==FALSE){
  cl <- makeCluster(num_cores, outfile="temp.txt")
  clusterExport(cl,c("N.sim","data2","es_MD","beta_x","refdat","all_JAGS_lists", "model_AD_NMA_CONT_RE",
                     "model_AD_NMR_CONT_RE_const","model_AD_NMR_CONT_RE_const_UME",
                     "model_AD_NMR_CONT_RE_independent","model_AD_NMR_CONT_RE_independent_UME",
                     "n.iter","n.burnin","seedval", "jagsseedval"))
  clusterEvalQ(cl, {.libPaths( c( "/h/kwarteng/mypkgs" , .libPaths() ) )}) 
  clusterEvalQ(cl, {library(dplyr)})
  clusterEvalQ(cl, {library(jagsUI)})#, lib.loc = "/h/kwarteng/mypkgs"
  clusterEvalQ(cl, {set.seed(jagsseedval)})# if the JAGS models are run separately, they require the same seed in order to obtain the same results
  }
  return(cl)
  }
# Model files iterations, and burnin ####################
model_AD_NMR_CONT_RE_independent<-"JAGS/model_AD_NMR_CONT_RE_independent.txt"
model_AD_NMR_CONT_RE_const<-"JAGS/model_AD_NMR_CONT_RE_constant.txt"
model_AD_NMA_CONT_RE<-"JAGS/model_AD_NMA_CONT_RE.txt"
model_AD_NMR_CONT_RE_const_UME<-"JAGS/model_AD_NMR_CONT_RE_constant_inconsistency_Idir.txt"
model_AD_NMR_CONT_RE_independent_UME<-"JAGS/model_AD_NMR_CONT_RE_independent_inconsistency_Idir.txt"
n.iter = 30000
n.burnin = 10000

# JAGS data prep ###################


all_JAGS_lists<-NULL
for(i in 1:N.sim){
  temp_mydata_wide<-
    data2[[i]]%>%
    mutate(
      treat_num=treat,
      # Ivalrand=ifelse(Ival!=0 & Ivalrand==0,-1,Ivalrand), # for JAGS it needs to be 1 and -1, unless neither ref nor comp is favored.
      seTE=sd/sqrt(n),
      prec=seTE^-2)%>%
    group_by(study)%>%
    mutate(nstudy=sum(n),
           treat_order=row_number())%>%
    ungroup()%>%
    tidyr::pivot_wider(names_from = treat_order,
                       values_from =c("yi","n","treat","treat_num","treatment","sd", "prec","seTE", ,starts_with("Ival_")))#,"Ivalrand"#"Ival_std" "Ival_gen"
  
  all_JAGS_lists[[i]]<-list("ns"=length(unique(temp_mydata_wide$study)),
                            "nt"=max(temp_mydata_wide%>%
                                       select(starts_with("treat_num")),na.rm=TRUE),
                            "na"=rowSums(is.na(select(temp_mydata_wide,
                                                      starts_with("treat_num")))==FALSE),
                            "x"=temp_mydata_wide$xcov,
                            "mx"=0,#mean(temp_mydata_wide$xcov),
                            ## need a column for each arm
                            "t"=temp_mydata_wide%>%select(starts_with("treat_num")),
                            "y"=temp_mydata_wide%>%select(starts_with("yi")),
                            ## need a column for each arm
                            "n"=temp_mydata_wide%>%select(starts_with("n_")),
                            "prec"= temp_mydata_wide%>%select(starts_with("prec_")),
                            "I_dir"= temp_mydata_wide%>%select(starts_with("Ival_std_"))
  )
}


# Define functions for models/output ####################



## NMA ################
run_NMA_JAGS<-function(data_filter){ # each individual dataset
  tryCatch({
  temp_nma<- jags(
    data = data_filter,
    inits = NULL,
    parameters.to.save = c("ddiff","sd"),
    n.chains = 2, n.iter = n.iter,
    n.burnin = n.burnin,DIC=F,
    model.file = model_AD_NMA_CONT_RE)
  
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
           paramtype=case_when(paramtype=="d"~"es",
                               paramtype=="ddiff"~"es",
                               paramtype=="betadiff"~"beta_x",
                               paramtype=="sd"~"tau",
                               TRUE~paramtype))%>%
    left_join(refdat,by=c("Comparator", "Reference","paramtype"))
  
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMA_res")==FALSE){
    NMA_res<-data.frame(error=TRUE)}
  
  )
  return(NMA_res)
  
}



## IND ################
run_ind_JAGS<-function(data_filter){ # each individual dataset
  tryCatch({
  temp_nmr<- jags(
    data = data_filter,
    inits = NULL,
    parameters.to.save = c("ddiff", "ddiff0", "betadiff","sd"),
    n.chains = 2, n.iter = n.iter,
    n.burnin = n.burnin,DIC=F,
    model.file = model_AD_NMR_CONT_RE_independent)
  
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
           paramtype=case_when(paramtype=="ddiff0"~"es",
                               paramtype=="betadiff"~"beta_x",
                               paramtype=="sd"~"tau",
                               TRUE~paramtype))%>%
    left_join(refdat,by=c("Comparator", "Reference","paramtype"))
  
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  return(NMR_res)
  
}

## COMMON ########################
run_common_JAGS<-function(data_filter){ # each individual dataset
  tryCatch({
  temp_nmr<- jags(
    data = data_filter,
    inits = NULL,
    parameters.to.save = c("ddiff","ddiff0", "betadiff","sd"),
    n.chains = 2, n.iter = n.iter,
    n.burnin = n.burnin,DIC=F,
    model.file = model_AD_NMR_CONT_RE_const)
  
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
           paramtype=case_when(paramtype=="ddiff0"~"es",
                               paramtype=="betadiff"~"beta_x",
                               paramtype=="sd"~"tau",
                               TRUE~paramtype))%>%
    left_join(refdat,by=c("Comparator", "Reference","paramtype"))
  
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  return(NMR_res)
  
}

## COMMON UMIE ##############

run_common_UMIE_JAGS<-function(data_filter){ # each individual dataset
  tryCatch({
  temp_nmr<- jags(
    data = data_filter,
    inits = NULL,
    parameters.to.save = c("ddiff0", "ddiff", "beta","tau"),
    n.chains = 2, n.iter = n.iter,
    n.burnin = n.burnin,DIC=F,
    model.file = model_AD_NMR_CONT_RE_const_UME)
  
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
           paramtype=case_when(paramtype=="ddiff0"~"es",
                               paramtype=="beta"~"beta_x",
                               TRUE~paramtype))%>%
    left_join(refdat,by=c("Comparator", "Reference","paramtype"))
  
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  return(NMR_res)
  
}

## IND UMIE #########################
run_ind_UMIE_JAGS<-function(data_filter){ # each individual dataset
  tryCatch({
  temp_nmr<- jags(
    data = data_filter,
    inits = NULL,
    parameters.to.save = c("ddiff0", "ddiff", "beta","tau"),
    n.chains = 2, n.iter = n.iter,
    n.burnin = n.burnin,DIC=F,
    model.file = model_AD_NMR_CONT_RE_independent_UME)
  
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
           paramtype=case_when(paramtype=="ddiff0"~"es",
                               paramtype=="beta"~"beta_x",
                               TRUE~paramtype))%>%
    left_join(refdat,by=c("Comparator", "Reference","paramtype"))
  
  }, 
  
  error=function(e){cat("ERROR in simulation: ",conditionMessage(e), "\n")},
  
  finally = if(exists("NMR_res")==FALSE){
    NMR_res<-data.frame(error=TRUE)}
  
  )
  return(NMR_res)
  
}



# run models in clusters  ###################
if(is_cluster_active()==TRUE) {
  stopCluster(cl)
}

cl<-startjagscl()
ws_snapshot<-ls()
## NMA #################
clusterEvalQ(cl, {set.seed(jagsseedval)})
clusterExport(cl, "run_NMA_JAGS")
NMA_res_JAGS=parLapply(cl,1:N.sim, function(x) run_NMA_JAGS(all_JAGS_lists[[x]]))
NMA_res_JAGS<-bind_rows(NMA_res_JAGS, .id = "simulation")%>%
   rename(es_MD=Estimate)%>%
   mutate(jagsseedval=jagsseedval,bias=TE_NMR -es_MD,
            coverage_TE=ifelse(es_MD>=lci_NMR &es_MD<=uci_NMR,TRUE,FALSE),
            mse=(bias^2)+(seTE_NMR^2),
          length=uci_NMR-lci_NMR )
write.csv(NMA_res_JAGS,paste(hardwd,"/NMAjagsres_continuous_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".csv",sep=""))
stopCluster(cl)
gc()
rm(list=ls()[!ls()%in%c(ws_snapshot,"ws_snapshot")])
cl<-startjagscl()
## IND ###############
clusterEvalQ(cl, {set.seed(jagsseedval)})
clusterExport(cl, "run_ind_JAGS")
ind_res_JAGS=parLapply(cl,1:N.sim, function(x) run_ind_JAGS(all_JAGS_lists[[x]]))
ind_res_JAGS<-bind_rows(ind_res_JAGS, .id = "simulation")%>%
  rename(es_MD=Estimate)%>%
  mutate(jagsseedval=jagsseedval, bias=TE_NMR -es_MD,
         coverage_TE=ifelse(es_MD>=lci_NMR &es_MD<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )
write.csv(ind_res_JAGS,paste(hardwd,"/Ijagsres_continuous_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".csv",sep=""))
stopCluster(cl)
gc()
rm(list=ls()[!ls()%in%c(ws_snapshot,"ws_snapshot")])
cl<-startjagscl()

## IND UMIE ###############

clusterEvalQ(cl, {set.seed(jagsseedval)})
clusterExport(cl, "run_ind_UMIE_JAGS")
ind_UMIE_res_JAGS=parLapply(cl,1:N.sim, function(x) run_ind_UMIE_JAGS(all_JAGS_lists[[x]]))

ind_UMIE_res_JAGS<-bind_rows(ind_UMIE_res_JAGS, .id = "simulation")%>%
  rename(es_MD=Estimate)%>%
  mutate(jagsseedval=jagsseedval,bias=TE_NMR -es_MD,
         coverage_TE=ifelse(es_MD>=lci_NMR &es_MD<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )
write.csv(ind_UMIE_res_JAGS,paste(hardwd,"/IUjagsres_continuous_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".csv",sep=""))
stopCluster(cl)
gc()
rm(list=ls()[!ls()%in%c(ws_snapshot,"ws_snapshot")])
cl<-startjagscl()

## COMMON ###############

clusterEvalQ(cl, {set.seed(jagsseedval)})
clusterExport(cl, "run_common_JAGS")
common_res_JAGS=parLapply(cl,1:N.sim, function(x) run_common_JAGS(all_JAGS_lists[[x]]))
common_res_JAGS<-bind_rows(common_res_JAGS, .id = "simulation")%>%
  rename(es_MD=Estimate)%>%
  mutate(jagsseedval=jagsseedval,bias=TE_NMR -es_MD,
         coverage_TE=ifelse(es_MD>=lci_NMR &es_MD<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )
write.csv(common_res_JAGS,paste(hardwd,"/Cjagsres_continuous_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".csv",sep=""))
stopCluster(cl)
gc()
rm(list=ls()[!ls()%in%c(ws_snapshot,"ws_snapshot")])
cl<-startjagscl()
## COMMON UMIE ###############
clusterEvalQ(cl, {set.seed(jagsseedval)})
clusterExport(cl, "run_common_UMIE_JAGS")
common_UMIE_res_JAGS=parLapply(cl,1:N.sim, function(x) run_common_UMIE_JAGS(all_JAGS_lists[[x]]))
common_UMIE_res_JAGS<-bind_rows(common_UMIE_res_JAGS, .id = "simulation")%>%
  rename(es_MD=Estimate)%>%
  mutate(jagsseedval=jagsseedval,
		  bias=TE_NMR -es_MD,
         coverage_TE=ifelse(es_MD>=lci_NMR &es_MD<=uci_NMR,TRUE,FALSE),
         mse=(bias^2)+(seTE_NMR^2),
         length=uci_NMR-lci_NMR )
write.csv(common_UMIE_res_JAGS,paste(hardwd,"/CUjagsres_continuous_",true_assumption,ifelse(true_consistency==FALSE,"_UMIE_","_"),true_arms,"_",seedval,".csv",sep=""))
stopCluster(cl)
gc()
rm(list=ls()[!ls()%in%c(ws_snapshot,"ws_snapshot")])


# combine res ###################################
outputdir<-paste(hardwd,"/",sep="")
myfiles<-list.files(outputdir,full.names = TRUE)
my_jagscsv<-myfiles[grepl("\\.csv", myfiles)==TRUE]
my_jagscsv<-my_jagscsv[grepl("jagsres", my_jagscsv)==TRUE]
my_jagscsv<-my_jagscsv[grepl(seedval, my_jagscsv)==TRUE]
my_jagscsv<-my_jagscsv[grepl(paste("_",true_assumption,"_",sep=""), my_jagscsv)==TRUE]
my_jagscsv<-my_jagscsv[as.numeric(grepl("UMIE", my_jagscsv))==((true_consistency-1)*-1)]
my_jagscsv<-my_jagscsv[grepl(true_arms, my_jagscsv)==TRUE]

all_jags_res_df<-NULL
for(csvi in 1:length(my_jagscsv)){
  all_jags_res_df<-read.csv(my_jagscsv[csvi])%>%
    mutate(filename=my_jagscsv[csvi])%>%
    bind_rows(all_jags_res_df)
}
all_jags_res_df<-all_jags_res_df%>%
  mutate(src="JAGS",
         comp_num=as.integer(comp_num),
         ref_num=as.integer(ref_num),
         model=case_when(grepl("CUjagsres",filename)~"common\nUMIE",
                         grepl("IUjagsres",filename)~"ind\nUMIE",
                         grepl("Cjagsres",filename)~"common",
                         grepl("Ijagsres",filename)~"ind",
                         grepl("NMAjagsres",filename)~"NMA",
                         TRUE~"other"))%>%
  rename(rowname=rowname.x)%>%
  select(-X, rowname.y,graphname,  graphname2)%>%
  mutate_if(is.character,function(x) as.factor(x))
