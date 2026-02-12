# shiny transformations to standardize output and minimize file size of inputs

# read from a specific seed and output to shiny directory
# gendir<-getwd()
# seedval<-545

# transformations ######################
# graph res ################
graph_res<-readRDS(file=paste(gendir, paste("seed_", seedval,sep=""),"filter_param_res.rds", sep="/"))

graph_res2<-graph_res%>%
  filter(!(model=="NMA" & paramtype=="beta_x"))%>%
  mutate(tau=factor(tau,levels=c("tau_none", "tau_moderate", "tau_moderate2","tau_moderate3", "tau_high")),
         tau_txt=factor(tau_txt,levels=c("none", "moderate", "moderate2","moderate3", "high")),
         density=factor(density, levels=c("connect_all","connect_full",
                                          "connect_moderate_dense", "connect_moderate",  "connect_moderate_sparse",
                                          "connect_low")))%>%
  mutate(true_model=trimws(paste(ifelse(true_assumption=="independent","ind",true_assumption), ifelse(true_consistency==FALSE, "UMIE",""))))%>%
  mutate(tau=as.factor(ifelse(tau_txt=="moderate3","tau_high", as.character(tau))),
         tau_txt=as.factor(ifelse(tau_txt=="moderate3","high", as.character(tau_txt))),
         tau_txt=factor(tau_txt, levels = c("none", "moderate", "high")))%>%
  mutate(densitybkup=density,
         density=case_when(densitybkup=="connect_full"~100,
                           densitybkup=="connect_moderate"~50,
                           densitybkup=="connect_moderate_dense"~75,
                           densitybkup=="connect_moderate_sparse"~35,
                           densitybkup=="connect_low"~25,
                           TRUE~as.numeric(NA)),
         keep_prop_txt=factor(keep_prop_txt, levels=c("moderate_sparse", "moderate", "moderate_dense", "full")),
         comp=ifelse(comp=="NA NA", NA, comp),
         reference_treat=substr(comp,3,3),
         comparator_treat=substr(comp,1,1))%>%
  mutate(density=factor(density, 
                        levels = rev(sort(unique(density)))),
         tau=factor(tau, levels=rev(sort(unique(tau)))))%>%
  mutate(src=case_when(src=="JAGS"~"Bayesian", 
                       src=="RMA_MV"~"Frequentist",
                       TRUE~src),
         paramtype=case_when(paramtype=="es"~"(Adjusted)\nTreatment Effect",
                             paramtype=="beta_x"~"Treatment:Covariate\nInteraction Effect",
                             TRUE~paramtype),
         model=case_when(model=="ind"~"independent",
                         model%in%c("ind\nUMIE","ind\nUMIE\nmin")~"independent\nUMIE",
                         TRUE~model),
         model=paste0(toupper(substr(model,1,1)),substr(model,2,nchar(model))),
         model=factor(model,
                      levels=c("NMA", "Common", "Independent", "Independent\nUMIE", "Common\nUMIE")),
         heterogeneity_num=case_when(tau_txt == "none"     ~ 0,
                                     tau_txt == "moderate"  ~ 0.28, #doi:10.1093/ije/dys041. Semi Objective Non Pharmacological
                                     tau_txt == "high" ~0.5, #doi:10.1093/ije/dys041. Subjective Non Pharmacological
                                     TRUE~as.numeric(NA)))%>%
  # select(-ends_with("mcse"))%>%
  select(-starts_with("median"))%>%
  select(-c(seed, outcome, interaction, events, tau, densitybkup, consistency_txt, interaction_direction, event_rate, n_groups,mse, rmse, mean_mse, mean_rmse))%>%
  rename(comparison=comp, parameter=paramtype, mse=man_mse,mean_mse=mean_man_mse, implemented_model=model,tau=heterogeneity_num, 
         density_txt=keep_prop_txt)%>%
  select(scenarionum, true_model, true_arms, true_consistency, true_assumption, density_txt, density, tau_txt, tau,
         comparison, comparator_treat, reference_treat, src, parameter, implemented_model, everything())%>%
  arrange(scenarionum, #true_consistency, true_assumption, density_txt,tau
          comparison, comparator_treat, reference_treat, src, parameter, implemented_model)%>%
  rename(Source=src)

# diag res ############
diag_res<-readRDS(file=paste(gendir, paste("seed_", seedval,sep=""),"diag_res.rds", sep="/"))
diag_res2<-graph_res2%>%count(density, density_txt, fileloc, tau, tau_txt)%>%select(-n)%>%
  left_join(diag_res%>%select(-c(density, tau, tau_txt)))%>%
  mutate(true_model=ifelse(true_model=="ind UMIE", "ind UMIE min", true_model))%>%
  filter(!model%in%c("ind\nUMIE\nmax"))%>%
  mutate(src=case_when(src=="JAGS"~"Bayesian", 
                       src=="RMA_MV"~"Frequentist",
                       TRUE~src),
         model=case_when(model=="ind"~"independent",
                         model%in%c("ind\nUMIE","ind\nUMIE\nmin")~"independent\nUMIE",
                         TRUE~model),
         model=paste0(toupper(substr(model,1,1)),substr(model,2,nchar(model))),
         model=factor(model,
                      levels=c("NMA", "Common", "Independent", "Independent\nUMIE", "Common\nUMIE")),
  )%>%
  select(-c(seed, outcome, interaction, events, consistency_txt, interaction_direction, event_rate, keep_prop_txt, filename, keep_prop_txt, filenum))%>%
  select(scenarionum, true_model, true_arms, true_consistency, true_assumption,density_txt, density, tau_txt, tau, src,
         parameter=paramtype, implemented_model=model,simulation, ML, REML,
         everything())%>%
  rename(Source=src)


# alldata ########
gc()
alldata<-readRDS(file=paste(gendir, paste("seed_", seedval,sep=""),"all_data.rds", sep="/"))

alldata2<-alldata%>%
    as.data.frame()%>%
    mutate(fileloc=gsub("^.*(?=seed)", "",fileloc, perl=TRUE),
           fileloc=gsub("/data2_binary.rds", "",fileloc),
           simulation=as.integer(sim),
           master_scenarionum=as.numeric(as.factor(fileloc)),
           scenario_sim=paste(master_scenarionum, sim),
           studlab=paste(study, sim, master_scenarionum))%>%
    select(-c(density, tau, tau_txt))%>%
    right_join(graph_res2%>%select(fileloc, scenarionum, true_arms, true_consistency, true_assumption,density_txt, density, tau_txt, tau)%>%distinct())%>%
  select(-c(comparison, study.y, 
            seed, outcome, interaction, output, consistency_txt, interaction_direction, event_rate, keep_prop_txt, keep_prop_txt,study_reference_treat=treat.ref,
            study, master_scenarionum, sim))%>%
    select(scenarionum, true_arms, true_consistency, true_assumption,density_txt, density, tau_txt, tau,network, networkid, simulation, studlab, xcov, treatment, arm_directionality=Ival_std, 
           events,everything())
    
  
# write #####
readr::write_rds(x = graph_res2, file = paste(gendir, paste("shiny_results/seed_", seedval,sep=""),"graph_res.rds", sep="/"),"xz", compression = 9L)
readr::write_rds(x = diag_res2, file = paste(gendir, paste("shiny_results/seed_", seedval,sep=""),"diag_res.rds", sep="/"),"xz", compression = 9L)
readr::write_rds(x = alldata2, file = paste(gendir, paste("shiny_results/seed_", seedval,sep=""),"input_data.rds", sep="/"),"xz", compression = 9L)

