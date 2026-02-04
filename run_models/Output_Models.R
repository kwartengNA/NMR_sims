## combine jags res ########

outputdir<-paste(loopdir,"/",sep="")
myfiles<-list.files(outputdir,full.names = TRUE)
my_jagscsv<-myfiles[grepl("\\.rds", myfiles)==TRUE]
my_jagscsv<-my_jagscsv[grepl("jagsres", my_jagscsv)==TRUE]
my_jagscsv<-my_jagscsv[grepl(seedval, my_jagscsv)==TRUE]
my_jagscsv<-my_jagscsv[grepl(paste("_",true_assumption,"_",sep=""), my_jagscsv)==TRUE]
my_jagscsv<-my_jagscsv[as.numeric(grepl("UMIE", my_jagscsv))==((true_consistency-1)*-1)]
my_jagscsv<-my_jagscsv[grepl(true_arms, my_jagscsv)==TRUE]

all_jags_res_df<-NULL
if(length(my_jagscsv)>0){
for(csvi in 1:length(my_jagscsv)){
  all_jags_res_df<-readRDS(my_jagscsv[csvi])%>%#read.csv(my_jagscsv[csvi])
    mutate(filename=my_jagscsv[csvi])%>%
    bind_rows(all_jags_res_df)
}
all_jags_res_df<-all_jags_res_df%>%
  mutate(src="JAGS",
         comp_num=as.integer(comp_num),
         ref_num=as.integer(ref_num),
         model=case_when(grepl("stdCUjagsres",filename)~"common\nUMIE\nstd",
                         grepl("stdIUjagsres",filename)~"ind\nUMIE\nstd",
                         grepl("CUjagsres",filename)~"common\nUMIE",
                         grepl("IUjagsres",filename)~"ind\nUMIE",
                         grepl("Cjagsres",filename)~"common",
                         grepl("Ijagsres",filename)~"ind",
                         grepl("NMAjagsres",filename)~"NMA",
                         TRUE~"other"))%>%
  rename(rowname=rowname.x)%>%
  select(-c(rowname.y, seed, outcome, true_arms, density, interaction, events, output,
            true_consistency, true_assumption, keep_prop_txt, consistency_txt, tau, tau_txt, interaction_direction, event_rate))%>%#information is in fileloc
  mutate_if(is.character,function(x) as.factor(x))
}
## ---- combres----------------------------------------------------------------------------
if(exists("all_jags_res_df")==FALSE | is.null(all_jags_res_df)){
  all_jags_res_df<-NULL
  } else{
    all_jags_res_df<-all_jags_res_df%>%
      mutate(simulation=as.integer(simulation))
  }
if(exists("all_freq_res_df")==FALSE){
  all_freq_res_df<-NULL
} else {
  all_freq_res_df<-all_freq_res_df%>%
    mutate(simulation=as.integer(simulation),
           filename=paste0(loopdir,"/all_freq_res_df.rds"))
}
if(exists("std_freq_res_df")==FALSE){
  std_freq_res_df<-NULL
} else {
  std_freq_res_df<-std_freq_res_df%>%
    mutate(simulation=as.integer(simulation),
           filename=paste0(loopdir,"/std_freq_res_df.rds"))
}

all_res<-all_jags_res_df%>%
  bind_rows(all_freq_res_df)%>%
  bind_rows(std_freq_res_df)%>%
  mutate(bias2=bias^2)
rm(all_jags_res_df)
gc()

