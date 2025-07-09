## ---- datagen----------------------------------------------------------------------------
seedval<-ifelse(exists("seedval")==FALSE,sample(1:1000,size=1),seedval)
jagsseedval<-ifelse(exists("jagsseedval")==FALSE,sample(1:1000,size=1),jagsseedval)
set.seed(seedval)
N.sim=ifelse(exists("N.sim")==FALSE,2,N.sim) # number of simulations
NT= 7 #### number of treatments in the network
NS= 3  #### number of studies per comparison (NS*NT for multi-arm)
keep_prop<-ifelse(exists("keep_prop")==FALSE,1,keep_prop)# proportion of studies to keep
tau=ifelse(exists("tau")==FALSE,sqrt(runif(n=1,min=0,max=1)),tau)

Npmin=30 #### minimum number of patients per arm
Npmax=60 #### maximum number of patients per arm
es_MD= c(0,seq(from =1/(NT-1), to = 1, by = 1/(NT-1))) #(logOR) ### logOR across studies when xcov is 0 (for each comparator treatment vs 1 normally i would choose a random value)

beta_range=0.75
# common/constant/consistent
true_consistency=ifelse(exists("true_consistency")==FALSE,TRUE,true_consistency)
true_assumption=ifelse(exists("true_assumption")==FALSE,"independent",true_assumption)
interaction_direction=ifelse(exists("interaction_direction")==FALSE,NA,interaction_direction)

# common/constant/consistent
if(true_consistency==TRUE){
  if(true_assumption=="common"){
    beta_x_old<-c(0,rep(runif(n=1,min=-1*beta_range,max=beta_range),NT-1))# common/constant/consistent
  } else if(true_assumption=="independent"){
    beta_x_old<-c(0,runif(n=NT-1,min=-1*beta_range,max=beta_range))# independent/consistent
  } else if(true_assumption=="NMA"){
    beta_x_old<-rep(runif(n=1,min=-1*beta_range,max=beta_range),NT)# NMA
  }
  if(is.null(interaction_direction)|is.na(interaction_direction)){
    beta_x_old<-beta_x_old# unchanged
  }else if(interaction_direction=="positive"){
    beta_x_old<-abs(beta_x_old)
  } else if(interaction_direction=="negative"){
    beta_x_old<-abs(beta_x_old)*-1
  }
  
  beta_x<-outer(beta_x_old,beta_x_old, FUN="-")*-1
} else {
  if(true_assumption=="common"){
    beta_x<-matrix(nrow=NT,ncol=NT,rep(x=runif(n=1,min=-1*beta_range,max=beta_range),times=NT^2))# common
  } else if(true_assumption=="independent"){
    beta_x<-matrix(nrow=NT,ncol=NT,runif(min=-1*beta_range,max=beta_range, n=NT^2))# independent
  } else if(true_assumption=="NMA"){
    beta_x<-matrix(nrow=NT,ncol=NT,rep(x=0,times=NT^2))
  }
  beta_x[lower.tri(beta_x)]<-beta_x[upper.tri(beta_x)]*-1 # not required, but for simplicity
}
diag(beta_x)<-0
beta_x


maxxrange<-5
minxrange<--1
true_arms=ifelse(exists("true_arms")==FALSE,"twoarm",true_arms)
if(true_arms=="multiarm"){
  source("simulate_data/sim_cont_multiarm.R")
} else if(true_arms=="mixarm"){
  source("simulate_data/sim_cont_mixarm.R")
} else {
  source("simulate_data/sim_cont2arm.R")
}
if(exists("keep_studies")==FALSE){
  keep_studies<-NULL 
}

if(is.null(keep_studies)){
  alltwoarm_studyids<-data2[[1]]%>%count(study)%>%filter(n==2)
  allmultiarm<-data2[[1]]%>%count(study)%>%filter(n>2)
  keep_studies<-sample(alltwoarm_studyids$study,#all 2 arm study numbers
                       floor(keep_prop*length(unique(data2[[1]]$study)))-nrow(allmultiarm), #number of 2-arm studies to keep
                       replace = FALSE)
  keep_studies<-c(keep_studies,allmultiarm$study)
} # otherwise use the existing keep studies

for(i in 1:N.sim){
  data2[[i]]<-data2[[i]]%>%
    filter(study%in%c(keep_studies))
}




## ---- refdat new -----------------------------------------------------------------------------------
refdat_orig<-data.frame(rowname=c(paste("es_A",LETTERS[1:NT],sep=""),
                                  # paste("beta_x_A",LETTERS[1:NT],sep=""),
                                  "tau","seed"),
                        graphname=c(paste("treat_name",LETTERS[1:NT],sep=""),
                                    # paste("`",LETTERS[1:NT],":xcov`",sep=""),
                                    "tau","seed"),
                        graphname2=c(paste("treat_name",LETTERS[1:NT],sep=""),
                                     # paste("`",LETTERS[1:NT],"_vs_A:Ivalrand:xcov`",sep=""),
                                     "tau","seed"),
                        Estimate=c(es_MD,#beta_x,
                                   tau,seedval),
                        Comparator=c(rep(LETTERS[1:NT],1),NA,NA),
                        Reference=c(rep("A",(1*NT)),NA,NA),
                        paramtype=c(rep("es",NT),#rep("beta_x",NT),
                                    "tau","seed"))

refdat<-refdat_orig
# all comparisons (UMIE models)
for(i in 1:(NT)){
  refdat<-refdat_orig%>%
    group_by(paramtype)%>%
    mutate(Comparator=Comparator,
           Reference=lag(Comparator,i),
           rowname=paste(paramtype,"_",Reference,Comparator,sep=""),
           paramtype=paramtype,
           Estimate=case_when(paramtype=="es"~Estimate-lag(Estimate,i),
                              TRUE~Estimate),
           graphname2=case_when(paramtype=="beta_x"~
                                  paste("`",Comparator,"_vs_",Reference,":Ivalrand:xcov`",sep=""),
                                paramtype=="es"~paste("treat_name",Comparator,sep=""),
                                TRUE~graphname))%>%
    ungroup()%>%
    filter(is.na(Reference)==FALSE)%>%
    bind_rows(refdat)
}

# reverse comparisons
refdat<-refdat%>%
  mutate(rowname=ifelse(!paramtype%in%c("tau","seed"),
                        paste(paramtype,"_",Comparator,Reference,sep=""),rowname),
         Estimate=ifelse(!paramtype%in%c("tau","seed","beta_x"),Estimate*-1,Estimate),
         Comparatornew=Reference,
         Reference=Comparator,
         graphname=case_when(paramtype=="beta_x"~
                               paste("`",Comparatornew,":xcov`",sep=""),
                             paramtype=="es"~paste("treat_name",Comparatornew,sep=""),
                             TRUE~graphname),
         graphname2=case_when(paramtype=="beta_x"~
                                paste("`",Comparatornew,"_vs_",Reference,":Ivalrand:xcov`",sep=""),
                              paramtype=="es"~paste("treat_name",Comparatornew,sep=""),
                              TRUE~graphname))%>%
  select(-Comparator)%>%
  rename(Comparator=Comparatornew)%>%
  bind_rows(refdat)%>%
  distinct()%>%
  mutate(comp=paste(Comparator, Reference))

refdat<-as.data.frame.table(beta_x)%>%
  rename(Estimate=Freq,
         Reference=Var1,
         Comparator=Var2)%>%
  mutate(rowname=paste(Comparator,Reference),
         graphname=paste("`",Comparator,":xcov`",sep=""),
         graphname2=paste("`",Comparator,"_vs_",Reference,":Ivalrand:xcov`",sep=""),,
         paramtype="beta_x",
         comp=paste(Comparator,Reference))%>%
  bind_rows(refdat)