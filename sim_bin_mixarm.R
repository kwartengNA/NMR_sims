# Combine 2 and multi-arm data binary ##################
# run simulations for two arm and multiarm #############
#2 arm studies 
source("sim_bin2arm.R")
arm2_data2<-data2
arm2_data1<-data1
source("sim_bin_multiarm.R")
armmulti_data2<-data2
armmulti_data1<-data1

# randomly define number of arms in study ##################################
if(NT>3){
  N.arms<-data.frame(study=1:(NS),
                     n_arms=rbinom(n=NS, size=NT-2, prob=c(0.843,0.132,rep((1-0.975)/(NT-4), NT-4)))+2) #nmadb
} else {
  N.arms<-data.frame(study=1:(NS),
                     n_arms=rbinom(n=NS, size=NT-2, prob=c(0.843,0.132))+2) #nmadb
}

# combine multi-arm and two-arm data ########################
# such that so that they fulfill the conditions of number of studies per comparison
rm(data1,data2)# just to ensure we only have 1 version of the final data
data2 <-NULL # retain results. this overwrites the other data2 results for compatibility with existing syntax

for(nsim in 1:N.sim){
  temp_armmulti_data2<-armmulti_data2[[nsim]]%>%left_join(N.arms, by=c("study"))%>%
    group_by(study, Ival_gen)%>%
    mutate(n_arms_temp=ifelse(Ival_gen==0,1,n_arms),
           keep_arms=case_when(Ival_gen==0~1,
                               Ival_gen==1~ifelse(row_number()%in% sample(x = n(), size=unique(n_arms_temp)-1),1,0)))%>%
    ungroup()%>%
    filter(keep_arms==1)%>%
    select(-c(n_arms_temp, treat.ref, n_arms, keep_arms))%>%
    mutate(studyname=paste("gen_multiarm", study, sep="_"))%>%
    group_by(study)%>%
    mutate(Ival_std=ifelse(treat==min(treat),0,1))%>%# we have to redefine the trt order Ival
    ungroup()
  
  # just to id unique pairs
  temp_multiarm<-pairwise(treat=treatment,n = n,event = events,studlab=study,
           data=temp_armmulti_data2%>%mutate(events=ifelse(events==n|(events==0),n/2,events)),
           sm="OR")
  
  expected<-as.data.frame(combn(x=LETTERS[1:NT],m=2))%>%
    pivot_longer(cols=everything())%>%
    group_by(name)%>%
    summarise(comp=paste(value, collapse=" "))%>%
    ungroup()%>%
    mutate(countn=3)
  
  is_missing<-temp_multiarm%>%
    mutate(comp=paste(treat1,treat2))%>%
    count(comp)%>%
    right_join(expected, by=c("comp"))%>%
    mutate(missing=countn-ifelse(is.na(n)==TRUE,0,n))%>%
    filter(missing>0)%>%
    select(comp, missing)
  
  
  data2[[nsim]]<-arm2_data2[[nsim]]%>%
    group_by(study)%>%
    summarise(comp=paste(treatment, collapse=" "))%>%
    ungroup()%>%
    right_join(is_missing, by=c("comp"))%>%
    group_by(comp)%>%
    filter(row_number()<=missing)%>%
    ungroup()%>%
    select(study)%>%
    left_join(arm2_data2[[nsim]], by=c("study"))%>%
    mutate(studyname=paste("gen_twoarm", study, sep="_"),
           Ival_gen=Ival_std)%>%
    bind_rows(temp_armmulti_data2)%>%
    mutate(study=as.numeric(as.factor(studyname)))
  
  rm(temp_multiarm,is_missing, expected)
}

# netgraph(netmeta(pairwise(treat=treatment,n = n,event = events,studlab=study,
#                           data=data2[[1]]%>%mutate(events=ifelse(events==n|(events==0),n/2,events)),
#                           sm="OR")),
#          plastic=FALSE, lwd=3, col="navyblue",
#          start.layout="circle",
#          cex=1.1,
#          offset=0.05,
#          points=TRUE,bg.points="#FAA613",cex.points=3,col.points="navyblue",
#          number.of.studies=TRUE,
#          pch.points=21, multiarm = TRUE, col.multiarm = "#FAA613")
#          
