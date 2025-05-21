# Simulate binary data with UMIE interactions
# Modified NAK to add interactions
# link: https://github.com/TEvrenoglou/PL-NMA/blob/main/Scenario%201-16/Scenario%201.R


library(dplyr)

dlong <- function(data){
  # prepare the data set for the long-format
  dat1 <- data %>%
    dplyr::select(studlab, t1, n1, r1,xcov) %>%
    mutate(treatment=LETTERS[t1],
           Ival_std=0)%>%
    rename("study" = "studlab", 
           "treat" = "t1", 
           "events" = "r1", 
           "n" = "n1")
  
  dat2 <- data %>%
    dplyr::select(studlab, t2, n2, r2,xcov) %>%
    mutate(treatment=LETTERS[t2],
           Ival_std=1)%>%
    rename("study" = "studlab", 
           "treat" = "t2", 
           "events" = "r2", 
           "n" = "n2")
  
  data <- rbind(dat1, dat2)
  data <- data %>%
    mutate(ID = paste(study, treat, sep ="")) %>%
    filter(!duplicated(ID)) %>%
    dplyr::select(-ID)
  
  data=data[order(data$study),]
  return(data)
}


# continue f(x) ###########################################
data1=list()
data2=list()
logOR=list()
OR=list()

### define treatment indices
t1=c()
t2=c()
for (i in 1:(NT-1)){
  for (k in (i+1):NT){
    for(j in 1:NS){
      t1=c(t1,i)
      t2=c(t2,k)      }}}
N.stud=length(t1)

### define patients per treatment arm
for (i in 1:N.sim){
  
  data1[[i]]=data.frame(t1,t2)
  data1[[i]]$studlab=c(1:(N.stud))
  data1[[i]]$n1=data1[[i]]$n2=round(runif(N.stud,Npmin,Npmax))
  # data1[[i]]$xcov=rbinom(n=N.stud,size=maxxrange,prob=0.5)# cat covariate
  data1[[i]]$xcov=runif(n = N.stud,min = minxrange,max=maxxrange)
  # data1[[i]]$xcov=rnorm(n = N.stud,mean = 0,sd=maxxrange) # If the observation gets too large this can lead to issues in the OR
  
  logOR[[i]]=es_lOR  
  OR[[i]]=exp(logOR[[i]])
  
  }


#### define probabilities per treatment, per study arm
for (i in 1:N.sim){  
  data1[[i]]$p.ref=runif(N.stud,minr,maxr) #### study-specific probability of an event in treatment 1 # modified to increase probability
  data1[[i]]$odds.ref=data1[[i]]$p.ref/(1-data1[[i]]$p.ref)
}

#### define probabilities per treatment, per study arm
#new# 
Sigma=matrix(c(tau^2,tau^2/2,tau^2/2,tau^2),nrow=2)

for (i in 1:N.sim){  
  for(j in 1:N.stud){
    test1=MASS::mvrnorm(n=1,mu=c(logOR[[i]][data1[[i]]$t1[j]],logOR[[i]][data1[[i]]$t2[j]]),Sigma=Sigma)
    data1[[i]]$st.sp.logOR.t1[j]=ifelse(data1[[i]]$t1[j]==1,0,
                                        test1[1])
    data1[[i]]$st.sp.logOR.t2[j]=ifelse(data1[[i]]$t1[j]==1,
                                        rnorm(n=1,mean=logOR[[i]][data1[[i]]$t2[j]],sd=tau),
                                        test1[2])
    data1[[i]]$odds.t1[j]=data1[[i]]$odds.ref[j]*exp(data1[[i]]$st.sp.logOR.t1[j]+0*data1[[i]]$xcov[j])
    data1[[i]]$odds.t2[j]=data1[[i]]$odds.ref[j]*exp(data1[[i]]$st.sp.logOR.t2[j]+beta_x[data1[[i]]$t1[j],data1[[i]]$t2[j]]*data1[[i]]$xcov[j])

	data1[[i]]$p.t1[j]=data1[[i]]$odds.t1[j]/(1+data1[[i]]$odds.t1[j])
    data1[[i]]$p.t2[j]=data1[[i]]$odds.t2[j]/(1+data1[[i]]$odds.t2[j])
  }}


#### generate the data
for (i in 1:N.sim){  
  for(j in 1:N.stud){
    data1[[i]]$r1[j]=rbinom(1,data1[[i]]$n1[j],data1[[i]]$p.t1[j]) 
    data1[[i]]$r2[j]=rbinom(1,data1[[i]]$n2[j],data1[[i]]$p.t2[j]) 
  }}

for (i in 1:N.sim){ 
  # data1[[i]]=data1[[i]][,-c(7:12)] # wide data
  data1[[i]]=data1[[i]][,c("t1","t2","studlab","xcov","n2","n1","r1","r2")] # wide data
  data2[[i]]=dlong(data1[[i]]) # long version
}

