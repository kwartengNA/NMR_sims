# Simulate continuous data with interactions
# Modified NAK to add interactions
# es_MD<-es_lOR

library(dplyr)

dlong <- function(data){
  # prepare the data set for the long-format
  dat1 <- data %>%
    dplyr::select(studlab, t1, n1, yi.t1,xcov, sd=sd.t2) %>%
    mutate(treat_num=t1,
           treatment=LETTERS[t1],
           Ival_std=0)%>%
    rename("study" = "studlab", 
           "treat" = "t1", 
           "yi" = "yi.t1", 
           "n" = "n1")
  
  dat2 <- data %>%
    dplyr::select(studlab, t2, n2, yi.t2,xcov, sd=sd.t2) %>%
    mutate(treat_num=t2,,
           treatment=LETTERS[t2],
           Ival_std=1)%>%
    rename(
           "study" = "studlab", 
           "treat" = "t2", 
           "yi" = "yi.t2", 
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
MD=list()
# OR=list()

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
  
  MD[[i]]=es_MD  
  
  }


#new# 
Sigma=matrix(c(tau^2,tau^2/2,tau^2/2,tau^2),nrow=2)

for (i in 1:N.sim){  
  for(j in 1:N.stud){
    test1=MASS::mvrnorm(n=1,mu=c(MD[[i]][data1[[i]]$t1[j]],MD[[i]][data1[[i]]$t2[j]]),Sigma=Sigma)
    data1[[i]]$st.sp.MD.t1[j]=ifelse(data1[[i]]$t1[j]==1,0,
                                        test1[1])
    data1[[i]]$st.sp.MD.t2[j]=ifelse(data1[[i]]$t1[j]==1,
                                        rnorm(n=1,mean=MD[[i]][data1[[i]]$t2[j]],sd=tau),
                                        test1[2])
    data1[[i]]$yi.t1[j]=data1[[i]]$st.sp.MD.t1[j]+0*data1[[i]]$xcov[j]
    data1[[i]]$yi.t2[j]=data1[[i]]$st.sp.MD.t2[j]+beta_x[data1[[i]]$t1[j],data1[[i]]$t2[j]]*data1[[i]]$xcov[j]
    data1[[i]]$sd.t1[j]=abs(rnorm(n=1,mean=1,sd=.1))
    data1[[i]]$sd.t2[j]=abs(rnorm(n=1,mean=1,sd=.1))
    
  }}

for (i in 1:N.sim){ 
  
  data2[[i]]=dlong(data1[[i]]) # long version
}