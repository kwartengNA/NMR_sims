# Simulate continuous data with consistent interactions. multi arm
# Modified NAK to add interactions

library(dplyr)
library(netmeta)

data1=list()
data2=list()

### define treatment indices
t=(rep(c(1:NT),NS)) # all treatments occur in all studies
stud=rep(1:NS, each=NT)
for (i in 1:N.sim){ 
  data1[[i]]=data.frame(stud,t)
}
### define patients per treatment arm
for (i in 1:N.sim){   
  data1[[i]]$n=rep(round(runif(NS,Npmin,Npmax)),each=NT)
  data1[[i]]$xcov=rep(runif(NS,minxrange,maxxrange),each=NT)
}

narms=length(data1[[1]]$stud)

#### define reference treatment for each study. not effective
for (i in 1:N.sim){  
  data1[[i]]$treat.ref=rep(sample(size = NS,x=1:NT,replace=TRUE),each=NT)
}

Sigma=matrix(rep((tau^2)/2,times=NT^2),nrow=NT)
diag(Sigma)<-tau^2

es_MD_matrix<-outer(es_MD,es_MD, FUN="-")*-1

for (i in 1:N.sim){   

  for(j in 1:narms){
    
    test2=MASS::mvrnorm(n=1,mu=es_MD_matrix[data1[[i]]$treat.ref[j],],Sigma=Sigma)
    
    data1[[i]]$st.sp.MD.t[j]=ifelse(data1[[i]]$t[j]==data1[[i]]$treat.ref[j],
                                       0,
                                       test2[data1[[i]]$t[j]])
    data1[[i]]$yi.t[j]=data1[[i]]$odds.ref[j]*exp(data1[[i]]$st.sp.logOR.t[j]+ beta_x[data1[[i]]$treat.ref[j],data1[[i]]$t[j]]*data1[[i]]$xcov[j]) #
    data1[[i]]$sd.t[j]=abs(rnorm(n=1,mean=1,sd=.1))
  }
  }

##########

for (i in 1:N.sim){
  data2[[i]]<-data1[[i]][,c("stud","t","n","st.sp.MD.t","sd.t","xcov","treat.ref")]
  names(data2[[i]])<-c("study","treat","n","yi","sd","xcov","treat.ref")
  data2[[i]]$treatment<-LETTERS[data2[[i]]$treat]
  
  data2[[i]]$Ival_std<-ifelse(data2[[i]]$treat==1,0,1)
  data2[[i]]$Ival_gen<-ifelse(data2[[i]]$treat==data2[[i]]$treat.ref,0,1)
  }
