# Simulate binary data with consistent interactions. multi arm
# Modified NAK to add interactions
# link: https://github.com/TEvrenoglou/PL-NMA/blob/main/Scenario%2017-32/Scenario%2020.R


data1=list()
data2=list()

### define treatment indices
NS_all=NS#*NT
t=(rep(c(1:NT),NS_all)) # all treatments occur in all studies
stud=rep(1:NS_all, each=NT)
for (i in 1:N.sim){ 
  data1[[i]]=data.frame(stud,t)
}
### define patients per treatment arm
for (i in 1:N.sim){   
  data1[[i]]$n=rep(round(runif(NS_all,Npmin,Npmax)),each=NT)
  data1[[i]]$xcov=rep(runif(NS_all,minxrange,maxxrange),each=NT)
}

narms=length(data1[[1]]$stud)

#### define odds per treatment, per study arm
for (i in 1:N.sim){  
  data1[[i]]$p.ref=rep(runif(NS_all,minr,maxr),each=NT)
  data1[[i]]$odds.ref=data1[[i]]$p.ref/(1-data1[[i]]$p.ref)
  # data1[[i]]$treat.ref=rep(sample(size = NS_all,x=1:NT,replace=TRUE),each=NT)
  data1[[i]]$treat.ref=rep(sample(size = NS_all,x=1:NT,replace=TRUE),each=NT)
}

Sigma=matrix(rep((tau^2)/2,times=NT^2),nrow=NT)
diag(Sigma)<-tau^2

es_lOR_matrix<-outer(es_lOR,es_lOR, FUN="-")*-1

for (i in 1:N.sim){   

  for(j in 1:narms){
    
    test2=MASS::mvrnorm(n=1,mu=es_lOR_matrix[data1[[i]]$treat.ref[j],],Sigma=Sigma)
    
    data1[[i]]$st.sp.logOR.t[j]=ifelse(data1[[i]]$t[j]==data1[[i]]$treat.ref[j],
                                       0,
                                       test2[data1[[i]]$t[j]])
    data1[[i]]$odds.t[j]=data1[[i]]$odds.ref[j]*exp(data1[[i]]$st.sp.logOR.t[j]+ beta_x[data1[[i]]$treat.ref[j],data1[[i]]$t[j]]*data1[[i]]$xcov[j]) #
    data1[[i]]$p.t[j]=data1[[i]]$odds.t[j]/(1+data1[[i]]$odds.t[j])
  }
  }


#### generate the events

for (i in 1:N.sim){  
  for(j in 1:narms){
    data1[[i]]$r[j]=rbinom(n=1,size=data1[[i]]$n[j],prob=data1[[i]]$p.t[j])
  }}

###########

for (i in 1:N.sim){
  data2[[i]]<-data1[[i]][,c("stud","t","n","r","xcov","treat.ref")]
  names(data2[[i]])<-c("study","treat","n","events","xcov","treat.ref")
  data2[[i]]$treatment<-LETTERS[data2[[i]]$treat]
  
  data2[[i]]$Ival_std<-ifelse(data2[[i]]$treat==1,0,1)
  data2[[i]]$Ival_gen<-ifelse(data2[[i]]$treat==data2[[i]]$treat.ref,0,1)
  }

##################
