# Applied example EXAMPLE 5. CERTOLIZUMAB: CONTINUOUS COVARIATE PMID: 27905717
# libs and data ############
library(readxl)
library(dplyr)
library(jagsUI)
library(ggplot2)
library(netmeta)
source("netmetareg.R") 
mydata <- read_excel("data/CERTOLIZUMAB_data.xlsx")
mydata_long<-mydata%>%
  select(number_arms,ID,Study,treat=t_1, n=n_1, events=r_1,x)%>%
  bind_rows(mydata%>%
              select(number_arms,ID,Study,treat=t_2, n=n_2, events=r_2,x))
# Model files iterations, and burnin ####################
model_AD_NMR_BIN_RE_independent<-"JAGS/model_AD_NMR_BIN_RE_independent.txt"
model_AD_NMR_BIN_RE_const<-"JAGS/model_AD_NMR_BIN_RE_constant.txt"
model_AD_NMA_BIN_RE<-"JAGS/model_AD_NMA_BIN_RE.txt"
model_AD_NMR_BIN_RE_const_UME<-"JAGS/model_AD_NMR_BIN_RE_constant_inconsistency_Idir.txt"
model_AD_NMR_BIN_RE_independent_UME<-"JAGS/model_AD_NMR_BIN_RE_independent_inconsistency_Idir.txt"
n.iter = 30000
n.burnin = 10000

# example freq NMA ########
mypw<-pairwise(treat = treat, event = events, n = n, studlab = ID,data=mydata_long)
mynet<-netmeta(mypw)
netgraph(mynet)
mynma<-netmetareg(mynet, covar =NULL, consistency = TRUE,assumption = "i")
mynmr_i<-netmetareg(mynet, covar = x, consistency = TRUE,assumption = "i")
mynmr_c<-netmetareg(mynet, covar = x, consistency = TRUE,assumption = "c")
mynmr_iumie<-netmetareg(mynet, covar = x, consistency = FALSE,assumption = "i")
mynmr_cumie<-netmetareg(mynet, covar = x, consistency = FALSE,assumption = "c")
# JAGS equivalents ######################
# data transformation #######
  temp_mydata_wide<-
    mydata_long%>%
    rename(treat_num=treat)%>%
    arrange(ID,treat_num)%>%
    group_by(ID)%>%
    mutate(nstudy=sum(n),
           treat_order=row_number(),
           Ival_std=ifelse(treat_order==min(treat_order),0,1))%>%
    ungroup()%>%
    tidyr::pivot_wider(names_from = treat_order,
                       values_from =c("events","n","treat_num","Ival_std"))
  
data_filter <- list("ns"=length(unique(temp_mydata_wide$ID)),
                      "nt"=max(temp_mydata_wide%>%
                                 select(starts_with("treat_num")),na.rm=TRUE),
                      "na"=rowSums(is.na(select(temp_mydata_wide,
                                                starts_with("treat_num")))==FALSE),
                      "x"=temp_mydata_wide$x,
                      "mx"=0,#mean(temp_mydata_wide$xcov),
                      ## need a column for each arm
                      "t"=temp_mydata_wide%>%select(starts_with("treat_num")),
                      ## need a column for each arm
                      "r"=temp_mydata_wide%>%select(starts_with("events_")),
                      "n"=temp_mydata_wide%>%select(starts_with("n_")),
                      ## need a column for each arm
                      "I_dir"= temp_mydata_wide%>%select(starts_with("Ival_"))
                    )

    
# NMA ########
jags_nma<- jags(
    data = data_filter,
    inits = NULL,
    parameters.to.save = c("lor","sd"),
    n.chains = 2, n.iter = n.iter,
    n.burnin = n.burnin,DIC=F,
    model.file = model_AD_NMA_BIN_RE)
  
  
NMA_res<-as.data.frame(jags_nma$summary)%>%
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
           paramtype=case_when(paramtype=="lor"~"es",
                               paramtype=="betadiff"~"beta_x",
                               paramtype=="sd"~"tau",
                               TRUE~paramtype))

# NMR ICIE ##############
jags_nmri<- jags(
  data = data_filter,
  inits = NULL,
  parameters.to.save = c("lor", "lor0", "betadiff","sd"),
  n.chains = 2, n.iter = n.iter,
  n.burnin = n.burnin,DIC=F,
  model.file = model_AD_NMR_BIN_RE_independent)

NMRi_res<-as.data.frame(jags_nmri$summary)%>%
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
         paramtype=case_when(paramtype=="lor0"~"es",
                             paramtype=="betadiff"~"beta_x",
                             paramtype=="sd"~"tau",
                             TRUE~paramtype))

# NMR CCIE ########################
jags_nmrc<- jags(
      data = data_filter,
      inits = NULL,
      parameters.to.save = c("lor", "lor0","betadiff","sd"),
      n.chains = 2, n.iter = n.iter,
      n.burnin = n.burnin,DIC=F,
      model.file = model_AD_NMR_BIN_RE_const)
    
NMRc_res<-as.data.frame(jags_nmrc$summary)%>%
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
             paramtype=case_when(paramtype=="lor0"~"es",
                                 paramtype=="betadiff"~"beta_x",
                                 paramtype=="sd"~"tau",
                                 TRUE~paramtype))
    


# NMR CUMIE ########
jags_nmrcumie<- jags(
  data = data_filter,
  inits = NULL,
  parameters.to.save = c("LOR","lor0", "beta","tau"),
  n.chains = 2, n.iter = n.iter,
  n.burnin = n.burnin,DIC=F,
  model.file = model_AD_NMR_BIN_RE_const_UME)

NMRcumie_res<-as.data.frame(jags_nmrcumie$summary)%>%
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
         paramtype=case_when(paramtype=="lor0"~"es",
                             paramtype=="beta"~"beta_x",
                             TRUE~paramtype))
# NMR IUMIE ########
jags_nmriumie<- jags(
  data = data_filter,
  inits = NULL,
  parameters.to.save = c("LOR","lor0", "beta","tau"),
  n.chains = 2, n.iter = n.iter,
  n.burnin = n.burnin,DIC=F,
  model.file = model_AD_NMR_BIN_RE_independent_UME)

NMRiumie_res<-as.data.frame(jags_nmriumie$summary)%>%
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
         paramtype=case_when(paramtype=="lor0"~"es",
                             paramtype=="beta"~"beta_x",
                             TRUE~paramtype))
# list of res ######
all_res<-NMA_res%>%
  mutate(model="mynma", src="JAGS")%>%
  bind_rows(NMRi_res%>%mutate(model="mynmr_i", src="JAGS"))%>%
  bind_rows(NMRc_res%>%mutate(model="mynmr_c", src="JAGS"))%>%
  bind_rows(NMRcumie_res%>%mutate(model="mynmr_cumie", src="JAGS"))%>%
  bind_rows(NMRiumie_res%>%mutate(model="mynmr_iumie", src="JAGS"))%>%
  bind_rows(mynma$results%>%mutate(model="mynma", src="RMA_MV", type="d"))%>%
  bind_rows(mynmr_i$results%>%mutate(model="mynmr_i", src="RMA_MV"))%>%
  bind_rows(mynmr_c$results%>%mutate(model="mynmr_c", src="RMA_MV"))%>%
  bind_rows(mynmr_iumie$results%>%mutate(model="mynmr_iumie", src="RMA_MV"))%>%
  bind_rows(mynmr_cumie$results%>%mutate(model="mynmr_cumie", src="RMA_MV"))%>%
  mutate(lower=coalesce(lower, lci_NMR),
         upper=coalesce(upper, uci_NMR),
         coef=coalesce(coef, TE_NMR),
         se=coalesce(seTE_NMR, se),
         converged=ifelse(Rhat<=1.1|is.na(Rhat),TRUE,FALSE),
         treat=gsub("X","",coalesce(treat,comp_num)),
         ref_num=coalesce(ref_num, "1"),
         comp=ifelse(grepl("comp",treat)==TRUE, 
                     gsub("\\.comp_","",treat),
                     paste(ref_num,treat)),
         comp=gsub("_vs_"," ",comp),
         type=case_when(paramtype=="es"~"d",
                        paramtype=="beta_x"~"beta",
                        TRUE~coalesce(type, paramtype)))%>%
  filter(type%in%c("d","beta") & (ref_num==1|is.na(ref_num)))
rownames(all_res)<-NULL
all_res%>%
  ggplot()+
  geom_pointrange(aes(x=coef, xmin=lower, xmax=upper, y=comp, color=src, shape=converged), position=position_dodge(width=0.5))+
  geom_vline(aes(xintercept = 0))+
  scale_shape_manual(values=c("TRUE"=19,"FALSE"=5))+
  facet_grid(type~model, scales="free")+
  theme_bw()

all_res%>%
  filter(model%in%c("mynma", "mynmr_c"))%>%
  select(src, model,comp_num, ref_num, Comparator, Reference, coef, se,lower, upper,Rhat)
