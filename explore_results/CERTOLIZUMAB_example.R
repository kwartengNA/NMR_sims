options(ggplot2.discrete.colour= c("#FAA613","#000080", "#D7263D","#FF33CC","#1C77C3", "#4C4C4C","#4AA61C","#FBCEB1"),# expand to 10 colors
        ggplot2.discrete.fill= c("#FAA613","#000080",  "#D7263D","#FFA3AF","#1C77C3","#4C4C4C", "#4AA61C","#FBCEB1"))

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
# mydata_long$treatment<-LETTERS[mydata_long$treat]
mypw<-pairwise(treat = treat, event = events, n = n, studlab = ID,data=mydata_long)
mynet<-netmeta(mypw)
netgraph(mynet, plastic=FALSE, lwd=0.75,
         points = TRUE, pch = 21, cex.points = 3, col.points = "black",
         col = "#000080",
         bg.points = "#FAA613", thickness = "se.fixed")
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
  mutate(model="NMA", src="JAGS")%>%
  bind_rows(NMRi_res%>%mutate(model="ICIE", src="JAGS"))%>%
  bind_rows(NMRc_res%>%mutate(model="CCIE", src="JAGS"))%>%
  bind_rows(NMRcumie_res%>%mutate(model="CUMIE", src="JAGS"))%>%
  bind_rows(NMRiumie_res%>%mutate(model="IUMIE", src="JAGS"))%>%
  bind_rows(mynma$results%>%mutate(model="NMA", src="RMA_MV", type="d"))%>%
  bind_rows(mynmr_i$results%>%mutate(model="ICIE", src="RMA_MV"))%>%
  bind_rows(mynmr_c$results%>%mutate(model="CCIE", src="RMA_MV"))%>%
  bind_rows(mynmr_iumie$results%>%mutate(model="IUMIE", src="RMA_MV"))%>%
  bind_rows(mynmr_cumie$results%>%mutate(model="CUMIE", src="RMA_MV"))%>%
  mutate(lower=coalesce(lower, lci_NMR),
         upper=coalesce(upper, uci_NMR),
         coef=coalesce(coef, TE_NMR),
         se=coalesce(seTE_NMR, se),
         converged=ifelse(Rhat<=1.1|is.na(Rhat),TRUE,FALSE),
         comp_treat=gsub("X","",coalesce(treat,comp_num)),
         comp_treat=LETTERS[as.numeric(comp_treat)],
         ref_num=coalesce(ref_num, "1"),
         ref_treat=LETTERS[as.numeric(ref_num)],
         comp=ifelse(grepl("comp",treat)==TRUE, gsub("\\.comp_","",treat),
                     paste(ref_treat,comp_treat)),
         comp=gsub("_vs_"," ",comp),
         comp=ifelse(grepl("comp",treat)==TRUE,
                     paste(LETTERS[as.numeric(as.character(gsub(" .*","",comp)))],LETTERS[as.numeric(as.character(gsub(".* ","",comp)))] ),
                     comp),
         type=case_when(paramtype=="es"~"d",
                        paramtype=="beta_x"~"beta",
                        TRUE~coalesce(type, paramtype)))%>%
  filter(type%in%c("d","beta"))
tau_res<-NMA_res%>%
  mutate(model="NMA", src="JAGS")%>%
  bind_rows(NMRi_res%>%mutate(model="ICIE", src="JAGS"))%>%
  bind_rows(NMRc_res%>%mutate(model="CCIE", src="JAGS"))%>%
  bind_rows(NMRcumie_res%>%mutate(model="CUMIE", src="JAGS"))%>%
  bind_rows(NMRiumie_res%>%mutate(model="IUMIE", src="JAGS"))%>%
  bind_rows(data.frame(TE_NMR=sqrt(mynma$tau2))%>%mutate(model="NMA", src="RMA_MV", type="tau"))%>%
  bind_rows(data.frame(TE_NMR=sqrt(mynmr_i$tau2))%>%mutate(model="ICIE", src="RMA_MV", type="tau"))%>%
  bind_rows(data.frame(TE_NMR=sqrt(mynmr_c$tau2))%>%mutate(model="CCIE", src="RMA_MV", type="tau"))%>%
  bind_rows(data.frame(TE_NMR=sqrt(mynmr_iumie$tau2))%>%mutate(model="IUMIE", src="RMA_MV", type="tau"))%>%
  bind_rows(data.frame(TE_NMR=sqrt(mynmr_cumie$tau2))%>%mutate(model="CUMIE", src="RMA_MV", type="tau"))%>%
  mutate(lower=lci_NMR,
         upper=uci_NMR,
         coef=TE_NMR,
         se=seTE_NMR,
         converged=ifelse(Rhat<=1.1|is.na(Rhat),TRUE,FALSE),
         type=case_when(paramtype=="es"~"d",
                        paramtype=="beta_x"~"beta",
                        TRUE~coalesce(type, paramtype)))%>%
  filter(type%in%c("tau"))


rownames(tau_res)<-NULL
rownames(all_res)<-NULL
# graphs #########
all_res%>%
  filter((ref_num==1|is.na(ref_num)))%>%
  ggplot()+
  geom_pointrange(aes(x=coef, xmin=lower, xmax=upper, y=comp, color=src, shape=converged), position=position_dodge(width=0.5))+
  geom_vline(aes(xintercept = 0))+
  scale_shape_manual(values=c("TRUE"=19,"FALSE"=5))+
  facet_grid(type~model, scales="free")+
  theme_bw()

tau_res%>%
  rename(tau=coef)%>%
  ggplot()+
  geom_col(aes(y=tau, x=model,fill=src, alpha=converged))+
  geom_pointrange(aes(y=tau, ymin=lower, ymax=upper, x=model, 
                      color=src, 
                      shape=converged))+
  geom_hline(aes(yintercept = 0))+
  scale_alpha_manual(values=c("TRUE"=1,"FALSE"=0.5))+
  scale_shape_manual(values=c("TRUE"=19,"FALSE"=5))+
  scale_fill_manual(values=c("JAGS"=NA,"RMA_MV"="#FAA613"))+
  scale_color_manual(values=c("JAGS"="#000080","RMA_MV"=NA))+
  theme_bw()

tau_res%>%filter(model=="NMA")%>%select(NMA_tau=coef,src)%>%
  right_join(tau_res)%>%
  rename(tau=coef)%>%
  ggplot()+
  geom_errorbar(aes(y=tau, ymin=lower, ymax=upper, x=model, 
                    color=src))+
  geom_point(aes(y=tau, x=model,color=src, shape=converged), size=2)+
  geom_hline(aes(yintercept = 0))+
  geom_hline(aes(yintercept = NMA_tau, color=src),lty=2)+
  scale_shape_manual(values=c("TRUE"=19,"FALSE"=5))+
  theme_bw()

all_res%>%
  filter(model%in%c("mynma", "mynmr_c"))%>%
  select(src, model,comp_num, ref_num, Comparator, Reference, coef, se,lower, upper,Rhat)
