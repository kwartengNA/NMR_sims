library(dplyr)
library(tidyr)
mydata <- readRDS("~/NMR_simulation_project/NMR_sims/seed_545/twoarm_consistency_all_resnew.rds")
all_conf_mat<-NULL
myfiles<-list.files(path = "C:/Users/kwarteng/Documents/NMR_simulation_project/NMR_sims/seed_545",full.names = TRUE)%>%as.data.frame()%>%rename(fileloc=1)%>%filter(grepl("all_res",fileloc))
for(i in myfiles$fileloc){
  mydata<-readRDS(i)

mydata<-mydata%>%
  filter(paramtype%in%c("tau", "es", "beta_x"))%>%
  filter(Comparator!=Reference)%>%
  filter(is.na(TE_NMR)==FALSE)%>%#failed est
  select(filename,model, Comparator, Reference, paramtype, es_NMR,es_logOR, TE_NMR, lci_NMR, uci_NMR, src)%>%
  mutate(es_NMR=coalesce(es_logOR, es_NMR),
         zero_res=ifelse(lci_NMR<=0 & uci_NMR>=0,1,0),
         zero_es=ifelse(es_NMR==0,1,0))

gc()

conf_mat<-mydata%>%
  filter(#model%in%c("ind","common") & src=="JAGS" & 
           paramtype=="beta_x")%>%
  select(filename,src, model, paramtype,zero_res,zero_es, lci_NMR,TE_NMR,uci_NMR, es_NMR)%>%
  count(filename, src, model, paramtype,zero_res, zero_es)%>%
  mutate(
    type = case_when(
      zero_res == 1 & zero_es == 1 ~ "TN",
      zero_res == 1 & zero_es == 0 ~ "FN",
      zero_res == 0 & zero_es == 1 ~ "FP",
      zero_res == 0 & zero_es == 0 ~ "TP"
    )
  ) %>%
  group_by(filename,src, model, paramtype,type) %>%
  summarise(n = sum(n), .groups = "drop")%>%
  ungroup()%>%
  group_by(filename,src, model, paramtype) %>%
  complete(
    type = c("TP", "FP", "FN", "TN"),
    fill = list(n = 0)
  ) %>%
  ungroup()

all_conf_mat<-all_conf_mat%>%bind_rows(conf_mat)
}
###################################################################
filenames<-all_conf_mat%>%select(filename)%>%
  distinct()%>%
  separate(col = filename, into=paste0("val",1:15), remove = FALSE, sep = "/")%>%
  rename(true_arms="val8", true_consistency=val9, true_assumption=val10, true_tau=val11, true_density=val12)

all_conf_mat %>%
  left_join(filenames%>%select(filename, contains("true")))%>%
  group_by(true_arms, #true_consistency, true_assumption, 
           true_tau, true_density,
    src, model, paramtype, type)%>%
  summarise(n=sum(n))%>%
  ungroup()%>%
  tidyr::pivot_wider(names_from = type, values_from = n)
  
t1<-all_conf_mat %>%
  left_join(filenames%>%select(filename, contains("true")))%>%
  filter(true_arms!="multiarm")%>%
  group_by(#true_arms, #true_consistency, true_assumption,true_tau, 
           true_density,
           src, model, paramtype, type)%>%
  summarise(n=sum(n))%>%
  ungroup()%>%
  tidyr::pivot_wider(names_from = type, values_from = n) %>%
  group_by(#true_arms, #true_consistency, true_assumption,true_tau, 
           true_density,
           src, model, paramtype)%>%
  summarise(
    accuracy    = (TP + TN) / (TP + TN + FP + FN),
    sensitivity = TP / (TP + FN), # also power
    specificity = TN / (TN + FP),    
    ppv         = TP / (TP + FP)
  )%>%
  ungroup()

t1%>%
  mutate(model=gsub("\n"," ",model),
         model=gsub(" min","",model),
         true_density=gsub("connect_","", true_density),
         true_density=factor(true_density, levels=c("full","moderate_dense","moderate","moderate_sparse","low")))%>%
  filter(!model%in%c("common_old", "ind UMIE max"))%>%
  filter(true_density!="low")%>%
  ggplot()+
  geom_tile(aes(x=true_density, y=interaction( src,model), fill=sensitivity))+
  theme_minimal()
