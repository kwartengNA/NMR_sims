# gendir<-getwd()
# seedval <- 120# seed value

# import filenames ########
outputfiles<-
  list.files(paste(gendir, paste("seed_", seedval,sep=""), sep="/"),#paste(gendir, paste("seed_", seedval, "_multi",sep=""), sep="/"),
             full.names = TRUE, recursive = TRUE)%>%
  lazy_dt()%>%
  as.data.frame()%>%
  rename(val=1)%>%
  mutate(relative_loc=gsub(gendir,"",val))%>%
  tidyr::separate(col=relative_loc, into = paste("vals", 1:11, sep=""), sep = "/", remove = FALSE)%>%
  rename(seed=vals2, outcome=vals3, true_arms=vals4, consistency_txt=vals5, true_assumption =vals6, tau =vals7, density=vals8, interaction=vals9, events =vals10, output=vals11)%>%
  mutate(output=ifelse(is.na(output) & outcome=="continuous" & grepl("_events", interaction), events, output),
         events=ifelse(outcome=="continuous" & grepl("_events", interaction), "NA", events),
         true_consistency=ifelse(consistency_txt=="consistency", TRUE, FALSE),
         tau_txt=gsub("tau_","",tau), 
         keep_prop_txt=gsub("connect_","",density),
         interaction_direction=gsub("_beta","", interaction),
         event_rate=gsub("_events", "", events),
         event_rate=ifelse(event_rate=="NA", NA, event_rate))%>%
  select(-starts_with("vals"))
# reconstruct/import network data #####################
data2files<-outputfiles%>%  filter(grepl("data2_.*rds", output))#%>%filter(grepl("csv", output))
starttime<-Sys.time()
alldata2<-vector("list", nrow(data2files))
for(i in 1:nrow(data2files)){
  tempdata<-readRDS(data2files[i,"val"])%>%bind_cols(data2files[i,])%>%rename(fileloc=val)#%>%filter(sim<=20)
  if(data2files[i,"outcome"]=="binary" & data2files[i,"true_arms"]=="twoarm"){
    tempdata<-tempdata%>%rename(r=events...5, events=events...19)
  } else  if(data2files[i,"outcome"]=="binary" & data2files[i,"true_arms"]%in%c("mixarm")){
    # gs data generation adds columns
    tempdata<-tempdata%>%rename(r=events...7, events=events...22)
  } else  if(data2files[i,"outcome"]=="binary" & !data2files[i,"true_arms"]%in%c("twoarm", "mixarm")){
    tempdata<-tempdata%>%rename(r=events...5, events=events...20)
  }
  
  alldata2[[i]] <- tempdata
  rm(tempdata)
  gc()
}
alldata2<-bind_rows(alldata2)

# assign network IDs for unique geometries
alldata2<-alldata2%>%
  group_by(fileloc, sim, study)%>%
  summarise(treatment=paste(treatment, collapse =";"))%>%
  ungroup()%>%
  count(fileloc, sim, treatment)%>%
  group_by(sim,fileloc)%>%
  summarise(network=paste(treatment, n, sep=",",collapse=";"))%>% # which unique network
  ungroup()%>%
  mutate(networkid=as.numeric(as.factor(network)))%>%
  right_join(alldata2, by=c("fileloc", "sim"))

alldata2<-alldata2%>%
  group_by(fileloc, sim, study)%>%
  mutate(Ival_order=ifelse(treat==min(treat, na.rm=TRUE),0,1))%>%
  ungroup()
gc()

# reconstruct/import reference data #####################
refdatfiles<-outputfiles%>%  filter(grepl("refdat_.*rds", output))#%>%filter(grepl("csv", output))
starttime<-Sys.time()
all_refdat<-vector("list", nrow(refdatfiles))
for(i in 1:nrow(refdatfiles)){
  tempdata<-readRDS(refdatfiles[i,"val"])%>%bind_cols(refdatfiles[i,])
  
  all_refdat[[i]] <- tempdata
  rm(tempdata)
  gc()
}
all_refdat<-bind_rows(all_refdat)


# export data ##################
readr::write_rds(x=alldata2, file=paste(gendir, paste("seed_", seedval,sep=""),"all_data.rds", sep="/"), compress="xz",compression = 9L)
readr::write_rds(x=all_refdat, file=paste(gendir, paste("seed_", seedval,sep=""),"all_refdat.rds", sep="/"), compress="xz",compression = 9L)
