
options(ggplot2.discrete.colour= c("#FAA613","#000080", "#D7263D","#FF33CC","#1C77C3", "#4C4C4C","#4AA61C","#FBCEB1"),# expand to 10 colors
        ggplot2.discrete.fill= c("#FAA613","#000080",  "#D7263D","#FFA3AF","#1C77C3","#4C4C4C", "#4AA61C","#FBCEB1"))


# libs ##############
# remove parameter grouping
library(shiny)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(devEMF)
# import data ###############
# setwd("C:/Users/kwarteng/Documents/NMR_sims")
gendir<-getwd()
seedval<-545
N.sim<-1000
outputdir<-paste(gendir, 
                 "man",
                 paste("seed_", seedval,sep=""), 
                 sep="/")
graph_res<-readRDS(file=paste(gendir, paste("seed_", seedval,sep=""),"filter_param_res.rds", sep="/"))%>%
  filter(!(model=="NMA" & paramtype=="beta_x"))
# graph_res<-readRDS(file=paste(gendir, paste("seed_", seedval,sep=""),"param_res.rds", sep="/"))%>%
#   filter(!(model=="NMA" & paramtype=="beta_x"))%>%
#   mutate(true_model=trimws(paste(ifelse(true_assumption=="independent","ind",true_assumption), ifelse(true_consistency==FALSE, "UMIE",""))))

diag_res<-readRDS(file=paste(gendir, paste("seed_", seedval,sep=""),"diag_res.rds", sep="/"))
diag_res<-diag_res%>%filter(model!="ind\nUMIE\nmax")
alldata2<-readRDS(file=paste(gendir, paste("seed_", seedval,sep=""),"all_data.rds", sep="/"))
alldata2<-alldata2%>%
  as.data.frame()%>%
  mutate(simulation=as.integer(sim),
         master_scenarionum=as.numeric(as.factor(fileloc)),
         scenario_sim=paste(master_scenarionum, sim),
         studyidnew=paste(study, sim, master_scenarionum))%>%
  right_join(graph_res%>%
               select(outcome, consistency_txt, true_assumption, true_arms, tau, density, scenarionum)%>%
               distinct())
N.sim<- max(alldata2$simulation)
gc()


## update labels ############
graph_res2<-graph_res%>%
  filter(keep_prop_txt!="low")%>%
  filter(grepl("std", model)==FALSE)%>%
  filter(model!="ind\nUMIE\nmax")%>%
  filter(model!="common_old")%>% # uniform heterogeneity prior
  # filter(true_arms=="mixarm")%>%
  # filter(true_arms=="twoarm")%>%
  mutate(true_arms=factor(true_arms, levels=c("twoarm", "mixarm", "multiarm")),
         RMSE=mean_rmse,
         paramtype=case_when(paramtype=="es"~"(Adjusted)\nTreatment Effect",
                             paramtype=="beta_x"~"Treatment:Covariate\nInteraction Effect",
                             TRUE~paramtype),
         model=case_when(model=="ind"~"ICIE",
                         model%in%c("ind\nUMIE","ind\nUMIE\nmin")~"IUMIE",
                         model%in%c("common")~"CCIE",
                         model%in%c("common\nUMIE")~"CUMIE",
                         TRUE~model),
         model=paste0(toupper(substr(model,1,1)),substr(model,2,nchar(model))),
         model=factor(model,
                      levels=c("NMA", "CCIE", "ICIE", "IUMIE", "CUMIE")),
         true_model=trimws(paste(ifelse(true_assumption=="independent","ind",true_assumption), ifelse(true_consistency==FALSE, "UMIE",""))),
         true_model  =case_when(true_model=="ind"~"True Model:\nICIE",
                                true_model=="ind UMIE"~"True Model:\nIUMIE",
                                true_model=="common UMIE"~"True Model:\nCUMIE",
                                true_model=="common"~"True Model:\nCCIE",
                                TRUE~paste0("True Model:\n",true_model)),
         true_model  =factor(true_model,
                             levels=c("True Model:\nNMA", "True Model:\nCCIE", "True Model:\nICIE", "True Model:\nIUMIE", "True Model:\nCUMIE")),
         keep_prop_num=case_when(keep_prop_txt == "all"      ~ 1,
                                 keep_prop_txt == "full"      ~ 1,
                                 keep_prop_txt == "moderate" ~ 0.5,
                                 keep_prop_txt == "low"      ~ 0.25,
                                 keep_prop_txt == "moderate_dense" ~ 0.75,
                                 keep_prop_txt == "moderate_sparse" ~ 0.35,
                                 TRUE~as.numeric(NA)),
         heterogeneity_num=case_when(tau_txt == "none"     ~ 0,
                                     tau_txt == "moderate" & true_arms=="mixarm" ~ 0.28, #doi:10.1093/ije/dys041. Semi Objective Non Pharmacological
                                     tau_txt == "moderate" & true_arms=="twoarm" ~ 0.28,# 0.25 was original default
                                     tau_txt == "moderate2" ~0.36, #doi:10.1093/ije/dys041. Subjective Non Pharmacological
                                     tau_txt == "moderate3" ~0.5, #doi:10.1093/ije/dys041. Subjective Non Pharmacological
                                     # tau_txt == "high" & true_arms=="twoarm" ~ 0.5, # 0.5 was original default. former high definition because 1 is not likely
                                     tau_txt == "high" & true_arms=="mixarm" ~ 1, # Seide
                                     TRUE~as.numeric(NA)),
         keep_prop_txt=factor(keep_prop_txt, levels=c("all", "full","moderate","low","moderate_dense","moderate_sparse" )),
         scenario=paste(paste(#"density:",
           keep_prop_num),
           paste(#"heterogeneity:",
             heterogeneity_num),  sep="\n"),
         scenario=factor(scenario, 
                         levels=c("1\n0",
                                  "1\n0.28",
                                  "1\n0.25",
                                  "1\n0.5",
                                  "1\n1",
                                  
                                  "0.75\n0",
                                  "0.75\n0.28",
                                  "0.75\n0.25",
                                  "0.75\n0.5",
                                  "0.75\n1",
                                  
                                  "0.5\n0",
                                  "0.5\n0.28",
                                  "0.5\n0.25",
                                  "0.5\n0.5",
                                  "0.5\n1",
                                  
                                  
                                  "0.35\n0",
                                  "0.35\n0.28",
                                  "0.35\n0.25",
                                  "0.35\n0.5",
                                  "0.35\n1",
                                  
                                  "0.25\n0",
                                  "0.25\n0.28",
                                  "0.25\n0.25",
                                  "0.25\n0.5",
                                  "0.25\n1")))%>%#mutate(model=paste(model,src))
  filter(src=="RMA_MV")

# Netgraph with export ################################

direct_evidence_col<-data.frame(sumdir=1:4,
                                sumdircol=factor(1:4,levels=c(1:4),labels=c(  "#89CFF0","#627DEC","#281CA6","#090057")))#,"#090057"
unique_networks<-alldata2%>%
  filter(keep_prop_txt!="low")%>%
  filter(true_arms=="mixarm")%>%
  group_by(networkid)%>%
  mutate(num_scenarios=length(unique(scenarionum)),
         which_scenarios=paste(unique(scenarionum), collapse=";"))%>%
  filter(scenario_sim==min(scenario_sim))%>%
  ungroup()


emf(file = paste(outputdir,"Networks.emf",sep="/"), height = 4,width = 15, emfPlus = FALSE )

plot.new()
n_facets <- length(unique(unique_networks$networkid))+1
facet_col <- n_facets#ceiling(sqrt(n_facets))
facet_row <- ceiling(n_facets / facet_col)
par(mfrow=c(facet_row,facet_col),
    mar = c(1, 1, 1, 1))

for(i in unique(unique_networks$networkid)){
  tmp_pairwise<-meta::pairwise(treat=treatment,n = n,event = events,studlab=studyidnew,data=unique_networks%>%filter(networkid==i),sm="OR")%>% # binary
    group_by(study)%>% # to get contributing direct comparisons in UMIE
    mutate(treat.ref=min(c(treat.orig1, treat.orig2)))%>%
    mutate(considered_direct=ifelse(treat.orig1==treat.ref|treat.orig2==treat.ref,TRUE,FALSE))%>%
    ungroup()%>%
    group_by(treat1,treat2)%>%
    mutate(sumdir=sum(considered_direct))%>%
    ungroup()
  
  mynet<-netmeta::netmeta(TE = TE,seTE = seTE,treat1 = treat1,treat2 = treat2,studlab = studlab, data=tmp_pairwise, method.tau = "REML", reference.group = "A")
  netmeta::netgraph(mynet, plastic=FALSE, 
                    lwd=3,
                    points=TRUE,
                    cex.points=3,
                    pch.points=21,
                    col.points="#090057",
                    bg.points="#FBE2B3",
                    col=tmp_pairwise%>%select(treat1,treat2,sumdir)%>%distinct()%>%arrange(treat1, treat2)%>%mutate(sumdir=as.character(factor(sumdir, levels=direct_evidence_col$sumdir, labels=direct_evidence_col$sumdircol)))%>%pull(sumdir),
                    main=paste(unique(paste(toupper(substr(tmp_pairwise$keep_prop_txt, 1, 1)), substr(tmp_pairwise$keep_prop_txt, 2, nchar(tmp_pairwise$keep_prop_txt)), sep="")), "Density"), 
                    col.multiarm="lightgrey", multiarm=TRUE)
  
}
plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), xaxt = "n", yaxt = "n", bty="n")
title(main="Reference Containing\nComparisons (UMIE)", line=-1)
legend(x="center",
       legend = direct_evidence_col$sumdir,
       # col = "black",
       fill = as.character(direct_evidence_col$sumdircol),
       # pch = 24,
       bty = "n",
       pt.cex = 1,
       text.col = "black",
       horiz = F ,
       inset = c(0.4,0))
dev.off()

## biasgraph trends in density (heterogeneity is constant) ##########

greylines<-graph_res2 %>%
  mutate(x = interaction((keep_prop_num),heterogeneity_num, true_arms, sep="\n"))%>%
  distinct(x)%>%
  filter(grepl("0.35", x) & grepl("0\\\n",x) & 
           grepl("multiarm",x )==FALSE)

# how she anticipates it ##########
coverage_graph1<-graph_res2%>%
  filter(true_arms!="multiarm")%>%
  filter(paramtype!="tau")%>%
  mutate(x = interaction((keep_prop_num),heterogeneity_num, sep="\n"))%>%
  mutate(paramtype=factor(paramtype, levels=c("(Adjusted)\nTreatment Effect", "Treatment:Covariate\nInteraction Effect", "tau")))%>%
  ggplot(aes(x = x,y = prop_coverage, color = model, group = interaction(heterogeneity_num,model, true_arms))) +
  geom_hline(aes(yintercept =yint), data=data.frame(yint = seq((0.95-qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)), #llimit,
                                                               (0.95+qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)), # ulimit,
                                                               0.0001)),color = "lightgrey")+# shaded confidence region
  geom_hline(aes(yintercept = 0.95), lty=2, linewidth=0.5)+
  geom_line(linewidth = 1) +
  # geom_vline(xintercept = greylines$x,linetype = "dashed",color = "grey40", position = position_nudge(-0.5)) + #NEW
  geom_point(size = 2) +
  facet_grid(true_model~paramtype+true_arms, scales = "free_y", switch="y") +
  labs(color="Implemented Model")+
  scale_color_manual(values=c("NMA"="#FAA613","CCIE"="#000080",  "ICIE"="#D7263D","IUMIE"="#FF33CC","CUMIE"="#1C77C3"))+
  #ggtitle("Bias Across Scenarios")+
  # ylab(label = "Mean Cov")+
  scale_y_continuous(labels = function(x) formatC(x, format = "f", digits = 2))+
  scale_x_discrete(labels=function(x)case_when(grepl("0.35\\\n0\\n",x)==TRUE~gsub("-","",x),
                                               grepl("0.35",x)==TRUE~gsub("-|twoarm|mixarm","",x),
                                               TRUE~sub("\\n.*","",x)))+
  xlab(label = "")+
  theme_bw() +
  theme(
    strip.background = element_rect(fill = c("#F5F8FC")),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    panel.spacing = unit(0.05, "lines")
  )



# how he anticipates it ##########
coverage_graph2<-graph_res2%>%
  filter(true_arms!="multiarm")%>%
  filter(paramtype!="tau")%>%
  mutate(x = interaction((keep_prop_num),heterogeneity_num,true_arms, sep="\n"))%>%
  mutate(paramtype=factor(paramtype, levels=c("(Adjusted)\nTreatment Effect", "Treatment:Covariate\nInteraction Effect", "tau")))%>%
  ggplot(aes(x = x,y = prop_coverage, color = model, group = interaction(heterogeneity_num,model, true_arms))) +
  geom_hline(aes(yintercept =yint), data=data.frame(yint = seq((0.95-qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)), #llimit,
                                                               (0.95+qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)), # ulimit,
                                                               0.0001)),color = "lightgrey")+# shaded confidence region
  geom_hline(aes(yintercept = 0.95), lty=2, linewidth=0.5)+
  geom_line(linewidth = 1) +
  geom_vline(xintercept = greylines$x,linetype = "dashed",color = "grey40", position = position_nudge(-0.5)) + #NEW
  geom_point(size = 2) +
  facet_grid(true_model~paramtype, scales = "free_y", switch="y") +
  labs(color="Implemented Model")+
  scale_color_manual(values=c("NMA"="#FAA613","CCIE"="#000080",  "ICIE"="#D7263D","IUMIE"="#FF33CC","CUMIE"="#1C77C3"))+
  #ggtitle("Bias Across Scenarios")+
  # ylab(label = "Mean Cov")+
  scale_y_continuous(labels = function(x) formatC(x, format = "f", digits = 2))+
  scale_x_discrete(labels=function(x)case_when(grepl("0.35\\\n0\\n",x)==TRUE~gsub("-","",x),
                                               grepl("0.35",x)==TRUE~gsub("-|twoarm|mixarm","",x),
                                               TRUE~sub("\\n.*","",x)))+
  xlab(label = "")+
  theme_bw() +
  theme(
    strip.background = element_rect(fill = c("#F5F8FC")),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    panel.spacing = unit(0.05, "lines")
  )


# old coverage #########
greylines<-graph_res2 %>%
  mutate(x = interaction(heterogeneity_num, desc(keep_prop_num),true_arms, sep="\n"))%>%
  distinct(x)%>%
  filter(grepl("0\\\n", x) &
           grepl("1\\\n",x) & 
           grepl("multiarm",x )==FALSE)

coveragegraph<-graph_res2%>%
  # filter(true_assumption=="common")%>%
  filter(true_arms!="multiarm")%>%
  filter(paramtype!="tau")%>%
  mutate(prop_coverage=ifelse(prop_coverage<0.90,0.90,prop_coverage))%>%
  mutate(x = interaction(heterogeneity_num, desc(keep_prop_num),true_arms, sep="\n"))%>%
  ggplot()+
  geom_hline(aes(yintercept =yint), data=data.frame(yint = seq((0.95-qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)), #llimit,
                                                               (0.95+qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)), # ulimit,
                                                               0.0001)),color = "lightgrey")+# shaded confidence region
  geom_hline(aes(yintercept = 0.95), lty=2, linewidth=0.5)+
  geom_col(aes(x=x, y=prop_coverage, fill=model, color=model),#color="white",
           position=position_dodge(width=0.7), width = 0.5)+
  scale_y_continuous(trans=scales::trans_new(name="sq",
                                             transform = function(x)x-0.95,
                                             inverse=function(y) y+0.95),
                     breaks=c(0.9,
                              (0.95+qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)),
                              (0.95-qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)),
                              seq(0.95,1,0.05)),
                     labels = scales::number_format(accuracy = 0.01))+
  geom_vline(xintercept = greylines$x,linetype = "dashed",color = "grey40", position = position_nudge(-0.5)) + #NEW
  facet_grid(true_model~paramtype, scales="free_y", switch="y")+
  scale_fill_manual(values=c("NMA"="#FAA613","CCIE"="#000080",  "ICIE"="#D7263D","IUMIE"="#FF33CC","CUMIE"="#1C77C3"))+
  scale_color_manual(values=c("NMA"="#FAA613","CCIE"="#000080",  "ICIE"="#D7263D","IUMIE"="#FF33CC","CUMIE"="#1C77C3"))+
  scale_x_discrete(labels=function(x)case_when(grepl("0\\\n-1",x)==TRUE~gsub("-","",x),
                                               grepl("0\\\n",x)==TRUE~gsub("-|twoarm|mixarm","",x),
                                               TRUE~sub("\\n.*","",x)))+
  theme_bw()+
  labs(fill="Implemented Model", color="Implemented Model")+
  ylab(label = "Coverage Proportion")+
  xlab(label = "")+
  theme(strip.background = element_rect(fill = c("#F5F8FC")),
        plot.title = element_text(hjust = 0.5),
        panel.spacing = unit(0.05, "lines"))
## export coverage ##########


cov_emfsize<-7.3
emf(file=paste(outputdir,"coverage_exwide1.emf", sep="/"), height=cov_emfsize, width=2.5*cov_emfsize, emfPlus = FALSE)
# png(file=paste(outputdir,"coverage_exwide.png", sep="/"), height=cov_emfsize, width=2.5*cov_emfsize, units = "in", res=300)
coveragegraph+
  facet_grid(true_model~paramtype, scales="free_y", switch="y")+
  theme(#strip.text = element_text(size = 11, face = "bold"),
  legend.position = "bottom",
  legend.box = "horizontal",
  legend.key.size = unit(0.3, "lines"),
  # legend.key.width = unit(0.3*2, "lines")#png
  legend.key.width = unit(0.3*5, "lines")#emf
)
dev.off()

emf(file=paste(outputdir,"coverage_exlong.emf", sep="/"), height=cov_emfsize*1.5, width=1.5*cov_emfsize, emfPlus = FALSE)
coveragegraph+
  facet_grid(true_model~paramtype, scales="free_y", switch="y")+
  theme(#strip.text = element_text(size = 11, face = "bold"),
  legend.position = "bottom",
  legend.box = "horizontal",
  legend.key.size = unit(0.3, "lines"),
  # legend.key.width = unit(0.3*2, "lines")#png
  legend.key.width = unit(0.3*5, "lines")#emf
)
dev.off()

# Bias trends in density (heterogeneity is constant) ##########

greylines<-graph_res2 %>%
  mutate(x = interaction((keep_prop_num),heterogeneity_num, true_arms, sep="\n"))%>%
  distinct(x)%>%
  filter(grepl("0.35", x) & grepl("0\\\n",x) & 
           grepl("multiarm",x )==FALSE)

biasgraph<-graph_res2%>%
  filter(true_arms!="multiarm")%>%
  mutate(x = interaction((keep_prop_num),heterogeneity_num, true_arms, sep="\n"))%>%
  mutate(paramtype=factor(paramtype, levels=c("(Adjusted)\nTreatment Effect", "Treatment:Covariate\nInteraction Effect", "tau")))%>%
  ggplot(aes(x = x,y = mean_bias, color = model, group = interaction(heterogeneity_num,model, true_arms))) +
  geom_line(linewidth = 1) +
  geom_hline(aes(yintercept = 0))+
  geom_vline(xintercept = greylines$x,linetype = "dashed",color = "grey40", position = position_nudge(-0.5)) + #NEW
  geom_point(size = 2) +
  facet_grid(true_model~paramtype, scales = "free_y", switch="y") +
  labs(color="Implemented Model")+
  scale_color_manual(values=c("NMA"="#FAA613","CCIE"="#000080",  "ICIE"="#D7263D","IUMIE"="#FF33CC","CUMIE"="#1C77C3"))+
  #ggtitle("Bias Across Scenarios")+
  ylab(label = "Mean Bias")+
  scale_y_continuous(labels = function(x) formatC(x, format = "f", digits = 2))+
  scale_x_discrete(labels=function(x)case_when(grepl("0.35\\\n0\\n",x)==TRUE~gsub("-","",x),
                                               grepl("0.35",x)==TRUE~gsub("-|twoarm|mixarm","",x),
                                               TRUE~sub("\\n.*","",x)))+
  xlab(label = "")+
  theme_bw() +
  theme(
    strip.background = element_rect(fill = c("#F5F8FC")),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    panel.spacing = unit(0.05, "lines")
  )


# MSE trends in density (bias2+var) (heterogeneity is constant) ######
mse_graph<-graph_res2%>%
  filter(true_arms!="multiarm")%>%
  filter(paramtype!="tau")%>%
  mutate(x = interaction((keep_prop_num),heterogeneity_num, true_arms, sep="\n"))%>%
  mutate(paramtype=factor(paramtype, levels=c("(Adjusted)\nTreatment Effect", "Treatment:Covariate\nInteraction Effect", "tau")))%>%
  ggplot(aes(x = x,y = mean_man_mse, color = model, group = interaction(heterogeneity_num,model, true_arms))) +
  geom_line(linewidth = 1) +
  geom_hline(aes(yintercept = 0))+
  geom_vline(xintercept = greylines$x,linetype = "dashed",color = "grey40", position = position_nudge(-0.5)) + #NEW
  geom_point(size = 2) +
  facet_grid(true_model~paramtype, scales = "free_y", switch="y") +
  labs(color="Implemented Model")+
  scale_color_manual(values=c("NMA"="#FAA613","CCIE"="#000080",  "ICIE"="#D7263D","IUMIE"="#FF33CC","CUMIE"="#1C77C3"))+
  #ggtitle("Bias Across Scenarios")+
  ylab(label = "Mean Bias")+
  scale_y_continuous(labels = function(x) formatC(x, format = "f", digits = 2))+
  scale_x_discrete(labels=function(x)case_when(grepl("0.35\\\n0\\n",x)==TRUE~gsub("-","",x),
                                               grepl("0.35",x)==TRUE~gsub("-|twoarm|mixarm","",x),
                                               TRUE~sub("\\n.*","",x)))+
  xlab(label = "")+
  theme_bw() +
  theme(
    strip.background = element_rect(fill = c("#F5F8FC")),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    panel.spacing = unit(0.05, "lines")
  )
## export bias and rmse #############
# size<-700 #if using viewer
emfsize<-8
  # png("temp.png",height=size, width=2*size)
emf(file=paste(outputdir,"mse_ex.emf", sep="/"), height=emfsize*2/3, width=2*emfsize, emfPlus = FALSE)
# png(file=paste(gendir, paste("seed_", seedval,sep=""),"temp.png", sep="/"), height=emfsize, width=emfsize*2, units="in", res=300, pointsize=50)
mse_graph+theme(#strip.text = element_text(size = 11, face = "bold"),
                    legend.position = "bottom",
                    legend.box = "horizontal",
                    legend.key.size = unit(0.3, "lines"),
                    # legend.key.width = unit(0.3*2, "lines")#png
                    legend.key.width = unit(0.3*5, "lines"),#emf
                    axis.text.x=element_text(hjust=0.5),
                    strip.text = element_text(size = 11),
                    text = element_text(size = 12),  
                    axis.text = element_text(size = 8), # need to adjust
                    axis.title = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 11)
                    )
dev.off()

emf(file=paste(outputdir,"bias_ex.emf", sep="/"), height=emfsize, width=2*emfsize, emfPlus = FALSE)
# png(file=paste(gendir, paste("seed_", seedval,sep=""),"temp.png", sep="/"), height=emfsize, width=emfsize*2, units="in", res=300, pointsize=50)
biasgraph+theme(#strip.text = element_text(size = 11, face = "bold"),
  legend.position = "bottom",
  legend.box = "horizontal",
  legend.key.size = unit(0.3, "lines"),
  # legend.key.width = unit(0.3*2, "lines")#png
  legend.key.width = unit(0.3*5, "lines"),#emf
  axis.text.x=element_text(hjust=0.5),
  strip.text = element_text(size = 11),
  text = element_text(size = 12),  
  axis.text = element_text(size = 8), # need to adjust
  axis.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 11)
)
dev.off()
