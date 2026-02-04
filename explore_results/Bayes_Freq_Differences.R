# Differences in Bias Freq vs Bayesian ##############
t1<-graph_res2%>%
  # filter(true_arms=="mixarm" &true_assumption=="NMA" & model=="ICIE" & src=="RMA_MV" & paramtype=="(Adjusted)\nTreatment Effect")%>%
  select(src, paramtype, model, true_model, true_arms,heterogeneity_num,keep_prop_num, starts_with("mean"))%>%
  distinct()%>%
  group_by(src, paramtype, model, true_model, true_arms)%>%#,heterogeneity_num,keep_prop_num
  summarise(pos_bias=sum(mean_bias>0),
            neg_bias=sum(mean_bias<0),
            n=n())%>%
  ungroup()%>%
  mutate(direction=case_when(pos_bias==n ~"positive",
                             neg_bias==n ~"negative",
                             TRUE~"neither"))%>%
  select(true_model, model, paramtype, src, everything())%>%
  arrange(true_model, true_arms, model, paramtype, src)%>%
  tidyr::pivot_wider(values_from = c(pos_bias,neg_bias, n, direction), names_from = src)

# write.csv(t1, "temp.csv")

t1%>%
  filter(true_arms!="multiarm")%>%
  count(direction_JAGS, direction_RMA_MV)%>%
  mutate(direction_JAGS=ifelse(direction_JAGS== direction_RMA_MV, "no change", direction_JAGS),
         direction_RMA_MV=ifelse(direction_JAGS== "no change", "no change", direction_RMA_MV))%>%
  group_by(direction_JAGS, direction_RMA_MV)%>%
  summarise(n=sum(n))%>%
  arrange(desc(n))

# Differences in Coverage Freq vs Bayesian ##############
t2<-graph_res2%>%
  filter(paramtype!="tau")%>%
  # filter(true_arms=="mixarm" &true_assumption=="NMA" & model=="ICIE" & paramtype=="(Adjusted)\nTreatment Effect")%>%
  select(src, paramtype, model, true_model, true_arms,heterogeneity_num,keep_prop_num, prop_coverage)%>%
  distinct()%>%
  mutate(lci=(0.95-qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)),
         uci=(0.95+qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)),
         coverage=case_when(prop_coverage>=uci~"overcoverage",
                            prop_coverage<=lci~"undercoverage",
                            prop_coverage>=lci&prop_coverage<=uci~"nominal"))%>%
  select(true_model, true_arms, model, paramtype, heterogeneity_num, keep_prop_num, src, everything())%>%
  arrange(true_model, true_arms, model, paramtype, heterogeneity_num, keep_prop_num, src)%>%
  tidyr::pivot_wider(values_from = c(prop_coverage,coverage), names_from = src)

# write.csv(t2, "temp1.csv")

t2%>%
  filter(true_arms!="multiarm")%>%
  count(coverage_JAGS, coverage_RMA_MV)%>%
  mutate(coverage_JAGS=ifelse(coverage_JAGS== coverage_RMA_MV, "no change", coverage_JAGS),
         coverage_RMA_MV=ifelse(coverage_JAGS== "no change", "no change", coverage_RMA_MV))%>%
  group_by(coverage_JAGS, coverage_RMA_MV)%>%
  summarise(n=sum(n))
# Differences in MSE Trends ##############
t3<-
  graph_res2%>%
  filter(paramtype!="tau")%>%
  # filter(true_arms=="mixarm" &true_assumption=="NMA" & model=="ICIE" & paramtype=="(Adjusted)\nTreatment Effect")%>%
  select(src, paramtype, model, true_model, true_arms,heterogeneity_num,keep_prop_num, mean_man_mse)%>%
  distinct()%>%
  arrange(paramtype, model, true_model, true_arms, heterogeneity_num, src)%>%
  group_by(paramtype, model, true_model, true_arms, heterogeneity_num, src)%>%#,heterogeneity_num,keep_prop_num
  mutate(mse_dir=mean_man_mse-lag(mean_man_mse))%>%
  summarise(pos_mse=sum(mse_dir>0, na.rm = TRUE),#shiny

            neg_mse=sum(mse_dir<0, na.rm = TRUE),#mse decreasing as density decreases
            n=sum(is.na(mse_dir)==FALSE),
            direction=case_when(pos_mse==n ~"positive",neg_mse==n ~"negative",TRUE~"neither"))%>%
  select(true_model, model, paramtype, src, everything())%>%
  arrange(true_model, true_arms, model, paramtype, src)%>%
  tidyr::pivot_wider(values_from = c(pos_mse,neg_mse, n, direction), names_from = src)

t4<-
  graph_res2%>%
  filter(paramtype!="tau")%>%
  # filter(true_arms=="mixarm" &true_assumption=="NMA" & model=="ICIE" & paramtype=="(Adjusted)\nTreatment Effect")%>%
  select(src, paramtype, model, true_model, true_arms,heterogeneity_num,keep_prop_num, mean_man_mse)%>%
  distinct()%>%
  arrange(paramtype, model, true_model, true_arms, keep_prop_num, src)%>%
  group_by(paramtype, model, true_model, true_arms, keep_prop_num, src)%>%#,heterogeneity_num,keep_prop_num
  mutate(mse_dir=mean_man_mse-lag(mean_man_mse))%>%
  summarise(pos_mse=sum(mse_dir>0, na.rm = TRUE),#mse increasing as heterogeneity decreases
            neg_mse=sum(mse_dir<0, na.rm = TRUE),#mse decreasing as heterogeneity decreases
            n=sum(is.na(mse_dir)==FALSE),
            direction=case_when(pos_mse==n ~"positive",neg_mse==n ~"negative",TRUE~"neither"))%>%
  select(true_model, model, paramtype, src, everything())%>%
  arrange(true_model, true_arms, model, paramtype, src)%>%
  tidyr::pivot_wider(values_from = c(pos_mse,neg_mse, n, direction), names_from = src)

write.csv(t3, "temp1.csv")

t3%>%
  filter(true_arms!="multiarm")%>%
  count(direction_JAGS, direction_RMA_MV)%>%
  mutate(direction_JAGS=ifelse(direction_JAGS== direction_RMA_MV, "no change", direction_JAGS),
         direction_RMA_MV=ifelse(direction_JAGS== "no change", "no change", direction_RMA_MV))%>%
  group_by(direction_JAGS, direction_RMA_MV)%>%
  summarise(n=sum(n))%>%
  arrange(desc(n))

t4%>%
  filter(true_arms!="multiarm")%>%
  count(direction_JAGS, direction_RMA_MV)%>%
  mutate(direction_JAGS=ifelse(direction_JAGS== direction_RMA_MV, "no change", direction_JAGS),
         direction_RMA_MV=ifelse(direction_JAGS== "no change", "no change", direction_RMA_MV))%>%
  group_by(direction_JAGS, direction_RMA_MV)%>%
  summarise(n=sum(n))%>%
  arrange(desc(n))
