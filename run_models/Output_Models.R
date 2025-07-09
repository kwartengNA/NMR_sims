## ---- combres----------------------------------------------------------------------------
if(exists("all_jags_res_df")==FALSE){
  all_jags_res_df<-NULL} else 
  {
    all_jags_res_df<-all_jags_res_df%>%
      mutate(simulation=as.integer(simulation))
  }
if(exists("all_freq_res_df")==FALSE){
  all_freq_res_df<-NULL
} else {
  all_freq_res_df<-all_freq_res_df%>%
    mutate(simulation=as.integer(simulation))
}

all_res<-all_jags_res_df%>%
  bind_rows(all_freq_res_df)%>%
  mutate(bias2=bias^2)
rm(all_jags_res_df)
gc()

