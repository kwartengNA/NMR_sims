# HEADER ##############################
#' Summary method for objects of class netmetareg
#'
#' @description
#' Summary method for objects of class \code{netmetareg}. This is a
#' wrapper function for the outputs of the R function 
#' \code{\link[netmeta]{netmetareg}} in the R
#' package \bold{netmeta} and the R function \code{\link[metafor]{rma.mv}} in the R
#' package \bold{metafor}(Viechtbauer 2010).
#'
#' @details
#' Note, results are not back-transformed in printouts of
#' network meta-analyses using summary measures with transformations, e.g.,
#' log risk ratios are printed instead of the risk ratio if argument
#' \code{sm = "RR"}.
#'
#' @param object An object of class \code{netmetareg}.
#'
#' @return
#' A list of class "summary.netmeta" is returned with the following elements:
#' \item {summary} {Results for all basic comparisons  based
#'   on random effects model (data frame with columns TE, seTE,
#'   rowname, parametertype, treatment).}
#' \item {all_combinations} {Results for all pairwise comparisons based
#'   on random effects model (data frame with columns comp, TE_d_NMR,
#'   seTE_d_NMR, TE_beta_NMR, seTE_beta_NMR, treat1, treat2, cov_NMR).}
#'   In addition, a list \code{.netmetareg} Please view the help page of 
#'   R function \code{\link[netmeta]{netmetaref}} for more details on the 
#'   output from this function.
#' \item {graphdata}
#' \item {graphdata_se_area} {Predictions for all pairwise comparisons based
#'   on random effects model and user specified covariate values 
#'   (data frame with columns comp, xval, has_direct, TE_d_NMR,
#'   seTE_d_NMR, TE_beta_NMR, seTE_beta_NMR, treat1, treat2, cov_NMR,
#'   pred, var_d_beta, se_d_beta, ciymin, ciymax).}
#' @author Nana-Adjoa Kwarteng
#'   \email{nana-adjoa.kwarteng@@uniklinik-freiburg.de},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{netmeta}}
#'
#' @references
#' Viechtbauer W (2010):
#' Conducting Meta-Analyses in R with the Metafor Package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#'
#' @keywords models regression summary
#' 
#' @examples
#' # Add examples
#' 
#' @rdname netmetareg
#' @method netmetareg netmeta
#' @export 
# SUMMARY FUNCTION ###########################################
mini_summary<-function(object){
  
  covar.name<-object$.netmeta$covar.name
  covar_numeric<-is.numeric(object$.netmeta$covar)|is.integer(object$.netmeta$covar)
  
  res_b_df <- cbind.data.frame(object$b, object$se, object$zval, object$pval, object$ci.lb, object$ci.ub)
  names(res_b_df)<-c("est","SE", "zval", "pval", "ci.lb", "ci.ub")
  res_b_df$rowname <- rownames(res_b_df)
  res_b_TE <- res_b_df[grepl(":", res_b_df$rowname)==FALSE, ]
  
  # get TE estimates
  res_b_TE$treat<-res_b_TE$rowname
  res_b_TE$covar<-NA
  res_b_TE$covar_lvl<-NA
  res_b_TE$covar_ref<-NA
  res_b_TE$param_type<-"d"
  
  if(covar.name!="NULL"){
  # get estimates covariate levels for each treatment
  res_b_interactions <- res_b_df[grepl(":", res_b_df$rowname), ]
  # Split 'rowname' using strsplit
  split_names <- do.call(rbind, strsplit(res_b_interactions$rowname, "[:]"))
  res_b_interactions$param_type<-"beta"
  res_b_interactions$treat <- split_names[, 1]
  res_b_interactions$covar <- covar.name#split_names[, 2]
  
  if(covar_numeric==FALSE){
    # get observed covariate levels for each treatment
    dat<-object$.netmeta$x$data[, c("treat1","treat2", covar.name)]
    # Convert to long format
    unique_lvls<-rbind(setNames(dat[,c("treat1",covar.name)],c("treat","lvl")),
        setNames(dat[,c("treat2",covar.name)],c("treat","lvl")))
    unique_lvls <- unique(unique_lvls)
    if(object$.netmeta$assumption=="common"){
      unique_lvls$treat<-ifelse(unique_lvls$treat==object$.netmeta$reference.group,"nonref",unique_lvls$treat)
    }
    # capture estimated factor levels
    res_b_interactions$lvl <- gsub(covar.name,"",split_names[, 2])
    
    # Get the non-included covariate levels aka covariate references
    ref_lvls<-merge(unique_lvls, res_b_interactions, by = c("treat","lvl"), all.x = TRUE)
    ref_lvls<-ref_lvls[is.na(ref_lvls$est)==TRUE,c("treat","lvl")]
    colnames(ref_lvls)[colnames(ref_lvls) == "lvl"] <- "covar_ref"
    res_b_interactions<-merge(res_b_interactions, ref_lvls, by = c("treat"), all.x = TRUE) # should be 1 to many. If there is an issue then it is likely in rma.mv
    colnames(res_b_interactions)[colnames(res_b_interactions) == "lvl"] <- "covar_lvl"
  } else{
    res_b_interactions$covar_lvl=1
    res_b_interactions$covar_ref=0
    if(is.logical(object$.netmeta$x$data[[covar.name]])){
      res_b_interactions$covar_lvl<-as.logical(res_b_interactions$covar_lvl)
      res_b_interactions$covar_ref<-as.logical(res_b_interactions$covar_ref)
    }
  }
  res<-rbind.data.frame(res_b_TE[,c(names(res_b_interactions))],
                        res_b_interactions[,c(names(res_b_interactions))])
  
  } else{
    res<-res_b_TE
  }
  
  res$ref_treat<-object$.netmeta$reference.group
  res<-res[,c("rowname","treat","ref_treat","param_type", "covar","covar_lvl", "covar_ref", "est","SE", "zval", "pval", "ci.lb", "ci.ub")]# order the columns
  rownames(res)<-NULL
  return(res)
}
summary.netmetareg_old<-function(object,
                             min_x=floor(min(object$.netmeta$covar,na.rm=TRUE)),
                             max_x=ceiling(max(object$.netmeta$covar,na.rm=TRUE)),
                             increment_by=0.5,
                             graphset=TRUE,
                             warn.deprecated = gs("warn.deprecated")){
  netmeta:::chkclass(object, "netmetareg")
  # object <- netmeta:::updateversion(object) # what does this do
  # is.bin <- inherits(object, "netmetabin") # what is this for
  
  # extract NMR attributes/assumptions
  assumption=object$.netmeta$assumption
  reference.group=object$.netmeta$reference.group
  consistency=object$.netmeta$consistency
  covar=object$.netmeta$covar
  covar_is_factor<-is.factor(covar)
  
  # mini_summary #######################
  summary_obj<-as_df(object)
  # extract rma mv estimates ########################
  all_consistent_covariance_df<-data.frame(.comp=NA,cov_NMR=NA)# otherwise NMA will not run
  rma_mv_res<-cbind.data.frame(object$b,
                               object$se)
  getvcov<-as.data.frame(vcov(object))
  names(rma_mv_res)<-c("TE","seTE")
  rma_mv_res$rowname<-rownames(rma_mv_res)
  rma_mv_res$parametertype=ifelse(grepl(":",rma_mv_res$rowname),"beta", "d")
  rma_mv_res$treatment<-ifelse(gsub(":.*|`","",rma_mv_res$rowname)=="","all",gsub(":.*|`","",rma_mv_res$rowname))
  rownames(rma_mv_res)<-NULL
  
  # generate matrix combinations for consistent parameters ########
  ## Treatment Effect Consistency ################################
  # if(treatment consistency is TRUE use the consistency equations for the adjusted treatment effects - which it always is for the current function){
    all_consistent_d_TEs<-outer(c(0,rma_mv_res[
      rma_mv_res$parametertype =="d","TE"]),
      c(0,rma_mv_res[rma_mv_res$parametertype =="d","TE"]), `-`)
    dimnames(all_consistent_d_TEs)<-list(c(reference.group,rma_mv_res[
      rma_mv_res$parametertype =="d","treatment"]),
      c(reference.group,rma_mv_res[
        rma_mv_res$parametertype =="d","treatment"]))
    
    #adjusted treatment effects (d TEs
    getvcovTEs<-getvcov[grepl(":",rownames(getvcov))==FALSE,grepl(":",colnames(getvcov))==FALSE]
    getvcovTEs<-cbind.data.frame(setNames(list(rep(0,nrow(getvcovTEs))), reference.group),getvcovTEs)
    getvcovTEs<-rbind(rep(0,ncol(getvcovTEs)),getvcovTEs)
    rownames(getvcovTEs)<-colnames(getvcovTEs)
    getvarTEs<-outer(diag(as.matrix(getvcovTEs)),diag(as.matrix(getvcovTEs)),"+")-(2*as.matrix(getvcovTEs))
    all_consistent_d_seTEs<-sqrt(getvarTEs)
    
    # convert to df
    all_consistent_d_TEs_df<-NULL
    all_consistent_d_seTEs_df<-NULL
    for(i in colnames(all_consistent_d_TEs)){
      temp_d_TEs_df<-as.data.frame(all_consistent_d_TEs)[i]
      names(temp_d_TEs_df)<-"TE_d_NMR"
      temp_d_TEs_df$.comp<-paste(rownames(temp_d_TEs_df),i,sep="_vs_")
      all_consistent_d_TEs_df<-rbind(temp_d_TEs_df,all_consistent_d_TEs_df)
      
      temp_d_seTEs_df<-as.data.frame(all_consistent_d_seTEs)[i]
      names(temp_d_seTEs_df)<-"seTE_d_NMR"
      temp_d_seTEs_df$.comp<-paste(rownames(temp_d_seTEs_df),i,sep="_vs_")
      all_consistent_d_seTEs_df<-rbind(temp_d_seTEs_df,all_consistent_d_seTEs_df)
    }
  # end if(treatment consistency is TRUE, which it always is for the current function)
  ## Interaction Effect Consistency ################################
  # if(interaction consistency is TRUE use the consistency equations for the interactions - which it always is for the current function)
    if(is.null(covar)==TRUE|covar_is_factor==TRUE){
      # we don't have a solution for categorical covariates
      all_consistent_beta_TEs_df<-all_consistent_beta_seTEs_df<-data.frame(.comp=NA)
    } else if(consistency==TRUE){
    if(assumption=="common"){ # the common assumption has 1 beta per treatment
      all_consistent_beta_TEs<-outer(c(0,rep(rma_mv_res[rma_mv_res$parametertype =="beta","TE"],
                                             length(rma_mv_res[rma_mv_res$parametertype =="d","treatment"]))),
                                     c(0,rep(rma_mv_res[rma_mv_res$parametertype =="beta","TE"],
                                             length(rma_mv_res[rma_mv_res$parametertype =="d","treatment"]))), `-`)
      dimnames(all_consistent_beta_TEs)<-list(c(reference.group,rma_mv_res[rma_mv_res$parametertype =="d","treatment"]),
                                              c(reference.group,rma_mv_res[rma_mv_res$parametertype =="d","treatment"]))
      
      getvcovseTEs<-getvcov[grepl(":",rownames(getvcov))==TRUE,grepl(":",colnames(getvcov))==TRUE]
      all_consistent_beta_seTEs<-outer(c(0,rep(getvcovseTEs,
                                               length(rma_mv_res[rma_mv_res$parametertype =="d","treatment"]))),
                                       c(0,rep(getvcovseTEs,
                                               length(rma_mv_res[rma_mv_res$parametertype =="d","treatment"]))), FUN = function(x,y){abs(x-y)})#variance is always positive anyway
      dimnames(all_consistent_beta_seTEs)<-list(c(reference.group,rma_mv_res[rma_mv_res$parametertype =="d","treatment"]),
                                                c(reference.group,rma_mv_res[rma_mv_res$parametertype =="d","treatment"]))
      all_consistent_beta_seTEs<-sqrt(all_consistent_beta_seTEs)
      
      # covariances ###############
      temp_vcov<-getvcov[grepl(":",colnames(getvcov)),grepl(":",rownames(getvcov))==FALSE]
      all_consistent_covariance_df<-data.frame(cov_NMR=as.numeric(temp_vcov),
                                               .comp=paste(colnames(temp_vcov),reference.group,sep="_vs_"))
      all_consistent_covariance_df<-rbind(all_consistent_covariance_df,
                                          data.frame(cov_NMR=as.numeric(temp_vcov),
                                                     .comp=paste(reference.group,colnames(temp_vcov),sep="_vs_")))# reverse direction
      missing_comps<-all_consistent_d_TEs_df$.comp[!all_consistent_d_TEs_df$.comp %in% all_consistent_covariance_df$.comp]
      all_consistent_covariance_df<-rbind.data.frame(all_consistent_covariance_df,data.frame(.comp=missing_comps,
                                                                                             cov_NMR=0))
      
    } else if(assumption=="independent"){ # independent assumption
      all_consistent_beta_TEs<-outer(c(0,rma_mv_res[rma_mv_res$parametertype =="beta","TE"]),
                                     c(0,rma_mv_res[rma_mv_res$parametertype =="beta","TE"]), `-`)
      dimnames(all_consistent_beta_TEs)<-list(c(reference.group,rma_mv_res[rma_mv_res$parametertype =="beta","treatment"]),
                                              c(reference.group,rma_mv_res[rma_mv_res$parametertype =="beta","treatment"]))
      
      getvcovBetas<-as.matrix(getvcov[grepl(":",rownames(getvcov))==TRUE,grepl(":",colnames(getvcov))==TRUE])
      # getvcovBetas<-cbind(reference.group=rep(0,nrow(getvcovBetas)),getvcovBetas)
      getvcovBetas<-cbind.data.frame(setNames(list(rep(0,nrow(getvcovBetas))), reference.group),getvcovBetas)
      getvcovBetas<-rbind(rep(0,ncol(getvcovBetas)),getvcovBetas)
      rownames(getvcovBetas)<-colnames(getvcovBetas)
      getvarBetas<-outer(diag(as.matrix(getvcovBetas)),diag(as.matrix(getvcovBetas)),"+")-(2*as.matrix(getvcovBetas))
      all_consistent_beta_seTEs<-sqrt(getvarBetas)
      
      dimnames(all_consistent_beta_seTEs)<-list(c(reference.group,rma_mv_res[rma_mv_res$parametertype =="beta","treatment"]),
                                                c(reference.group,rma_mv_res[rma_mv_res$parametertype =="beta","treatment"]))
      
      # covariances ###############
      # make the interaction vs treatment matrix a proper square matrix, interactions in rows, treats in cols
      temp_vcov<-getvcov[grepl(":",colnames(getvcov)),grepl(":",rownames(getvcov))==FALSE]
      rownames(temp_vcov)<-gsub(":.*|`","",rownames(temp_vcov))
      temp_vcov_interactions<-rownames(temp_vcov)
      temp_vcov[names(temp_vcov)[!names(temp_vcov)%in%temp_vcov_interactions],]<-0
      temp_vcov[reference.group]<-0
      temp_vcov[reference.group,]<-0
      temp_vcov<-temp_vcov[order(rownames(temp_vcov)),order(colnames(temp_vcov))]
      
      # all_consistent_covariance_df
      temp_vcov2<-temp_vcov
      temp_vcov2<-as.data.frame(temp_vcov2)
      temp_vcov2<-reshape(temp_vcov2,
                          idvar = "interaction",
                          ids=rownames(temp_vcov2),
                          direction = "long", 
                          varying = list(names(temp_vcov2)),
                          v.names = "value", 
                          timevar = "treat", 
                          times = names(temp_vcov2))
      temp_vcov2$type<-ifelse(temp_vcov2$treat==temp_vcov2$interaction,"ref","interim")
      
      interim<-temp_vcov2[,!names(temp_vcov2)%in%c("type")]
      names(interim)<-gsub("value","interim",names(interim))
      
      ref<-temp_vcov2[temp_vcov2$type=="ref",!names(temp_vcov2)%in%c("treat","type")]
      names(ref)<-gsub("value","ref",names(ref))
      
      comb1<-merge(ref,interim, by="interaction")
      comb1<-merge(comb1,interim, by.x=c("treat","interaction"),by.y=c("interaction","treat"))
      comb1<-merge(comb1,ref, by.x="treat",by.y="interaction")
      comb1$cov_NMR<-ifelse(comb1$interaction==comb1$treat,comb1$ref.x,comb1$ref.x+comb1$ref.y-comb1$interim.x-comb1$interim.y)
      comb1$.comp<-ifelse(comb1$interaction==comb1$treat,
                         paste(comb1$interaction,reference.group,sep="_vs_"),
                         paste(comb1$interaction,comb1$treat,sep="_vs_"))
      all_consistent_covariance_df<-unique(comb1[,c(".comp","cov_NMR")])
    }
    
    all_consistent_beta_TEs_df<-NULL
    all_consistent_beta_seTEs_df<-NULL
    for(i in colnames(all_consistent_beta_TEs)){
      temp_beta_TEs_df<-as.data.frame(all_consistent_beta_TEs)[i]
      names(temp_beta_TEs_df)<-"TE_beta_NMR"
      temp_beta_TEs_df$.comp<-paste(rownames(temp_beta_TEs_df),i,sep="_vs_")
      all_consistent_beta_TEs_df<-rbind(temp_beta_TEs_df,all_consistent_beta_TEs_df)
      
      temp_beta_seTEs_df<-as.data.frame(all_consistent_beta_seTEs)[i]
      names(temp_beta_seTEs_df)<-"seTE_beta_NMR"
      temp_beta_seTEs_df$.comp<-paste(rownames(temp_beta_seTEs_df),i,sep="_vs_")
      all_consistent_beta_seTEs_df<-rbind(temp_beta_seTEs_df,all_consistent_beta_seTEs_df)
    }
    
    
    
    
  }
  # end if(interaction consistency is TRUE use the consistency equations for the interactions - which it always is for the current function)
  ## Interaction UMIE ################################
  # if(interaction consistency is FALSE only use observed values for the interactions - which is not currently available for the current function) # NOT TESTED
  if(consistency==FALSE){
    if(assumption=="common"){# common
      all_consistent_beta_TEs<-matrix(nrow = length(rma_mv_res[rma_mv_res$parametertype =="d","treatment"])+1,
                                      ncol=length(rma_mv_res[rma_mv_res$parametertype =="d","treatment"])+1,
                                      data=rma_mv_res[rma_mv_res$parametertype =="beta","TE"])
      diag(all_consistent_beta_TEs)<-0
      dimnames(all_consistent_beta_TEs)<-list(c(reference.group,rma_mv_res[rma_mv_res$parametertype =="d","treatment"]),
                                              c(reference.group,rma_mv_res[rma_mv_res$parametertype =="d","treatment"]))
      
      all_consistent_beta_seTEs<-matrix(nrow = length(rma_mv_res[rma_mv_res$parametertype =="d","treatment"])+1,
                                        ncol=length(rma_mv_res[rma_mv_res$parametertype =="d","treatment"])+1,
                                        data=rma_mv_res[rma_mv_res$parametertype =="beta","seTE"])
      diag(all_consistent_beta_seTEs)<-0
      dimnames(all_consistent_beta_seTEs)<-list(c(reference.group,rma_mv_res[rma_mv_res$parametertype =="d","treatment"]),
                                                c(reference.group,rma_mv_res[rma_mv_res$parametertype =="d","treatment"]))
      
      all_consistent_beta_TEs_df<-NULL
      all_consistent_beta_seTEs_df<-NULL
      for(i in colnames(all_consistent_beta_TEs)){
        temp_beta_TEs_df<-as.data.frame(all_consistent_beta_TEs)[i]
        names(temp_beta_TEs_df)<-"TE_beta_NMR"
        temp_beta_TEs_df$.comp<-paste(rownames(temp_beta_TEs_df),i,sep="_vs_")
        all_consistent_beta_TEs_df<-rbind(temp_beta_TEs_df,all_consistent_beta_TEs_df)
        
        temp_beta_seTEs_df<-as.data.frame(all_consistent_beta_seTEs)[i]
        names(temp_beta_seTEs_df)<-"seTE_beta_NMR"
        temp_beta_seTEs_df$.comp<-paste(rownames(temp_beta_seTEs_df),i,sep="_vs_")
        all_consistent_beta_seTEs_df<-rbind(temp_beta_seTEs_df,all_consistent_beta_seTEs_df)
      }
      
      # covariances ###############
      temp_vcov<-getvcov[grepl(":",colnames(getvcov)),grepl(":",rownames(getvcov))==FALSE]
      rownames(temp_vcov)<-gsub(":.*|`","",rownames(temp_vcov))
      temp_vcov[reference.group]<-0
      temp_vcov2<-as.numeric(temp_vcov)
      names(temp_vcov2)<-names(temp_vcov)
      temp_vcov2<-outer(temp_vcov2,temp_vcov2,"-")
      temp_vcov2<-temp_vcov2[order(rownames(temp_vcov2)),order(colnames(temp_vcov2))]
      temp_vcov2<-as.data.frame(temp_vcov2)
      comb1<-reshape(temp_vcov2,
                     idvar = "interaction",
                     ids=rownames(temp_vcov2),
                     direction = "long", 
                     varying = list(names(temp_vcov2)),
                     v.names = "cov_NMR", 
                     timevar = "treat", 
                     times = names(temp_vcov2))
      comb1$.comp<-paste(comb1$interaction,comb1$treat,sep="_vs_")# the correct order depends on which combination was used in model generation
      all_consistent_covariance_df<-unique(comb1[,c(".comp","cov_NMR")])
      
    }else if(assumption=="independent"){
      all_consistent_beta_TEs_df<-rma_mv_res[rma_mv_res$parametertype =="beta",c("TE","treatment")]
      all_consistent_beta_seTEs_df<-rma_mv_res[rma_mv_res$parametertype =="beta",c("seTE","treatment")]
      #reverse direction in case reverse direction occurs in original data
      all_consistent_beta_TEs_df<-rbind.data.frame(all_consistent_beta_TEs_df,
                                                   data.frame("TE"=all_consistent_beta_TEs_df$TE,#*-1,# reversal depends on interpretation of ivals
                                                              "treatment"=paste(gsub(".*_vs_","",all_consistent_beta_TEs_df$treatment), # without the vs this can get messy
                                                                                gsub("_vs_.*","",all_consistent_beta_TEs_df$treatment),
                                                                                sep=" ")))
      all_consistent_beta_TEs_df$treatment<-gsub(".comp_", "", all_consistent_beta_TEs_df$treatment) # since i added comp
      all_consistent_beta_TEs_df$treatment<-gsub(" ", "_vs_", all_consistent_beta_TEs_df$treatment) 
      all_consistent_beta_seTEs_df<-rbind.data.frame(all_consistent_beta_seTEs_df,
                                                     data.frame("seTE"=all_consistent_beta_seTEs_df$seTE,
                                                                "treatment"=paste(gsub(".*_vs_","",all_consistent_beta_seTEs_df$treatment),
                                                                                  gsub("_vs_.*","",all_consistent_beta_seTEs_df$treatment),sep="_vs_")))
      # all_consistent_beta_seTEs_df<-rbind.data.frame(all_consistent_beta_seTEs_df,
      #                                              data.frame("seTE"=all_consistent_beta_seTEs_df$seTE,#*-1,# reversal depends on interpretation of ivals
      #                                                         "treatment"=paste(gsub(".* ","",all_consistent_beta_seTEs_df$treatment),#gsub(".*_vs_","",all_consistent_beta_seTEs_df$treatment), # without the vs this can get messy
      #                                                                           gsub(" .*","",all_consistent_beta_seTEs_df$treatment), #gsub("_vs_.*","",all_consistent_beta_seTEs_df$treatment),
      #                                                                           sep=" ")))
      all_consistent_beta_seTEs_df$treatment<-gsub(".comp_", "", all_consistent_beta_seTEs_df$treatment) # since i added comp
      # all_consistent_beta_seTEs_df$treatment<-gsub(" ", "_vs_", all_consistent_beta_seTEs_df$treatment) 
      names(all_consistent_beta_TEs_df)<-c("TE_beta_NMR",".comp")
      names(all_consistent_beta_seTEs_df)<-c("seTE_beta_NMR",".comp")
      
      # covariances ###############
      # make the interaction vs treatment matrix a proper square matrix, interactions in rows, treats in cols
      temp_vcov<-getvcov[grepl(":",colnames(getvcov)),grepl(":",rownames(getvcov))==FALSE]
      rownames(temp_vcov)<-gsub("\\.comp_|:.*|`","",rownames(temp_vcov))
      temp_vcov_interactions<-rownames(temp_vcov)
      temp_vcov[names(temp_vcov)[!names(temp_vcov)%in%temp_vcov_interactions],]<-0
      temp_vcov[reference.group]<-0
      temp_vcov[reference.group,]<-0
      temp_vcov<-temp_vcov[order(rownames(temp_vcov)),order(colnames(temp_vcov))]
      temp_vcov2<-temp_vcov
      temp_vcov2<-as.data.frame(temp_vcov2)
      temp_vcov2<-reshape(temp_vcov2,
                          idvar = "interaction",
                          ids=rownames(temp_vcov2),
                          direction = "long", 
                          varying = list(names(temp_vcov2)),
                          v.names = "value", 
                          timevar = "treat", 
                          times = names(temp_vcov2))
      temp_vcov2$int_treat1<-gsub(".*_vs_","",temp_vcov2$interaction)
      temp_vcov2$int_treat2<-gsub("_vs_.*","",temp_vcov2$interaction)
      
      temp_vcov2$type<-ifelse(temp_vcov2$treat==temp_vcov2$int_treat1&temp_vcov2$treat==temp_vcov2$int_treat2,"ref","interim")
      
      interim<-temp_vcov2[temp_vcov2$type=="interim"&
                            (temp_vcov2$treat==temp_vcov2$int_treat1|
                               temp_vcov2$treat==temp_vcov2$int_treat2),
                          !names(temp_vcov2)%in%c("type")]
      names(interim)<-gsub("value","interim",names(interim))
      
      ref<-temp_vcov2[temp_vcov2$type=="ref",!names(temp_vcov2)%in%c("treat","type","int_treat1", "int_treat2")]
      names(ref)<-gsub("value","ref",names(ref))
      
      comb1<-merge(interim,ref, by.x="int_treat2",by.y="interaction")
      comb1<-merge(comb1,ref, by.x="int_treat1",by.y="interaction")
      comb1<-merge(comb1,interim, by.x=c("interaction","int_treat1","int_treat2"),by.y=c("interaction","int_treat1","int_treat2"))
      comb1<-comb1[comb1$treat.x==comb1$int_treat2,]# only in the same order as the covariate
      
      comb1$cov_NMR<-comb1$ref.x+comb1$ref.y+comb1$interim.x-comb1$interim.y
      comb1$.comp<-paste(comb1$treat.x,comb1$treat.y,sep="_vs_")
      combrev<-comb1
      combrev$.comp<-paste(combrev$treat.y,combrev$treat.x,sep="_vs_")#reverse, does not change covariance
      comb1<-rbind(comb1,combrev)
      all_consistent_covariance_df<-unique(comb1[,c(".comp","cov_NMR")])
    }
  }
  # end if(interaction consistency is FALSE)
  # combine elements for summary and all combinations #########################
  all_combinations<-merge(all_consistent_d_TEs_df,all_consistent_d_seTEs_df,by=".comp",all.x=TRUE)
  
  all_combinations<-merge(all_combinations,all_consistent_beta_TEs_df,by=".comp",all.x=TRUE)
  all_combinations<-merge(all_combinations,all_consistent_beta_seTEs_df,by=".comp",all.x=TRUE)
  
  all_combinations$treat1<-substr(all_combinations$.comp,1,regexpr("_",all_combinations$.comp)-1)
  all_combinations$treat2<-substr(all_combinations$.comp,regexpr("vs_",all_combinations$.comp)+3,nchar(all_combinations$.comp))# treat2 is the reference treatment or second term. TE=treat1-treat2
  
  all_combinations<-merge(all_combinations,all_consistent_covariance_df,by=".comp",all.x=TRUE)#covariances
  all_combinations$cov_NMR[all_combinations$treat1==all_combinations$treat2]<-0 # comparisons with themselves as reference have a cov of 0
  # create plot data #####################
  if(graphset==TRUE){ # FALSE would save disk space
    ## scatter plot ##########
    pair_data<-object$.netmeta$x$data# original data
    pair_data_rev<-pair_data
    # do the reversal
    nameroots<-names(pair_data_rev)[grepl("1$|2$",names(pair_data_rev))]
    nameroots<-table(gsub("1$|2$","",nameroots))
    rev_variables<-names(nameroots[nameroots==2])
    for(i in rev_variables) {
      assign(paste(i,"1temp",sep=""),pair_data_rev[[paste(i,1,sep="")]])
      pair_data_rev[[paste(i,1,sep="")]]<-pair_data_rev[[paste(i,2,sep="")]]
      pair_data_rev[[paste(i,2,sep="")]]<-get(paste(i,"1temp",sep=""))
    }
    pair_data_rev$TE<-pair_data_rev$TE*-1
    pair_data_rev$.pairwise_src<-"reverse"
    pair_data$.pairwise_src<-"original"
    # combine 
    pair_data<-rbind(pair_data,pair_data_rev)
    pair_data$.comp<-paste(pair_data$treat1,pair_data$treat2,sep="_vs_")
    
    # add summary results
    # for the NMA results, users should really specify a minx and maxx because there is no intuitive region. but if they don't, then 0,1
    if(is.null(all_combinations$TE_beta_NMR)==TRUE){
      all_combinations$TE_beta_NMR<-NA
      all_combinations$seTE_beta_NMR<-NA
      min_x=ifelse(abs(min_x)==Inf,0, min_x)
      max_x=ifelse(abs(max_x)==Inf,1, max_x)
    }
    graphdata<-merge(all_combinations,pair_data,
                     by.x=c(".comp","treat1","treat2"),
                     by.y=c(".comp","treat1","treat2"),all.x=TRUE)
    
    graphdata$TE_beta_NMR<-ifelse(is.na(graphdata$TE_beta_NMR),0,graphdata$TE_beta_NMR)
    graphdata$treat1<-ifelse(is.na(graphdata$treat1),gsub("_vs.*","",graphdata$.comp),graphdata$treat1)# contrast based might lead to missing treat1/2
    graphdata$treat2<-ifelse(is.na(graphdata$treat2),gsub(".*_vs_","",graphdata$.comp),graphdata$treat2)
    treat_order<-unique(c(graphdata$treat1,graphdata$treat2))# users can reorder outside the function
    graphdata$treat1<-factor(graphdata$treat1,levels=treat_order)
    graphdata$treat2<-factor(graphdata$treat2,levels=treat_order)
    
    graphdata$has_direct<-ifelse(is.na(graphdata$studlab)==FALSE,TRUE,FALSE)
    graphdata<-unique(graphdata)
    
    ## confidence region ################################
    graphdata_se_area<-data.frame(.comp=unique(graphdata$.comp),
                                  xval=rep(seq(min_x,max_x,by=increment_by),
                                           each=length(unique(graphdata$.comp))))
    graphdata_se_area<-merge(graphdata_se_area,unique(graphdata[,c(".comp","has_direct","TE_d_NMR","TE_beta_NMR","seTE_d_NMR","seTE_beta_NMR","cov_NMR","treat1","treat2")]),by=".comp",all.x=TRUE)
    graphdata_se_area$pred<-graphdata_se_area$TE_d_NMR+(graphdata_se_area$TE_beta_NMR*graphdata_se_area$xval)
    graphdata_se_area$var_d_beta<-(graphdata_se_area$seTE_d_NMR^2)+((graphdata_se_area$seTE_beta_NMR^2)*(graphdata_se_area$xval^2))+(2*graphdata_se_area$xval*ifelse(is.na(graphdata_se_area$cov_NMR),0,graphdata_se_area$cov_NMR))
    
    graphdata_se_area$se_d_beta<-sqrt(graphdata_se_area$var_d_beta)
    
    graphdata_se_area$ciymin=(graphdata_se_area$pred-(qnorm(1-(0.05/2))*graphdata_se_area$se_d_beta))
    graphdata_se_area$ciymax=(graphdata_se_area$pred+(qnorm(1-(0.05/2))*graphdata_se_area$se_d_beta))
    
  } else {
    graphdata<-FALSE
    graphdata_se_area<-FALSE}
  # return ######################
  res<-NULL
  # res$.netmeta<-object$.netmeta # alternative
  res$summary<-summary_obj
  res$all_combinations<-all_combinations
  res$graphdata<-graphdata
  res$graphdata_se_area<-graphdata_se_area
  res$.netmetareg<-object
  return(res)
}
# GRAPH FUNCTIONS ####################
light_plot.netmetareg<-function(graphset=NULL, graphdata=NULL, graphdata_se_area=NULL, xcov="xcov", confidence_color=NULL, diamond_color=NULL,se=TRUE,interaction_color=NULL , study_point_color=NULL, study_point_shape=NULL, diamond_width=0.2, 
                                plot_type="flashlight", unique_comparisons=FALSE){
  if(is.null(graphset) & (is.null(graphdata)|is.null(graphdata_se_area))){
    print("some error about missing reqs")
  } else if((is.null(graphdata)|is.null(graphdata_se_area))){
    graphdata<-graphset$graphdata
    graphdata_se_area<-graphset$graphdata_se_area
    xcov<-graphset$.netmetareg$.netmeta$covar.name
  }
  
  if(unique_comparisons==TRUE){
    graphdata<-graphdata[as.numeric(graphdata$treat2)<as.numeric(graphdata$treat1),]
    graphdata_se_area<-graphdata_se_area[as.numeric(graphdata_se_area$treat2)<as.numeric(graphdata_se_area$treat1),]
  }
  # Adjusted TEs as intercepts # make a diamond instead, adj_TE_color
  temp_graph<- graphdata%>% # optional filters
    mutate(lTE_d_NMR=TE_d_NMR - (qnorm(1-(0.05/2))*seTE_d_NMR),
           uTE_d_NMR=TE_d_NMR + (qnorm(1-(0.05/2))*seTE_d_NMR),
           TE_d_NMR1=TE_d_NMR,
           TE_d_NMR2=TE_d_NMR)%>%
    pivot_longer(cols=c(TE_d_NMR1,lTE_d_NMR,TE_d_NMR2,uTE_d_NMR),
                 values_to = "Treatment_Effect", names_to = "parametertype")%>%# named for default axis
    mutate(verticesx=case_when(parametertype=="TE_d_NMR1"~-1*diamond_width,
                               parametertype=="TE_d_NMR2"~diamond_width,
                               TRUE~0))%>%
    mutate(parametertype=factor(parametertype,ordered = TRUE,
                                levels=c("TE_d_NMR1","lTE_d_NMR",
                                         "TE_d_NMR2","uTE_d_NMR")))%>%
    arrange(parametertype)%>% # optional filters
    ggplot()
  
  # add confidence region
  if(se==TRUE){
    if(is.null(substitute(confidence_color))==TRUE|is.character(substitute(confidence_color))==TRUE){
      
      temp_graph<-temp_graph+
        geom_ribbon(data=graphdata_se_area,# optional filters
                    aes(x=xval,ymin=ciymin,ymax=ciymax, group=.comp), fill=ifelse(is.null(confidence_color),"black",confidence_color),alpha=0.2)
      
    } else if(is.name(substitute(confidence_color))) {
      temp_graph<-temp_graph+
        geom_ribbon(data=graphdata_se_area,# optional filters
                    aes(x=xval,ymin=ciymin,ymax=ciymax, group=.data[[deparse(substitute(confidence_color))]], fill=.data[[deparse(substitute(confidence_color))]]),alpha=0.2)
    }
  }
  # add diamonds
  if(is.null(substitute(diamond_color))==TRUE|is.character(substitute(diamond_color))==TRUE){
    
    temp_graph<-temp_graph+
      geom_polygon(aes(x=verticesx,y=Treatment_Effect, group=.comp), fill=ifelse(is.null(diamond_color),"black",diamond_color), alpha=0.3)
    
  } else if(is.name(substitute(diamond_color))) {
    temp_graph<-temp_graph+
      geom_polygon(aes(x=verticesx,y=Treatment_Effect, group=.comp, fill=.data[[deparse(substitute(diamond_color))]]), color="black", alpha=0.3)
  }
  
  # add study level observations
  if(plot_type%in%c("flashlight","firefly")){
    if(is.null(substitute(study_point_color))==TRUE|is.character(substitute(study_point_color))==TRUE){
      # temp_graph<-temp_graph+
      #   geom_pointrange(aes(x=xcov,y=TE,ymin=TE-(1.96*seTE),ymax=TE+(1.96*seTE)), color=ifelse(is.null(study_point_color),"black",study_point_color))
      # insert
      if(is.null(substitute(study_point_shape))==TRUE|is.numeric(substitute(study_point_shape))==TRUE){
        temp_graph<-temp_graph+
          geom_pointrange(aes(x=!!sym(xcov),y=TE,ymin=TE-(1.96*seTE),ymax=TE+(1.96*seTE)), color=ifelse(is.null(study_point_color),"black",study_point_color), shape=ifelse(is.null(study_point_shape),16,study_point_shape))
      } else if(is.name(substitute(study_point_shape))) {
        
        temp_graph<-temp_graph+
          geom_pointrange(aes(x=!!sym(xcov),y=TE,ymin=TE-(1.96*seTE),ymax=TE+(1.96*seTE), shape=.data[[deparse(substitute(study_point_shape))]]), color=ifelse(is.null(study_point_color),"black",study_point_color))
      }
      # end insert
      
    } else if(is.name(substitute(study_point_color))) {
      # temp_graph<-temp_graph+
      #   geom_pointrange(aes(x=!!sym(xcov),y=TE,ymin=TE-(1.96*seTE),ymax=TE+(1.96*seTE), color=.data[[deparse(substitute(study_point_color))]]))
      
      # insert
      if(is.null(substitute(study_point_shape))==TRUE|is.numeric(substitute(study_point_shape))==TRUE){
        temp_graph<-temp_graph+
          geom_pointrange(aes(x=!!sym(xcov),y=TE,ymin=TE-(1.96*seTE),ymax=TE+(1.96*seTE), color=.data[[deparse(substitute(study_point_color))]]), shape=ifelse(is.null(study_point_shape),16,study_point_shape))
      } else if(is.name(substitute(study_point_shape))) {
        
        temp_graph<-temp_graph+
          geom_pointrange(aes(x=!!sym(xcov),y=TE,ymin=TE-(1.96*seTE),ymax=TE+(1.96*seTE), shape=.data[[deparse(substitute(study_point_shape))]], color=.data[[deparse(substitute(study_point_color))]]))
      }
      # end insert
      
    }
  }
  
  # add interactions
  if(is.null(substitute(interaction_color))==TRUE|is.character(substitute(interaction_color))==TRUE){
    temp_graph<-temp_graph+
      geom_abline(aes(slope=TE_beta_NMR,intercept=TE_d_NMR), color=ifelse(is.null(interaction_color),"black",interaction_color),lty=2) # Interactions as slopes
  } else if(is.name(substitute(interaction_color))) {
    temp_graph<-temp_graph+
      geom_abline(aes(slope=TE_beta_NMR,intercept=TE_d_NMR,color=.data[[deparse(substitute(interaction_color))]]),lty=2) # Interactions as slopes
  }
  
  if(plot_type=="flashlight"){
    temp_graph<-temp_graph+
      facet_grid(treat1~treat2)+# can be overwritten outside of function
      scale_y_continuous(name="Treatment Effect", sec.axis = sec_axis(~ . , name = "Comparator", breaks = NULL, labels = NULL))+
      scale_x_continuous(name="Covariate Value", sec.axis = sec_axis(~ . , name = "Reference", breaks = NULL, labels = NULL))
  } else if(plot_type %in% c("searchlight","firefly")){
    temp_graph<-temp_graph+
      facet_grid(.~treat2)+# can be overwritten outside of function
      scale_x_continuous(name="Covariate Value", sec.axis = sec_axis(~ . , name = "Reference", breaks = NULL, labels = NULL))+
      scale_y_continuous(name="Treatment Effect")
  } 
  
  
  return(temp_graph)
}
forest.nmr.netmetareg<-function(graphset=NULL, combinations=NULL){
  if(is.null(combinations)){
    combinations<-graphset$all_combinations
  }
  
  temp_graph<-combinations%>%
    select(comp,treat1,treat2,cov_NMR,Est=TE_d_NMR,seEst=seTE_d_NMR)%>%
    mutate(parametertype="d")%>%
    bind_rows(combinations%>%
                select(comp,treat1,treat2,cov_NMR,Est=TE_beta_NMR,seEst=seTE_beta_NMR)%>%
                mutate(parametertype="beta"))%>%
    ggplot()+
    geom_vline(xintercept = 0, color="grey")+
    geom_pointrange(aes(x = Est ,y = treat1, 
                        xmin = Est -(qnorm(1-(0.05/2))*seEst),
                        xmax = Est +(qnorm(1-(0.05/2))*seEst)))+
    facet_grid(treat2~parametertype)+# can be overwritten outside of function
    scale_x_continuous(name="Estimate")
  
  return(temp_graph)
}

multimodel_light_plot.netmetareg<-function(graphsets, reference=NULL, comparator=NULL,plot_type="fairylight",xcov="xcov",
                                has_direct=FALSE,unique_comparisons=TRUE, # or effect =  interaction
                                effect = "treatment"){
  if(is.null(comparator)==FALSE & is.null(reference)==TRUE){
    reference<-comparator
    comparator<-NULL
    print(paste(reference, "used as reference"))
  }
  if(is.null(names(graphsets))| sum(names(graphsets)=="")>0){
    print("At least one graphset is unnamed. Please provide a named list")
  } else{
    # combine list data
    for(i in 1:length(graphsets)){
      tempgraphdata<-graphsets[[i]]$graphdata
      tempgraphdata$.model<-names(graphsets)[[i]]
      tempsedata<-graphsets[[i]]$graphdata_se_area
      tempsedata$.model<-names(graphsets)[[i]]
      if(i==1){
        allgraphdata<-tempgraphdata
        allsedata<-tempsedata
      } else {
        allgraphdata<-bind_rows(tempgraphdata,allgraphdata)# yes it's dplyr
        allsedata<-bind_rows(tempsedata,allsedata) # yes it's dplyr
      }
    }
    
    # general transformations 
    allgraphdata$xcov<-allgraphdata[[xcov]]
    allgraphdata$treat2<-as.factor(allgraphdata$treat2)
    allgraphdata$treat1<-as.factor(allgraphdata$treat1)
    if(is.null(reference)==FALSE){
      allgraphdata$treat1<-relevel(allgraphdata$treat1,ref = reference)
      allgraphdata$treat2<-relevel(allgraphdata$treat2,ref = reference)
      allsedata$treat1<-relevel(as.factor(allsedata$treat1),ref = reference)
      allsedata$treat2<-relevel(as.factor(allsedata$treat2),ref = reference)
    }
    
    allgraphdata$.comp<-gsub("_vs_","\n",allgraphdata$.comp)
    allsedata$.comp<-gsub("_vs_","\n",allsedata$.comp)
    allsedata$.model_comp<-paste(allsedata$.model, allsedata$.comp, sep="\n")
    allgraphdata$.model_comp<-paste(allgraphdata$.model, allgraphdata$.comp, sep="\n")
    
    
    # filter for comparisons containing direct evidence if specified
    if(has_direct==TRUE){
      allgraphdata<-allgraphdata[allgraphdata$has_direct==TRUE,]
      allsedata<-allsedata[allsedata$has_direct==TRUE,]
    }
    # filter if a specific reference and comparator is specified
    if(is.null(reference)==FALSE){
      allgraphdata<-allgraphdata[allgraphdata$treat2%in%reference,]
      allsedata<-allsedata[allsedata$treat2%in%reference,]
      
      # only filter for comparator if reference is filtered
      if(is.null(comparator)==FALSE){
        allgraphdata<-allgraphdata[allgraphdata$treat1%in%comparator,]
        allsedata<-allsedata[allsedata$treat1%in%comparator,]
      }
      
    } 
    
    if(unique_comparisons==TRUE){
      allgraphdata<-allgraphdata[as.numeric(allgraphdata$treat2)<as.numeric(allgraphdata$treat1),]
      allsedata<-allsedata[as.numeric(allsedata$treat2)<as.numeric(allsedata$treat1),]
    }
    
    # create graphs
    if(plot_type%in%c("fairylight") & effect=="treatment"){
      temp_graph<-
        allgraphdata%>%
        mutate(.model=factor(.model,
                             ordered = TRUE,
                             levels=names(graphsets)))%>%
        ggplot()+
        geom_vline(aes(xintercept = 0), lty=2)+
        geom_pointrange(aes(y=.model,
                            x=TE_d_NMR,
                            xmin=TE_d_NMR-(1.96*seTE_d_NMR),
                            xmax=TE_d_NMR+(1.96*seTE_d_NMR),
                            color=as.character(.model)))+
        facet_grid(treat1~treat2,scales = "free_x",switch = "y")+ # all below can change later
        scale_y_discrete(name=NULL)+
        scale_x_continuous(name="Treatment Effect Estimates")+
        scale_color_discrete(name="Model")+
        theme_bw()+
        theme(strip.placement = "outside",
              strip.text.y.left = element_text(angle = 0),
              axis.text.y = element_blank(),
              panel.spacing = unit(0, "in"))
      
    } else if(plot_type%in%c("fairylight") & effect=="interaction"){
      temp_graph<-
        allgraphdata%>%
        mutate(.model=factor(.model,
                             ordered = TRUE,
                             levels=names(graphsets)))%>%
        ggplot()+
        geom_vline(aes(xintercept = 0), lty=2)+
        geom_pointrange(aes(y=.model,
                            x=TE_beta_NMR,
                            xmin=TE_beta_NMR-(1.96*seTE_beta_NMR),
                            xmax=TE_beta_NMR+(1.96*seTE_beta_NMR),
                            color=as.character(.model)))+
        facet_grid(treat1~treat2,scales = "free_x",switch = "y")+ # all below can change later
        scale_y_discrete(name=NULL)+
        scale_x_continuous(name="Treatment Effect Estimates")+
        scale_color_discrete(name="Model")+
        theme_bw()+
        theme(strip.placement = "outside",
              strip.text.y.left = element_text(angle = 0),
              axis.text.y = element_blank(),
              panel.spacing = unit(0, "in"))
      
    } else if(plot_type%in%c("flashlight")){
      temp_graph<-allgraphdata%>%
        ggplot()+
        geom_ribbon(data=allsedata,
                    aes(x=xval,ymin=ciymin,
                        ymax=ciymax,fill=.model),alpha=0.2)+
        # geom_ribbon(data=allsedata,
        #             aes(x=xval,ymin=ciymin,
        #                 ymax=ciymax,fill=.model_comp),alpha=0.2)+
        geom_pointrange(aes(x=xcov,y=TE,ymin=TE-(1.96*seTE),ymax=TE+(1.96*seTE)))+
        geom_pointrange(aes(x=0,y=TE_d_NMR ,
                            ymin=TE_d_NMR -(1.96*seTE_d_NMR),
                            ymax=TE_d_NMR +(1.96*seTE_d_NMR),
                            color=.model),
                        shape=23)+
        geom_abline(aes(slope=TE_beta_NMR,intercept=TE_d_NMR,
                        color=.model),lty=2)+
        # geom_abline(aes(slope=TE_beta_NMR,intercept=TE_d_NMR,
        #                 color=.model_comp),lty=2)+
        scale_fill_discrete(name="Model")+
        scale_color_discrete(name="Model")+
        facet_grid(treat1~treat2, scales = "free_x")+# can be overwritten outside of function
        scale_y_continuous(name="Treatment Effect", sec.axis = sec_axis(~ . , name = "Comparator", breaks = NULL, labels = NULL))+
        scale_x_continuous(name="Covariate Value", sec.axis = sec_axis(~ . , name = "Reference", breaks = NULL, labels = NULL))+
        theme_bw()
    }
    # return(allgraphdata)
    return(temp_graph)
  }
}



light_plot.compadj.netmetareg<-function(netmeta_obj=NULL,
                      graphset=NULL, all_combinations=NULL, graphdata_se_area=NULL, # optional and only if we need interactions
                      confidence_color=NULL, 
                      #diamond_color=NULL,diamond_width=0.2,
                      se=TRUE,
                      interaction_color=NULL ,
                      Ival1=NULL,
                      Ival2=NULL,
                      add_interactions=FALSE,
                      xcov=NULL,
                      study_point_color=NULL, study_point_shape=NULL, 
                      plot_type="Independent Assymetry",
                      reference=NULL){
  if(is.null(netmeta_obj)){
    print("some error about missing reqs")
  } 
  
  if(!is.null(graphset)&(is.null(all_combinations)|is.null(graphdata_se_area))){
    all_combinations<-graphset$all_combinations
    graphdata_se_area<-graphset$graphdata_se_area
  }
  
  # Data Prep and Transformations ##################
  # get observation data in all directions
  all_obs<-netmeta_obj$data%>%mutate(.prec=seTE^-2)
  rev_obs<-all_obs
  
  nameroots<-names(rev_obs)[grepl("1$|2$",names(rev_obs))]
  nameroots<-table(gsub("1$|2$","",nameroots))
  rev_variables<-names(nameroots[nameroots==2])
  for(i in rev_variables) {
    assign(paste(i,"1temp",sep=""),rev_obs[[paste(i,1,sep="")]])
    rev_obs[[paste(i,1,sep="")]]<-rev_obs[[paste(i,2,sep="")]]
    rev_obs[[paste(i,2,sep="")]]<-get(paste(i,"1temp",sep=""))
  }
  rev_obs$TE<-rev_obs$TE*-1
  all_obs<-rbind(rev_obs,all_obs)
  all_obs<-netmeta_obj$TE.direct.random%>% # extract direct evidence estimates from NETMETA
    as.data.frame()%>%
    tibble::rownames_to_column(var = "treat1") %>%  # Add rownames as a column
    tidyr::pivot_longer(cols = -treat1, names_to = "treat2", values_to = "TE_d_NMA_direct")%>%
    filter(is.na(TE_d_NMA_direct)==FALSE)%>%
    left_join(all_obs, by=c("treat1", "treat2"))%>%
    mutate(adjTE=TE-TE_d_NMA_direct,
           .comp=paste(treat1,treat2,sep="_vs_"))
  
  if(grepl("UMIE",plot_type)){
    all_obs<-all_obs%>%
      filter(as.numeric({{Ival1}})-as.numeric({{Ival2}})==1) # comparator favored relative to reference (0 would mean no direction)
  }
  if(is.null(reference)==FALSE){
    all_obs<-all_obs%>%
      filter(treat2==reference) # comparator favored relative to reference (0 would mean no direction)
  }
  
  if(add_interactions==TRUE){
    graphdata_se_area<-all_obs%>%#removes self comparisons
      select(TE_d_NMA_direct,treat1,treat2)%>%
      distinct()%>%
      inner_join(graphdata_se_area, by=c("treat1","treat2"))%>%
      mutate(adjpred=pred-TE_d_NMA_direct,adjciymin=ciymin-TE_d_NMA_direct,adjciymax=ciymax-TE_d_NMA_direct)
    
    all_combinations<-all_obs%>%#removes self comparisons
      select(TE_d_NMA_direct,treat1,treat2)%>%
      distinct()%>%
      inner_join(all_combinations, by=c("treat1","treat2"))%>%
      mutate(adjTE_d_NMR=TE_d_NMR-TE_d_NMA_direct)
  }
  
  # Graph Elements ###############################
  temp_graph<-all_obs%>%
    ggplot()+
    geom_hline(aes(yintercept = 0),color="grey")
  
  # transform seTE area (graphdata_se_area)
  # add confidence region
  if(se==TRUE & add_interactions==TRUE){
    if(is.null(substitute(confidence_color))==TRUE|is.character(substitute(confidence_color))==TRUE){
      
      temp_graph<-temp_graph+
        geom_ribbon(data=graphdata_se_area,# optional filters
                    aes(x=xval,ymin=adjciymin,ymax=adjciymax, group=.comp), fill=ifelse(is.null(confidence_color),"black",confidence_color),alpha=0.2)
      
    } else if(is.name(substitute(confidence_color))) {
      temp_graph<-temp_graph+
        geom_ribbon(data=graphdata_se_area,# optional filters
                    aes(x=xval,ymin=adjciymin,ymax=adjciymax, group=.comp, fill=.data[[deparse(substitute(confidence_color))]]),alpha=0.2)
    }
  }
  
  # add study level observations
  if(is.null(substitute(study_point_color))==TRUE|is.character(substitute(study_point_color))==TRUE){
    if(is.null(substitute(study_point_shape))==TRUE|is.numeric(substitute(study_point_shape))==TRUE){
      temp_graph<-temp_graph+
        geom_point(aes(x={{xcov}},y=adjTE,size=.prec), color=ifelse(is.null(study_point_color),"black",study_point_color), shape=ifelse(is.null(study_point_shape),16,study_point_shape))
      
      temp_graph+
        geom_point(aes(x={{xcov}},y=adjTE,size=.prec), color=ifelse(is.null(study_point_color),"black",study_point_color), shape=ifelse(is.null(study_point_shape),16,study_point_shape))
      
    } else if(is.name(substitute(study_point_shape))) {
      
      temp_graph<-temp_graph+
        geom_point(aes(x={{xcov}},y=adjTE,size=.prec, shape=.data[[deparse(substitute(study_point_shape))]]), color=ifelse(is.null(study_point_color),"black",study_point_color))
    }
    
  } else if(is.name(substitute(study_point_color))) {
    if(is.null(substitute(study_point_shape))==TRUE|is.numeric(substitute(study_point_shape))==TRUE){
      temp_graph<-temp_graph+
        geom_point(aes(x={{xcov}},y=adjTE,size=.prec, color=.data[[deparse(substitute(study_point_color))]]), shape=ifelse(is.null(study_point_shape),16,study_point_shape))
    } else if(is.name(substitute(study_point_shape))) {
      temp_graph<-temp_graph+
        geom_point(aes(x={{xcov}},y=adjTE,size=.prec, shape=.data[[deparse(substitute(study_point_shape))]], color=.data[[deparse(substitute(study_point_color))]]))
    }
    
  }
  
  # if(is.null(substitute(study_point_color))==TRUE|is.character(substitute(study_point_color))==TRUE){
  #   if(is.null(substitute(study_point_shape))==TRUE|is.numeric(substitute(study_point_shape))==TRUE){
  #     temp_graph<-temp_graph+
  #       geom_pointrange(aes(x={{xcov}},y=adjTE,ymin=adjTE-(1.96*seTE),ymax=adjTE+(1.96*seTE)), color=ifelse(is.null(study_point_color),"black",study_point_color), shape=ifelse(is.null(study_point_shape),16,study_point_shape))
  #     
  #     temp_graph+
  #       geom_pointrange(aes(x={{xcov}},y=adjTE,ymin=adjTE-(1.96*seTE),ymax=adjTE+(1.96*seTE)), color=ifelse(is.null(study_point_color),"black",study_point_color), shape=ifelse(is.null(study_point_shape),16,study_point_shape))
  #     
  #   } else if(is.name(substitute(study_point_shape))) {
  #     
  #     temp_graph<-temp_graph+
  #       geom_pointrange(aes(x={{xcov}},y=adjTE,ymin=adjTE-(1.96*seTE),ymax=adjTE+(1.96*seTE), shape=.data[[deparse(substitute(study_point_shape))]]), color=ifelse(is.null(study_point_color),"black",study_point_color))
  #   }
  #   
  # } else if(is.name(substitute(study_point_color))) {
  #   if(is.null(substitute(study_point_shape))==TRUE|is.numeric(substitute(study_point_shape))==TRUE){
  #     temp_graph<-temp_graph+
  #       geom_pointrange(aes(x={{xcov}},y=adjTE,ymin=adjTE-(1.96*seTE),ymax=adjTE+(1.96*seTE), color=.data[[deparse(substitute(study_point_color))]]), shape=ifelse(is.null(study_point_shape),16,study_point_shape))
  #   } else if(is.name(substitute(study_point_shape))) {
  #     temp_graph<-temp_graph+
  #       geom_pointrange(aes(x={{xcov}},y=adjTE,ymin=adjTE-(1.96*seTE),ymax=adjTE+(1.96*seTE), shape=.data[[deparse(substitute(study_point_shape))]], color=.data[[deparse(substitute(study_point_color))]]))
  #   }
  #   
  # }
  
  
  # add interactions
  if(add_interactions==TRUE){
    if(is.null(substitute(interaction_color))==TRUE|is.character(substitute(interaction_color))==TRUE){
      temp_graph<-temp_graph+
        geom_abline(aes(slope=TE_beta_NMR,intercept=adjTE_d_NMR), color=ifelse(is.null(interaction_color),"black",interaction_color),lty=2, data=all_combinations) # Interactions as slopes
    } else if(is.name(substitute(interaction_color))) {
      temp_graph<-temp_graph+
        geom_abline(aes(slope=TE_beta_NMR,intercept=adjTE_d_NMR,color=.data[[deparse(substitute(interaction_color))]]),lty=2, data=all_combinations) # Interactions as slopes
    }
  }
  
  # faceting ###########
  if(plot_type%in%c("Independent UMIE Assymetry","Independent Assymetry")){# optional #UMIE or consistency (only for ind)
    temp_graph<-temp_graph+
      facet_grid(treat1~treat2)+# can be overwritten outside of function
      scale_y_continuous("reference favored <- y_{i}-unadjTE_{i} -> comparator favored",#name="Treatment Effect", 
                         sec.axis = sec_axis(~ . , name = "Comparator", breaks = NULL, labels = NULL))+
      scale_x_continuous(#name="Covariate Value", 
        name=deparse(substitute(xcov)),
        sec.axis = sec_axis(~ . , name = "Reference", breaks = NULL, labels = NULL))
  } else if(plot_type %in% c("Consistent Assymetry")){ #consistency (color coded for ind, black/white for common)
    temp_graph<-temp_graph+
      facet_grid(.~treat2)+# can be overwritten outside of function
      scale_x_continuous(#name="Covariate Value",
        name=deparse(substitute(xcov)),
        sec.axis = sec_axis(~ . , name = "Reference", breaks = NULL, labels = NULL))+
      scale_y_continuous("reference favored <- y_{i}-unadjTE_{i} -> comparator favored")
  } else if(plot_type=="UMIE Assymetry"){#UMIE (color coded by comparison!! for ind, black/white for common)
    facet_var_str <- as_label(enquo(Ival1))
    temp_graph<-temp_graph+
      facet_grid(.~.data[[facet_var_str]])+
      scale_y_continuous(name="Ival!=facetlabel favored <- y_{i}-unadjTE_{i} -> Ival=facetlabel favored")+ # revisit this label
      scale_x_continuous(#name="Covariate Value", 
        name=deparse(substitute(xcov)))
  }
  
  
  return(temp_graph)
}


ggfunnel.netmetareg<-function(netmeta_obj=NULL,
                   graphset=NULL, all_combinations=NULL, graphdata_se_area=NULL, # optional and only if we need interactions
                   confidence_color=NULL, 
                   se=TRUE,
                   interaction_color=NULL ,
                   Ival1=NULL,
                   Ival2=NULL,
                   add_interactions=FALSE,
                   xcov=NULL,
                   study_point_color=NULL, study_point_shape=NULL, study_point_size=1,
                   plot_type="Independent Assymetry",#c(Independent Assymetry, Independent UMIE Assymetry, Consistent Assymetry)
                   reference=NULL){
  if(is.null(netmeta_obj)){
    print("some error about missing reqs")
  } 
  
  if(!is.null(graphset)&(is.null(all_combinations)|is.null(graphdata_se_area))){
    all_combinations<-graphset$all_combinations
    graphdata_se_area<-graphset$graphdata_se_area
  }
  
  # Data Prep and Transformations ##################
  # get observation data in all directions
  all_obs<-netmeta_obj$data
  rev_obs<-all_obs
  
  nameroots<-names(rev_obs)[grepl("1$|2$",names(rev_obs))]
  nameroots<-table(gsub("1$|2$","",nameroots))
  rev_variables<-names(nameroots[nameroots==2])
  for(i in rev_variables) {
    assign(paste(i,"1temp",sep=""),rev_obs[[paste(i,1,sep="")]])
    rev_obs[[paste(i,1,sep="")]]<-rev_obs[[paste(i,2,sep="")]]
    rev_obs[[paste(i,2,sep="")]]<-get(paste(i,"1temp",sep=""))
  }
  rev_obs$TE<-rev_obs$TE*-1
  all_obs<-rbind(rev_obs,all_obs)
  all_obs<-netmeta_obj$TE.direct.random%>% # extract direct evidence estimates from NETMETA
    as.data.frame()%>%
    tibble::rownames_to_column(var = "treat1") %>%  # Add rownames as a column
    tidyr::pivot_longer(cols = -treat1, names_to = "treat2", values_to = "TE_d_NMA_direct")%>%
    filter(is.na(TE_d_NMA_direct)==FALSE)%>%
    left_join(all_obs, by=c("treat1", "treat2"))%>%
    mutate(adjTE=TE-TE_d_NMA_direct,
           ORadjTE=exp(adjTE), # if OR but not built yet
           # comp=paste(treat1,treat2,sep="_vs_"),
           comp=paste(pmin(treat1, treat2), pmax(treat1, treat2), sep="_vs_"))
  
  if(grepl("UMIE",plot_type)){
    all_obs<-all_obs%>%
      filter(as.numeric({{Ival1}})-as.numeric({{Ival2}})==1) # comparator favored relative to reference (0 would mean no direction)
  }
  if(is.null(reference)==FALSE){
    all_obs<-all_obs%>%
      filter(treat2==reference) # comparator favored relative to reference (0 would mean no direction)
  }
  if(add_interactions==TRUE){
    graphdata_se_area<-all_obs%>%#removes self comparisons
      select(TE_d_NMA_direct,treat1,treat2)%>%
      distinct()%>%
      inner_join(graphdata_se_area, by=c("treat1","treat2"))%>%
      mutate(adjpred=pred-TE_d_NMA_direct,adjciymin=ciymin-TE_d_NMA_direct,adjciymax=ciymax-TE_d_NMA_direct)
    
    all_combinations<-all_obs%>%#removes self comparisons
      select(TE_d_NMA_direct,treat1,treat2)%>%
      distinct()%>%
      inner_join(all_combinations, by=c("treat1","treat2"))%>%
      mutate(adjTE_d_NMR=TE_d_NMR-TE_d_NMA_direct)
  }
  
  # Graph Elements ###############################
  temp_graph<-all_obs%>%
    ggplot()+
    geom_vline(aes(xintercept = 0),color="grey")#exp(0) if OR
  
  # add confidence region ########
  if(se==TRUE & add_interactions==TRUE){
    if(is.null(substitute(confidence_color))==TRUE|is.character(substitute(confidence_color))==TRUE){
      
      temp_graph<-temp_graph+
        geom_ribbon(data=graphdata_se_area,# optional filters
                    aes(y=xval,xmin=adjciymin,xmax=adjciymax, group=comp), fill=ifelse(is.null(confidence_color),"black",confidence_color),alpha=0.2)
      
    } else if(is.name(substitute(confidence_color))) {
      temp_graph<-temp_graph+
        geom_ribbon(data=graphdata_se_area,# optional filters
                    aes(y=xval,xmin=adjciymin,xmax=adjciymax, group=comp, fill=.data[[deparse(substitute(confidence_color))]]),alpha=0.2)
    }
  }
  
  # add study level observations #######################
  if(is.null(substitute(study_point_color))==TRUE|is.character(substitute(study_point_color))==TRUE){
    if(is.null(substitute(study_point_shape))==TRUE|is.numeric(substitute(study_point_shape))==TRUE){
      temp_graph<-temp_graph+
        geom_point(aes(y={{xcov}},x=adjTE), size = study_point_size, color=ifelse(is.null(study_point_color),"black",study_point_color), shape=ifelse(is.null(study_point_shape),16,study_point_shape))
      
      temp_graph+
        geom_point(aes(y={{xcov}},x=adjTE), size = study_point_size, color=ifelse(is.null(study_point_color),"black",study_point_color), shape=ifelse(is.null(study_point_shape),16,study_point_shape))
      
    } else if(is.name(substitute(study_point_shape))) {
      
      temp_graph<-temp_graph+
        geom_point(aes(y={{xcov}},x=adjTE, shape=.data[[deparse(substitute(study_point_shape))]]), size = study_point_size, color=ifelse(is.null(study_point_color),"black",study_point_color))
    }
    
  } else if(is.name(substitute(study_point_color))) {
    if(is.null(substitute(study_point_shape))==TRUE|is.numeric(substitute(study_point_shape))==TRUE){
      temp_graph<-temp_graph+
        geom_point(aes(y={{xcov}},x=adjTE, color=.data[[deparse(substitute(study_point_color))]]), size = study_point_size, shape=ifelse(is.null(study_point_shape),16,study_point_shape))
    } else if(is.name(substitute(study_point_shape))) {
      temp_graph<-temp_graph+
        geom_point(aes(y={{xcov}},x=adjTE, shape=.data[[deparse(substitute(study_point_shape))]], size = study_point_size, color=.data[[deparse(substitute(study_point_color))]]))
    }
    
  }
  
  # if(is.null(substitute(study_point_color))==TRUE|is.character(substitute(study_point_color))==TRUE){
  #   if(is.null(substitute(study_point_shape))==TRUE|is.numeric(substitute(study_point_shape))==TRUE){
  #     temp_graph<-temp_graph+
  #       geom_pointrange(aes(x={{xcov}},y=adjTE,ymin=adjTE-(1.96*seTE),ymax=adjTE+(1.96*seTE)), color=ifelse(is.null(study_point_color),"black",study_point_color), shape=ifelse(is.null(study_point_shape),16,study_point_shape))
  #     
  #     temp_graph+
  #       geom_pointrange(aes(x={{xcov}},y=adjTE,ymin=adjTE-(1.96*seTE),ymax=adjTE+(1.96*seTE)), color=ifelse(is.null(study_point_color),"black",study_point_color), shape=ifelse(is.null(study_point_shape),16,study_point_shape))
  #     
  #   } else if(is.name(substitute(study_point_shape))) {
  #     
  #     temp_graph<-temp_graph+
  #       geom_pointrange(aes(x={{xcov}},y=adjTE,ymin=adjTE-(1.96*seTE),ymax=adjTE+(1.96*seTE), shape=.data[[deparse(substitute(study_point_shape))]]), color=ifelse(is.null(study_point_color),"black",study_point_color))
  #   }
  #   
  # } else if(is.name(substitute(study_point_color))) {
  #   if(is.null(substitute(study_point_shape))==TRUE|is.numeric(substitute(study_point_shape))==TRUE){
  #     temp_graph<-temp_graph+
  #       geom_pointrange(aes(x={{xcov}},y=adjTE,ymin=adjTE-(1.96*seTE),ymax=adjTE+(1.96*seTE), color=.data[[deparse(substitute(study_point_color))]]), shape=ifelse(is.null(study_point_shape),16,study_point_shape))
  #   } else if(is.name(substitute(study_point_shape))) {
  #     temp_graph<-temp_graph+
  #       geom_pointrange(aes(x={{xcov}},y=adjTE,ymin=adjTE-(1.96*seTE),ymax=adjTE+(1.96*seTE), shape=.data[[deparse(substitute(study_point_shape))]], color=.data[[deparse(substitute(study_point_color))]]))
  #   }
  #   
  # }
  
  
  # add interactions #################
  if(add_interactions==TRUE){
    if(is.null(substitute(interaction_color))==TRUE|is.character(substitute(interaction_color))==TRUE){
      temp_graph<-temp_graph+
        geom_abline(aes(slope=TE_beta_NMR,intercept=adjTE_d_NMR), color=ifelse(is.null(interaction_color),"black",interaction_color),lty=2, data=all_combinations) # Interactions as slopes
    } else if(is.name(substitute(interaction_color))) {
      temp_graph<-temp_graph+
        geom_abline(aes(slope=TE_beta_NMR,intercept=adjTE_d_NMR,color=.data[[deparse(substitute(interaction_color))]]),lty=2, data=all_combinations) # Interactions as slopes
    }
  }
  
  # faceting ###########
  if(plot_type%in%c("Independent UMIE Assymetry","Independent Assymetry")){# optional #UMIE or consistency (only for ind)
    temp_graph<-temp_graph+
      facet_grid(treat1~treat2)+# can be overwritten outside of function
      scale_y_reverse(deparse(substitute(xcov)),
                      sec.axis = sec_axis(~ . , name = "Comparator", breaks = NULL, labels = NULL))+
      scale_x_continuous("reference favored <- y_{i}-unadjTE_{i} -> comparator favored",#name="Treatment Effect", 
                         sec.axis = sec_axis(~ . , name = "Reference", breaks = NULL, labels = NULL))
  } else if(plot_type %in% c("Consistent Assymetry")){ #consistency (color coded for ind, black/white for common)
    temp_graph<-temp_graph+
      facet_grid(.~treat2)+# can be overwritten outside of function
      scale_x_continuous(name="reference favored <- y_{i}-unadjTE_{i} -> comparator favored",
                         sec.axis = sec_axis(~ . , name = "Reference", breaks = NULL, labels = NULL))+
      scale_y_reverse(name=deparse(substitute(xcov)))
  } else if(plot_type=="UMIE Assymetry"){#UMIE (color coded by comparison!! for ind, black/white for common)
    facet_var_str <- as_label(enquo(Ival1))
    temp_graph<-temp_graph+
      facet_grid(.~.data[[facet_var_str]])+
      scale_y_reverse(name=deparse(substitute(xcov)))+ # revisit this label
      scale_x_continuous(
        name="Ival!=facetlabel favored <- y_{i}-unadjTE_{i} -> Ival=facetlabel favored"
      )
  }
  
  
  return(temp_graph)
}
