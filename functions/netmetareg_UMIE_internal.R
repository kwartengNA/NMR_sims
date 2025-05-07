# result output ###########

as_df <- function(x,
                  covar = x$.netmeta$covar,
                  covar.name = x$.netmeta$covar.name) {
  
  treat_key<-data.frame("treat_long"=x$.netmeta$trts, "treat"=x$.netmeta$trts.abbr) # connect abbreviated trts to the original treatment names
  treat_key2<-treat_key
  treat_key2$treat <- ifelse(treat_key2$treat== x$.netmeta$reference.group, "nonref", treat_key2$treat) # Not intuitive. Common with consistency only incorporates covariates from reference containing comparisons, so all non-ref contain the observed covariate levels of the reference
  treat_key<-unique(rbind(treat_key2,treat_key))
  
  res_b_df <- data.frame(treat = "", est = x$b, se = x$se,
                         lower = x$ci.lb, upper = x$ci.ub,
                         z = x$zval, pval = x$pval)
						 
  rownames(res_b_df)<-gsub("^`|`$","",rownames(res_b_df)) # remove backticks which are created in the manual model matrix
  
  #
  res_b_df$treat <- rnam <- rownames(res_b_df)
  #
  sel_TE <- !grepl(":", rnam)
  res_b_TE <- res_b_df[sel_TE, , drop = FALSE]
  rnam_d <- paste0("d[", rnam[sel_TE], "]")
  #
  if (is.null(covar)) {
    rownames(res_b_TE) <- rnam_d
    return(res_b_TE)
  }
  
  rnam_beta <- paste0("beta[", rnam[!sel_TE], "]")
  #
  res_b_TE$type <- "d"
  res_b_TE$covar <- NA
  res_b_TE$cov_lvl <- NA
  res_b_TE$cov_ref <- NA
  res_b_TE$treat_long<-NA # retain information on original treatment name
  
  # Get estimates covariate cov_lvls for each treatment
  res_b_interactions <- res_b_df[!sel_TE, , drop = FALSE]
  # Split 'rowname' using strsplit
  split_names <- do.call(rbind, strsplit(rownames(res_b_interactions), "[:]"))
  #
  res_b_interactions$type <- "beta"
  res_b_interactions$treat <- split_names[ , 1]
  res_b_interactions$covar <- covar.name#split_names[, 2]
  res_b_interactions<-merge(res_b_interactions, treat_key, by=c("treat"),all.x = TRUE) # add original treatment names
  
  
  is_num <- is.numeric(covar)
  #
  if (!is_num) {
    # Get observed covariate levels for each treatment
    #
    dat <- x$.netmeta$x$data[, c("treat1","treat2", covar.name)]
    #
    # Convert to long format
    #
    unique_lvls <- rbind(setNames(dat[,c("treat1",covar.name)],c("treat_long","lvl")), # since this data frame comes from the original dataset, it relies on the long/original treatment names
                         setNames(dat[,c("treat2",covar.name)],c("treat_long","lvl")))
    unique_lvls <- unique(unique_lvls)
    unique_lvls<-merge(unique_lvls, data.frame(treat_long=x$.netmeta$trts, treat=x$.netmeta$trts.abbr), by=c("treat_long"), all.x = TRUE) # add corresponding abbreviations
    #
    if (x$.netmeta$assumption == "common") {
      unique_lvls$treat <- ifelse(unique_lvls$treat == x$.netmeta$reference.group,"nonref",unique_lvls$treat)
    }
    # capture estimated factor levels
    res_b_interactions$lvl <- gsub(covar.name,"",split_names[, 2])
    
    # Get the non-included covariate levels aka covariate references
    # aggregate(lvl ~ treat_long+treat, data = unique_lvls, FUN = min)
    ref_lvls <-
      merge(unique_lvls, res_b_interactions, by = c("treat_long","lvl", "treat"),
            all.x = TRUE)
    ref_lvls <- ref_lvls[is.na(ref_lvls$est) == TRUE,c("treat_long","lvl", "treat")]
    colnames(ref_lvls)[colnames(ref_lvls) == "lvl"] <- "cov_ref"
    res_b_interactions <-
      merge(res_b_interactions, ref_lvls, by = c("treat_long", "treat"), all.x = TRUE) # should be many to 1. If there is an issue then it is likely in rma.mv
    colnames(res_b_interactions)[colnames(res_b_interactions) == "lvl"] <- "cov_lvl"
    res_b_interactions<-res_b_interactions[,names(res_b_interactions)!="treat_long"] # remove treat_long since we just used it for internal merging
  }
  #
  res <- rbind(res_b_TE[, names(res_b_interactions)], res_b_interactions)
  #
  rownames(res) <- c(rnam_d, rnam_beta)
  
  attr(res, "reference.group") <- x$.netmeta$reference.group
  attr(res, "covar.name") <- covar.name
  #
  res <- res[, c("type", "treat",
                 if (!is_num) "cov_lvl", if (!is_num) "cov_ref",
                 "est", "se", "lower", "upper", "z", "pval")]
  #
  res
}
# create directionality parameter as treatment order ################
mkIval<-function(x, ref){ # treatment order as direction
  if(ref %in% c(x$treat1,x$treat2)){
    x$Ival1<-ifelse(x$treat1==ref, 0, 1)
    x$Ival2<-ifelse(x$treat2==ref, 0, 1)
  } else {
    x$Ival1<-0
    x$Ival2<-1
  }
  return(x)
} 

# check Ival. maybe add ability to specifiy 2 custom variables #############
# check categorical covariate distribution #############
chk_covarlvl<-function(x, consistency, covar.name){
  if(consistency==TRUE){
    x_long<-rbind(setNames(x[, c("treat1", covar.name)],c("treat", covar.name)),#.comp_
        setNames(x[, c("treat2", covar.name)],c("treat", covar.name)))
    x_long<-unique(x_long)
    unavail_lvls<-colSums((table(x_long$treat, x_long[[covar.name]])==0))==0
    
  } else if(consistency==FALSE){
    x_long<-x[, c(".comp_", covar.name)]
    x_long<-unique(x_long)
    unavail_lvls<-colSums((table(x_long$treat, x_long[[covar.name]])==0))==0
    }
  
  if(sum(unavail_lvls)==0){
    stop("Error: No common levels found in the categorical covariate across groups. Please check for missing levels across groups")
  }
  
  if(names(unavail_lvls[unavail_lvls!=0]) ==min(levels(x_long[[covar.name]]))){
    stop("Error: The reference level of the covariate is missing in some comparisons. Please ensure the variable is a factor and define an appropriate reference level.")
  }
  # otherwise proceed
}
