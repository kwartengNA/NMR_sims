#' Network meta-regression with a single continuous or binary covariate
#'
#' @description
#' Network meta-regression with a single continuous or binary covariate for
#' objects of class \code{netmeta}. This is a wrapper function for the R
#' function \code{\link[metafor]{rma.mv}} in the R package \bold{metafor}
#' (Viechtbauer 2010).
#'
#' @details
#' This R function is a wrapper function for R function
#' \code{\link[metafor]{rma.mv}} in the R package \bold{metafor}
#' (Viechtbauer 2010).
#'
#' Note, results are not back-transformed in printouts of
#' network meta-analyses using summary measures with transformations, e.g.,
#' log risk ratios are printed instead of the risk ratio if argument
#' \code{sm = "RR"}.
#'
#' Argument '\dots{}' can be used to pass additional arguments to R
#' function \code{\link[metafor]{rma.mv}}. For example, argument
#' \code{control} to provide a list of control values for the
#' iterative estimation algorithm. See help page of R function
#' \code{\link[metafor]{rma.mv}} for more details.
#'
#' @param x An object of class \code{netmeta}.
#' @param covar Continuous or binary covariate.
#' @param consistency A logical indicating whether a consistency or
#'   inconsistency model should be assumed.
#' @param assumption A character string indicating which assumption is done
#'   for the covariate; either "independent" or "common" (can be abbreviated).
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance tau-squared. Either
#'   \code{"FE"}, \code{"REML"}, or \code{"ML"}.
#' @param level The level used to calculate confidence intervals for regression
#'   coefficients.
#' @param reference.group Reference treatment.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard
#'   errors.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-value, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param print.se A logical specifying whether standard errors should be
#'   printed.
#' @param details.methods A logical specifying whether details on statistical
#'   methods should be printed.
#' @param \dots Additional arguments passed to R function
#'   \code{\link[metafor]{rma.uni}}.
#'
#' @return
#' An object of class \code{c("netmetareg", "rma.mv", "rma")}. Please
#' look at the help page of R function \code{\link[metafor]{rma.mv}}
#' for more details on the output from this function.
#'
#' In addition, a list \code{.netmeta} is added to the output containing
#' the following components:
#' \item{x, covar, method.tau}{As defined above.}
#' \item{dots}{Information provided in argument '\dots{}'.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{netmeta} used to create
#'   object.}
#' \item{version.metafor}{Version of R package \bold{metafor} used to
#'   create object.}
#'
#' @author Nana-Adjoa Kwarteng
#'   \email{nana-adjoa.kwarteng@uniklinik-freiburg.de},
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
#' @keywords models regression
#' 
#' @examples
#' # Add examples
#' 
#' @rdname netmetareg
#' @method netmetareg netmeta
#' @export 

netmetareg.netmeta <- function(x, covar = NULL,
                               consistency = TRUE,
                               assumption = "independent",
                               method.tau = if (!x$random) "FE" else "REML",
                               max_interactions=FALSE,
                               level = x$level.ma,
                               reference.group = x$reference.group,
                               nchar.trts = x$nchar.trts, ...) {
  
  chkclass(x, "netmeta")
  #
  x <- updateversion(x)
  
  #
  # Checks and assignments
  #
  chklogical(consistency)
  assumption <- setchar(assumption, c("independent", "common"))
  #
  reference.group <- setref(reference.group, x$trts)
  method.tau <- setchar(method.tau, c("REML", "ML", "FE"))
  chklevel(level)
  #
  sm <- x$sm
  
  #
  # Return network meta-analysis object if covariate is missing
  #
  
  if (missing(covar)) {
    warning("No network meta-regresssion conducted as argument 'covar'",
            "is missing.")
    return(x)
  }
  #
  data <- x$data
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  if (is.null(data))
    data <- sfsp
  #
  # Catch covariate
  #
  covar.name <- gsub("\"", "", deparse(substitute(covar)), fixed = TRUE)
  covar.name <- gsub("$", ".", covar.name, fixed = TRUE)
  covar.name <- gsub("+", "_", covar.name, fixed = TRUE)
  covar.name <- gsub("-", "_", covar.name, fixed = TRUE)
  covar.name <- gsub(":", "_", covar.name, fixed = TRUE)
  #
  covar <- catch("covar", mc, data, sfsp)
  #
  trts <- x$trts
  trts.abbr <- treats(x$trts, nchar.trts)
  #
  reference.group <- trts.abbr[trts == reference.group]
  
  dat <- data.frame(studlab = x$data$.studlab,
                    treat1 = x$data$.treat1,
                    treat2 = x$data$.treat2,
                    TE = x$data$.TE,
                    seTE = x$data$.seTE)
  
  
  #
  dat$treat1 <-
    as.character(factor(dat$treat1, levels = trts, labels = trts.abbr))
  dat$treat2 <-
    as.character(factor(dat$treat2, levels = trts, labels = trts.abbr))
  #
  dat$treat1 <- gsub(" ", "_", dat$treat1, fixed = TRUE)
  dat$treat2 <- gsub(" ", "_", dat$treat2, fixed = TRUE)
  #
  if (length(unique(c(dat$treat1, dat$treat2))) != length(x$trts))
    stop("Number of treatments changes after replacing whitespaces with '_'.",
         call. = FALSE)
  # #
  # if (is.factor(covar) && length(levels(covar)) > 2)
  #   stop("Network meta-regression not available for a covariate of ",
  #        "class 'factor' with more than two categories.",
  #        call. = FALSE)
  # #
  # if (is.logical(covar)) # deprecated because model matrix is now constructed manually
  #  covar <- as.numeric(covar)
  #
  dat[[covar.name]] <- covar
  #
  if (!is.null(x$data$.n1) & !is.null(x$data$.n2)) { #netmeta is not exporting n1 and n2 for ORs
    dat$n1 <- x$data$.n1
    dat$n2 <- x$data$.n2
  }
  #
  if (!is.null(x$data$.mean1) & !is.null(x$data$.mean2)) {
    dat$mean1 <- x$data$.mean1
    dat$mean2 <- x$data$.mean2
  }
  #
  if (!is.null(x$data$.sd1) & !is.null(x$data$.sd2)) {
    dat$sd1 <- x$data$.sd1
    dat$sd2 <- x$data$.sd2
  }
  #
  if (!is.null(x$data$.time1) & !is.null(x$data$.time2)) {
    dat$time1 <- x$data$.time1
    dat$time2 <- x$data$.time2
  }

  # the logic above does not work because the following variables do not have .variable equivalents exported from netmeta. They are needed for UMIE models and/or the calcV function
  if (!is.null(x$data$event1) & !is.null(x$data$event2)) {
    dat$event1 <- ifelse(x$data$treat1==x$data$.treat1, x$data$event1, x$data$event2)
    dat$event2 <- ifelse(x$data$treat2==x$data$.treat2, x$data$event2, x$data$event1)
  }
  if (!is.null(x$data$p1) & !is.null(x$data$p2)) {
    dat$p1 <- ifelse(x$data$treat1==x$data$.treat1, x$data$p1, x$data$p2)
    dat$p2 <- ifelse(x$data$treat2==x$data$.treat2, x$data$p2, x$data$p1)
  }
  if (!is.null(x$data$n1) & !is.null(x$data$n2)) {
    dat$n1 <- ifelse(x$data$treat1==x$data$.treat1, x$data$n1, x$data$n2)
    dat$n2 <- ifelse(x$data$treat2==x$data$.treat2, x$data$n2, x$data$n1)
  }
  if (!is.null(x$data$incr)) {
    dat$incr <- x$data$incr
  }

  # RETURN need more comprehensive logic/ability to supply Ival for UMIE
  # if Ival not supplied default to treatment order
  
  
  # if Ival supplied prior to netmeta. adapt if we allow custom variable
  if (consistency==FALSE & !is.null(x$data$Ival1) & !is.null(x$data$Ival2)){ 
    # must be coded 0.1. need to add chk for Ival consistency v
    dat$Ival1 <- ifelse(x$data$treat1==x$data$.treat1, x$data$Ival1, x$data$Ival2) # currently Ival must be defined before using netmeta. I think all multi arm studies might need to have the same reference treatment. Investigate. or should we allow users to provide individual names/ adapt the internal function mkIval
    dat$Ival2 <- ifelse(x$data$treat2==x$data$.treat2, x$data$Ival2, x$data$Ival1)
    } else if (consistency==FALSE & (is.null(x$data$Ival1) & is.null(x$data$Ival2))){
      
      mkIval(dat, ref = reference.group) # in netmetareg internal
      dat<-do.call(rbind, lapply(split(dat, dat$studlab), mkIval, ref=reference.group))
      rownames(dat)<-NULL
      warning("No Directionality Parameter Specified. Defaulting to treatment order")
    }
  #
  dat <- dat[order(dat$studlab, dat$treat1, dat$treat2), , drop = FALSE]
  
  #
  keep <- logical(0)
  wo <- logical(0)
  #
  if(consistency==TRUE){
    # check for constant covariate in reference treatment
    if(length(unique(dat[dat$treat1==reference.group|dat$treat2==reference.group,covar.name]))==1){
      stop("Invalid reference treatment for interaction. Insufficient variation in observed covariate values for reference treatment")
      }
    for (i in unique(dat$studlab)) {
      d.i <- dat[dat$studlab == i, , drop = FALSE]
      trts.i <- unique(sort(c(d.i$treat1, d.i$treat2)))
      if (reference.group %in% trts.i)
        ref.i <- reference.group
      else
        ref.i <- rev(trts.i)[1]
      #
      keep.i <- !(d.i$treat1 != ref.i & d.i$treat2 != ref.i)
      wo.i <- d.i$treat1 == ref.i
      #
      keep <- c(keep, keep.i)
      wo <- c(wo, wo.i)
    }
  } else  if(consistency==FALSE){
    # calcV is sensitive to the choice of reference treatment. but if we have study specific references we had to update the keep and wo statements
    dat$Ival_diff<-dat$Ival1-dat$Ival2
    keep.i <- dat$Ival_diff!=0 # whether a row contains a reference. equivalent !(d.i$treat1 != ref.i & d.i$treat2 != ref.i)
    wo.i <- dat$Ival1==0 # whether the reference treatment is treat1. equivalent to  d.i$treat1 == ref.i
    #
    keep <- c(keep, keep.i)
    wo <- c(wo, wo.i)
  }
  #
  dat <- dat[keep, , drop = FALSE] 
  #
  wo <- wo[keep] 
  
  if (sum(wo) > 0) {
    dat$TE[wo] <- -dat$TE[wo]
    #
    t2.i <- dat$treat2
    e2.i <- dat$event2
    n2.i <- dat$n2
    mean2.i <- dat$mean2
    sd2.i <- dat$sd2
    time2.i <- dat$time2
    Ival2.i <- dat$Ival2 
    #
    dat$treat2[wo] <- dat$treat1[wo]
    dat$event2[wo] <- dat$event1[wo]
    dat$n2[wo] <- dat$n1[wo]
    dat$mean2[wo] <- dat$mean1[wo]
    dat$sd2[wo] <- dat$sd1[wo]
    dat$time2[wo] <- dat$time1[wo]
    dat$Ival2[wo] <- dat$Ival1[wo]
    #
    dat$treat1[wo] <- t2.i[wo]
    dat$event1[wo] <- e2.i[wo]
    dat$n1[wo] <- n2.i[wo]
    dat$mean1[wo] <- mean2.i[wo]
    dat$sd1[wo] <- sd2.i[wo]
    dat$time1[wo] <- time2.i[wo]
    dat$Ival1[wo] <- Ival2.i[wo]
  }
  #
  ncols1 <- ncol(dat)
  dat <- contrmat(dat, grp1 = "treat1", grp2 = "treat2")
  ncols2 <- ncol(dat)
  varnames <- names(dat)[(ncols1 + 1):ncols2]
  #
  dat <- dat[order(dat$studlab), ]
  #
  trts.all <- varnames
  trts <- trts.all[trts.all!=reference.group] # consistency
  #
  # Calculate Variance-Covariance matrix
  #
  V <- bldiag(lapply(split(dat, dat$studlab), calcV, sm = sm))
  
  #
  
  if (consistency) {
    # trts <- trts.all[trts.all!=reference.group] # consistency
    if (assumption == "independent") {
      formula.nmr_default <- # renamed for construction of manual matrix below
        as.formula(paste("~ 0 + ",
                         paste(trts, collapse = " + "),
                         if (!is.null(covar))
                           paste0( " + ",
                             paste(paste0(trts, ":", covar.name),
                                 collapse = " + "))))
    } else {
      #
      # The interaction is all non-reference treatments vs reference,
      # so we can define a variable indicating whether it is reference
      # or non-reference.
      # Then get the interaction using the colon(:). The output has a
      # colon which is in the same format as the independent. Helpful
      # for later extraction
      #
      dat$nonref <- as.numeric(dat[[reference.group]] != 0) # is the reference group set above?
      #
      formula.nmr_default <- # renamed for construction of manual matrix below
        as.formula(paste("~ 0 + ",
                         paste(trts, collapse = " + "),
                         if (!is.null(covar))
                           paste0( " + ",
                                   paste(paste0("nonref:", covar.name),
                                         collapse = " + "))))
    }
  } else {
    # trts <- varnames[-length(varnames)]
       dat$Ival<-dat$Ival1-dat$Ival2 #basically Ivaldiff updated after "wo"
    if (assumption == "independent") {
      dat$.comp_<-paste(pmin(dat$treat1, dat$treat2), pmax(dat$treat1, dat$treat2), sep="_vs_")
      count_comp<-table(dat$.comp_)
      if(max_interactions==FALSE){ # for dev use
        dat$.comp_<-ifelse(dat$.comp_%in% names(count_comp[count_comp>=2]),dat$.comp_, "insufficient_data") # requires at least 2 observations
        }
      formula.nmr_default <-as.formula(paste("~ 0 + ",
                                             paste(trts, collapse = " + "),
                                             if (!is.null(covar)){
                                               paste0("+ .comp_:Ival:",covar.name)}
      ))
      
    } else {
      formula.nmr_default <-as.formula(paste("~ 0 + ",
                                             paste(trts, collapse = " + "),
                                             if (!is.null(covar)){
                                               paste0("+ Ival:",covar.name)}
                                             ))
    }
    
    
  }
  
  
  man_matrix<-model.matrix(formula.nmr_default, data=dat) # the default model matrix
  man_matrix<-man_matrix[,grepl("insufficient_data", colnames(man_matrix))==FALSE]
  # checks for non-continuous covariates
  if((is.factor(covar)|is.character(covar)) && suppressWarnings(length(unique(covar))>2)){
    # chk_covarlvl(x=dat, consistency = consistency, covar.name = covar.name) # in netmetareg internal
  }
  if(is.factor(covar)){# check type so that we know how to handle the specific model matrix
    man_matrix<-man_matrix[,grepl(paste(min(levels(covar)),"$", sep=""), colnames(man_matrix))==FALSE] # drop interactions which contain the reference covariate level if covariate is factor, character, logical
  } else if(is.character(covar)){ 
    man_matrix<-man_matrix[,grepl(paste(min(covar),"$", sep=""), colnames(man_matrix))==FALSE] 
  } else if(is.logical(covar)){
    man_matrix<-man_matrix[,grepl(paste(as.logical(min(covar)),"$", sep=""), colnames(man_matrix))==FALSE] # 
  }
  dat[colnames(man_matrix)]<-man_matrix ## update model matrix columns
  formula.nmr<- as.formula(paste("~ 0 ",
                                 paste0( " + ",paste("`",colnames(man_matrix),"`",sep="", ## the back ticks are necessary if the covariate levels/values contains spaces or special characters
                                                     collapse = " + "))))
  
  # Get rid of warning 'Undefined global functions or variables'
  treat1 <- treat2 <-  comparison <- NULL
  
  #
  # Covariate 'x' makes problems without removing network meta-analysis object x
  #
  ..x <- x
  rm(x)
  #
  dat.TE <- dat$TE
  dat$comparison <- paste(dat$treat1, dat$treat2, sep = " vs ")
  #
  suppressWarnings(res <-
    runNN(rma.mv,
          list(yi = dat.TE, V = V, 
               data = dat,
               mods = formula.nmr, 
               random = as.call(~ factor(comparison) | studlab),
               rho = 0.5,
               method = method.tau,
               level = 100 * level,
               ...)))
  #
  X <- res$X.f
  rownames(X) <- dat$studlab
  X <- X[rev(do.call(order, as.data.frame(X))), , drop = FALSE]
  #
  res$.netmeta <- list(x = ..x,
                       covar.name = covar.name,
                       covar = covar,
                       X = X,
                       consistency = consistency,
                       assumption = assumption,
                       method.tau = method.tau,
                       level = level,
                       reference.group = reference.group,
                       trts = ..x$trts,  # keep the long name of the treatment
                       trts.abbr = trts.abbr,
                       nchar.trts = nchar.trts,
                       dots = list(...),
                       call = match.call(),
                       version = packageDescription("netmeta")$Version,
                       version.metafor = packageDescription("metafor")$Version)
  # #
  res$results <- as_df(res)
  # #
  class(res) <- c("netmetareg", class(res))
  res$call <- NULL
  #
  res
}


#' @rdname netmetareg
#' @export netmetareg

netmetareg <- function(x, ...)
  UseMethod("netmetareg")


#' @rdname netmetareg
#' @method netmetareg default
#' @export

netmetareg.default <- function(x, ...)
  stop("Network meta-regression not available for an object of class '",
       class(x)[1], "'.",
       call. = FALSE)


#' @rdname netmetareg
#' @method print netmetareg
#' @export

print.netmetareg <- function(x,
                             digits = gs("digits"),
                             digits.se = gs("digits.se"),
                             digits.stat = gs("digits.stat"),
                             digits.pval = gs("digits.pval"),
                             print.se = FALSE,
                             details.methods = TRUE, ...) {
  
  chkclass(x, "netmetareg")
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 0, length = 1)
  chklogical(print.se)
  chklogical(details.methods)
  #
  covar <- x$.netmeta$covar
  #
  if (is.null(covar))
    cat("Network Meta-Analysis\n\n")
  else
    cat("Network Meta-Regression\n\n")
  #
  x.orig <- x
  #
  consistency <- x$.netmeta$consistency
  assumption <- x$.netmeta$assumption
  #
  # class(x) <- class(x)[class(x) != "netmetareg"]
  ##
  #print(x, ...)
  #
  lower <- NULL
  #
  if (is.numeric(covar)|is.integer(covar)|is.null(covar)) # update print output for as_df
    dat <- x$results[ , c("coef", if (print.se) "se",
                          "lower", "upper", "z", "pval")]
  else
    dat <- x$results[ , c("cov_lvl", # we need the covariate level for all assumptions if the covariate is not continuous. theoretically, we might not need them for logical vectors, because logical vectors can be coerced to continuous
                          "cov_ref",
                          "coef", if (print.se) "se",
                          "lower", "upper", "z", "pval")]
    #
  dat$coef <- formatN(dat$coef, digits = digits)
  if (print.se)
  dat$se <- formatN(dat$se, digits = digits.se, ...)
  dat$lower <- formatCI(formatN(dat$lower, digits = digits),
                        formatN(dat$upper, digits = digits))
  dat$upper <- NULL
  #
  dat <- replaceNA(dat, ".")
  #
  names(dat)[names(dat) == "lower"] <-
    paste0(round(100 * x$.netmeta$level, 1), "%-CI")
  #
  dat$z <- formatN(dat$z, digits = digits.stat)
  dat$pval <- formatPT(dat$pval, digits = digits.pval)
  names(dat)[names(dat) == "pval"] <- "p-value"
  #
  prmatrix(dat, quote = FALSE, right = TRUE, ...)
  #
  # Get rid of warning 'Undefined global functions or variables'
  comparison <- NULL
  #
  if (!is.null(covar) & details.methods) {
    cat("\nDetails on network meta-regression methods:\n")
    #
    if (x$.netmeta$method.tau == "FE")
      cat("- Common effects model\n")
    else
      cat("- Random effects model\n")
    #
    if (consistency)
      cat(paste0("- ", if (consistency) "C" else "Inc", "onsistency model\n"))
    #
    if (assumption == "independent")
      cat("- Independent slopes\n")
    else
      cat("- Common slope\n")
    #
    cat(paste0("- Reference group: ", x$.netmeta$reference.group, "\n"))
  }
  else if (details.methods) {
    cat("\nDetails on network meta-analysis methods:\n")
    #
    if (x$.netmeta$method.tau == "FE")
      cat("- Common effects model\n")
    else
      cat("- Random effects model\n")
  }
  #
  invisible(NULL)
}
