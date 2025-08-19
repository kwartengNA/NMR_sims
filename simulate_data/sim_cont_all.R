# Data Generation Setup  ############################

# Check if input values are assigned (i.e. from loop), or assign defaults   ############################
seedval      <- if (!exists("seedval")) sample(1:1000, 1) else seedval
jagsseedval  <- if (!exists("jagsseedval")) sample(1:1000, 1) else jagsseedval
set.seed(seedval)

# Assumptions
true_consistency <- if (!exists("true_consistency")) TRUE else true_consistency
true_assumption <- if (!exists("true_assumption")) "independent" else true_assumption
interaction_direction <- if (!exists("interaction_direction")) NA else interaction_direction

# Arms
true_arms <- if (!exists("true_arms")) "twoarm" else true_arms

# Define other simulation parameters   ############################
N.sim <- if (!exists("N.sim")) 2 else N.sim   # number of simulations
keep_prop <- if (!exists("keep_prop")) 1 else keep_prop
tau <- if (!exists("tau")) sqrt(runif(1, 0, 1)) else tau

NT <- 7   # number of treatments
NS <- 3   # number of studies per comparison (NS*NT for multi-arm)

# Patient allocation
Npmin <- 30
Npmax <- 60

# Effect sizes
es_MD <- c(0, seq(1/(NT-1), 1, length.out = NT-1))

# Interaction Size/assumptions
beta_range <- 0.75

# Covariate range
maxxrange <- 5
minxrange <- -1


# Generate Interaction Matrix   ############################

generate_beta_x <- function(NT, assumption, consistency, interaction_direction, beta_range) {
  if (consistency) {
    # Consistent interactions
    if (assumption == "common") {
      beta_x_old <- c(0, rep(runif(1, -beta_range, beta_range), NT-1))
    } else if (assumption == "independent") {
      beta_x_old <- c(0, runif(NT-1, -beta_range, beta_range))
    } else if (assumption == "NMA") {
      beta_x_old <- rep(runif(1, -beta_range, beta_range), NT)
    }
    
    if (!is.na(interaction_direction)) {
      if (interaction_direction == "positive") beta_x_old <- abs(beta_x_old)
      if (interaction_direction == "negative") beta_x_old <- -abs(beta_x_old)
    }
    
    beta_x <- outer(beta_x_old, beta_x_old, FUN = "-") * -1
    
  } else {
    # UMIE interactions
    if (assumption == "common") {
      vals <- rep(runif(1, -beta_range, beta_range), NT^2)
    } else if (assumption == "independent") {
      vals <- runif(NT^2, -beta_range, beta_range)
    } else if (assumption == "NMA") {
      vals <- rep(0, NT^2)
    }
    beta_x <- matrix(vals, nrow = NT, ncol = NT)
    beta_x[lower.tri(beta_x)] <- -beta_x[upper.tri(beta_x)]
  }
  
  diag(beta_x) <- 0
  beta_x
}

beta_x <- generate_beta_x(NT, true_assumption, true_consistency, interaction_direction, beta_range)

# Generate Network  ############################
if (true_arms == "multiarm") {
  source("simulate_data/sim_cont_multiarm.R") # all generations have the same structure
  dat.studies_long<-data2[[1]]%>%
    select(study, treat, treatment)
  
} else if (true_arms == "twoarm") {
  source("simulate_data/sim_cont2arm.R")
  dat.studies_long<-data2[[1]]%>% # all generations have the same structure
    select(study, treat, treatment)
  
} else if (true_arms == "mixarm") {
  
  if (!exists("dat.studies_merge")) {source("simulate_data/sim_mixarm_gs.R")} # define structure. creates dat.studies long in first mixarm/density configuration
  
  # Two-arm values
  source("simulate_data/sim_cont2arm.R")
  arm2_data2 <- data2
  arm2_data1 <- data1
      
  # Multi-arm values
  NS_bkup <- NS
  NS <- round(ctmixarm(NT)$nmulti)
  source("simulate_data/sim_cont_multiarm.R")
  armmulti_data2 <- data2
  armmulti_data1 <- data1
  NS <- NS_bkup
  
  # combine values. Apply structure in next step with keep studies subset
  data2 <- mapply(function(df1, df2) {
    dplyr::bind_rows(
      dplyr::mutate(df1, arms = "twoarm"),
      dplyr::mutate(df2, arms = "multiarm"))%>%
      group_by(study, arms) %>%
      mutate(comparison2 = ifelse(arms == "multiarm", "multiarm", paste(unique(treat), collapse = " vs "))) %>%
      ungroup() %>%
      arrange(comparison2, study) %>%
      group_by(comparison2, study) %>%
      mutate(rank = cur_group_id()) %>%
      ungroup() %>%
      group_by(comparison2) %>%
      mutate(rank = rank - min(rank) + 1) %>%
      ungroup()
    }, arm2_data2, armmulti_data2, SIMPLIFY = FALSE)
  
  }

# Keep Studies Subset ############################

if (!exists("dat.studies_merge")) dat.studies_merge <- NULL

if (is.null(dat.studies_merge)) {
  while (sum(dat.studies_merge$treat == 1) <= 1) {
    alltwoarm_studyids <- dat.studies_long %>% count(study) %>% filter(n == 2)
    allmultiarm <- dat.studies_long %>% count(study) %>% filter(n > 2)
    
    keep_studies <- sample(
      alltwoarm_studyids$study,
      floor(keep_prop * length(unique(dat.studies_long$study))) - nrow(allmultiarm),
      replace = FALSE
    )
    keep_studies <- c(keep_studies, allmultiarm$study)
    dat.studies_merge <- dat.studies_long %>% filter(study %in% keep_studies)
  }
  rm(allmultiarm, alltwoarm_studyids, keep_studies)
}

if (true_arms == "mixarm") {
  # can't join on study id
  for (i in 1:N.sim) {
    data2[[i]] <- dat.studies_merge %>%
      left_join(data2[[i]], by = c("treat","comparison2","rank")) %>%
      arrange(study.x, treat) %>%
      select(study = study.x, study.y, treat, comparison = comparison2,
             n, #events, 
             xcov, treatment, Ival_std, treat.ref)
  }

} else {
  for (i in 1:N.sim) {
    data2[[i]] <- dat.studies_merge %>%
      left_join(data2[[i]], by = c("treat","study", "treatment")) %>%
      arrange(study, treat) %>%
      select(study, treat, 
             n, yi, sd, 
             xcov, treatment, Ival_std, any_of("treat.ref"))
  }
}

# Reference Data   ############################

refdat_basic <- data.frame(
  rowname    = c(paste0("es_A", LETTERS[1:NT]), "tau", "seed"),
  graphname  = c(paste0("treat_name", LETTERS[1:NT]), "tau", "seed"),
  graphname2 = c(paste0("treat_name", LETTERS[1:NT]), "tau", "seed"),
  Estimate   = c(es_MD, tau, seedval),
  Comparator = c(LETTERS[1:NT], NA, NA),
  Reference  = c(rep("A", NT), NA, NA),
  paramtype  = c(rep("es", NT), "tau", "seed")
)

refdat <- refdat_basic

# Add all pairwise comparisons
for (i in 1:NT) {
  shifted <- refdat_basic %>%
    group_by(paramtype) %>%
    mutate(
      Reference = lag(Comparator, i),
      rowname   = paste(paramtype, Reference, Comparator, sep = "_"),
      Estimate  = ifelse(paramtype == "es", Estimate - lag(Estimate, i), Estimate),
      graphname2 = case_when(
        paramtype == "es" ~ paste0("treat_name", Comparator),
        TRUE ~ graphname
      )
    ) %>%
    filter(!is.na(Reference)) %>%
    ungroup()
  refdat <- bind_rows(refdat, shifted)
  gc()
}

# Reverse comparisons
refdat <- refdat %>%
  mutate(
    rowname = ifelse(!paramtype %in% c("tau","seed"),
                     paste(paramtype, Comparator, Reference, sep = "_"), rowname),
    Estimate = ifelse(paramtype == "es", -Estimate, Estimate),
    Comparatornew = Reference,
    Reference = Comparator,
    graphname = ifelse(paramtype == "es", paste0("treat_name", Comparatornew), graphname),
    graphname2 = ifelse(paramtype == "es",
                        paste0("treat_name", Comparatornew),
                        graphname2)
  ) %>%
  select(-Comparator) %>%
  rename(Comparator = Comparatornew) %>%
  bind_rows(refdat) %>%
  distinct() %>%
  mutate(comp = paste(Comparator, Reference))

# Add interaction
refdat <- as.data.frame.table(beta_x) %>%
  rename(Estimate = Freq, Reference = Var1, Comparator = Var2) %>%
  mutate(
    rowname    = paste(Comparator, Reference),
    graphname  = paste0("`", Comparator, ":xcov`"),
    graphname2 = paste0("`", Comparator, "_vs_", Reference, ":Ivalrand:xcov`"),
    paramtype  = "beta_x",
    comp       = paste(Comparator, Reference)
  ) %>%
  bind_rows(refdat)

