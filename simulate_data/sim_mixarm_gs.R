########################
trts <- seq_len(NT)

# All possible two-arm studies
#
combs2 <- t(combn(trts, 2))
#
dat.combs2 <- data.frame(combs2, row.names = seq_len(nrow(combs2))) %>%
  rename(treat1 = X1, treat2 = X2)

# All possible three-arm studies
#
combs3 <- t(combn(trts, 3))
#
dat.combs3 <- data.frame(combs3, row.names = seq_len(nrow(combs3))) %>%
  rename(treat1 = X1, treat2 = X2, treat3 = X3)

dat.studies <- vector("list", N.sim)
# creates structure of network ############
for (i in 1) { 
    n.multi<-round(ctmixarm(NT, NS)$nmulti)
    #
    j <- 1
    k <- 0 # stop after 3000 repetitions
    #
    while (j <= n.multi) {
      k <- k + 1
      #
      trts.ij <- dat.combs3[sample(seq_len(nrow(dat.combs3)), 1), ]
      combs2.ij <- t(combn(as.vector(t(trts.ij)), 2))
      #
      dat.ij <- data.frame(combs2.ij, row.names = seq_len(nrow(combs2.ij))) %>%
        rename(treat1 = X1, treat2 = X2) %>%
        mutate(study = j) %>%
        relocate(study, .before = treat1) %>%
        mutate(comparison = paste(treat1, "vs", treat2))
      #
      dat.studies.ij <- rbind(dat.studies, dat.ij)
      #
      if (!any(table(dat.studies.ij$comparison) > 3)) {
        dat.studies <- dat.studies.ij
        j <- j + 1
      }
      #
      if (j > 3000)
        stop("Max reached (j = 3000)")
    }
    #
    combs_multi <- dat.studies %>%
      select(treat1, treat2) %>% unique() %>%
      mutate(multi = "yes")
    #
    # Combine all possible pairwise comparisons with created three-arm studies
    #
    dat.merge <- merge(dat.combs2, combs_multi,
                       by = c("treat1", "treat2"), all.x = TRUE, all.y = TRUE)
    #
    # Add pairwise comparisons only consisting of two-arm studies
    #
    dat.2arm <- dat.merge %>% filter(is.na(multi)) %>% select(-multi)
    dat.2arm <- do.call(rbind, replicate(NS, dat.2arm, simplify = FALSE))
    
    dat.2arm %<>% mutate(study = row_number() + n.multi) %>%
      relocate(study, .before = treat1) %>%
      mutate(comparison = paste(treat1, "vs", treat2))
    #
    dat.studies <- rbind(dat.studies, dat.2arm)
    #
    # Add two-arm studies to pairwise comparisons with three-arm studies
    #
    dat.3arm <- dat.merge %>% filter(!is.na(multi)) %>% select(-multi) %>%
      mutate(comparison = paste(treat1, "vs", treat2)) %>%
      unique() # unique() necessary?
    #
    ttab <- table(dat.studies$comparison)
    #
    for (j in dat.3arm$comparison) {
      for (k in seq_len(NS - ttab[j])) {
        dat.jk <- dat.3arm %>% filter(comparison == j) %>%
          mutate(study = max(dat.studies$study) + 1)
        #
        dat.studies <- bind_rows(dat.studies, dat.jk)
      }
    }
  
  
}

dat.studies

dat.studies_long<-
  dat.studies%>%
  select(study, treat=treat1, comparison)%>%
  bind_rows(dat.studies%>%
              select(study, treat=treat2, comparison))%>%
  add_count(study)%>% # arms per study
  mutate(comparison2=ifelse(n>2, "multiarm", comparison))%>%
  select(-c(n, comparison))%>% # discard old comparison info and temp count variable
  distinct()%>%
  arrange(comparison2, study) %>% # creates ids for joins, it needs to start at 1 for each type of comparison
  group_by(comparison2, study) %>%
  mutate(rank = cur_group_id())%>%
  ungroup()%>%
  group_by(comparison2)%>%
  mutate(rank=rank-min(rank)+1)%>%
  ungroup()#%>%mutate(events=10, n=100) # dummy events to check network geometry


# cleanup
rm(combs_multi, combs2, combs2.ij, combs3, dat.2arm, dat.3arm, dat.combs2, dat.combs3, dat.ij, dat.jk, 
   dat.merge, dat.studies, dat.studies.ij, i, j, k, n.multi, trts,
   trts.ij,  ttab)