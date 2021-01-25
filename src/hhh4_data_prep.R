# prepare data for hhh4 data with age
week_age_region_data = aggregate(Count_CDR_OPIE_ID ~ Week_Key + PHE_Region_Cd + nw_age_grp, data=synth_noro_data, sum)
colnames(week_age_region_data) <- c("Yw", "Region", "AgeGrp", "Count")

# combine age group and region columns
week_age_region_data$Region.AgeGroup <- apply(week_age_region_data[, c("AgeGrp", "Region")],
                                              1, paste, collapse="")
# remove unnecessary columns
week_age_region_data <- week_age_region_data[, !(names(week_age_region_data) %in% c("Region","AgeGrp"))]

# generate spatial weights
eng_adjmap <- poly2adjmat(eng_map)
if (requireNamespace("spdep")) {
  eng_nbOrder <- nbOrder(eng_adjmap, maxlag=Inf)
}

# we can remove the eng_adjmap
rm(eng_adjmap)

## 1. Create initial data for sts object for hhh4 model -------------------------------------------

# set the time frame
date_from <- 201426 # ~"2014-07-06"
date_to <-  201927 #~"2019-06-30"

# split the data
#TRAIN <- 2:(4*52)
#TEST <- max(TRAIN) + 1:52


hhh4_data <- week_age_region_data[week_age_region_data$Yw > date_from & week_age_region_data$Yw < date_to, ]
nov_observed = data.matrix(as.data.frame.matrix(xtabs( Count ~ Yw + Region.AgeGroup, data=hhh4_data)))


### 2. Getting the right matrix shapes, expanding and normalising ----------------------------------------

# make population matrices for hhh4
pop_m <- data.matrix(pop_rg_df[rep(seq_len(nrow(pop_rg_df)), each=nrow(nov_observed)),])

# fix col names to match
colnames(nov_observed) <- colnames(pop_m)

# row-normalized contact matrix
uk_con_norm <- uk_contacts / rowSums(uk_contacts)

# spatial weights
eng_nbOrder_expand0 = matrix(rep(t(eng_nbOrder), 6), ncol=ncol(eng_nbOrder), byrow = TRUE)
eng_nbOrder_expand = matrix(rep(eng_nbOrder_expand0, 6), nrow = nrow(eng_nbOrder_expand0))
colnames(eng_nbOrder_expand) <- colnames(nov_observed)

rm(eng_nbOrder_expand0)


# list of region-specific predictors
regional_predictors = lapply(regional_predictors_v, expandRegionalPreditor, 
                             my_Matrix = nov_observed,
                             n = 6)
