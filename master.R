# master file

# read in source files
source('./src/required_packages.R')
source('./src/custom_functions.R')
source('./src/import_data.R')
source('./src/hhh4_data_prep.R')

### CREATE THE STS OBJECT FOR HHH4 ---------------------------------------------

start_w <- 27
start_Y <- 2014 

# split the data
TRAIN <- 2:(4*52)
TEST <- max(TRAIN) + 1:52

# make hhh4 data object
nov_sts <- sts(nov_observed, 
               start = c(as.numeric(start_Y), as.numeric(start_w)), 
               frequency = 52,
               population = pop_m, 
               neighbourhood = eng_nbOrder_expand) 

shared_data <- list(t= epoch(nov_sts) - 1,
                    pop = population(nov_sts)/rowSums(population(nov_sts)))

# helper variables
REGIONS <- unique(stratum(nov_sts, 1))
NREGIONS <- length(REGIONS)
GROUPS <- unique(stratum(nov_sts, 2))
NGROUPS <- length(GROUPS)

# setup a model matrix with group indicators
group_indicators <- sapply(GROUPS, function (g) {
  index <- which(stratum(nov_sts, which = 2) == g)
  res <- col(nov_sts)
  res[] <- res %in% index
  res
}, simplify = FALSE, USE.NAMES = TRUE)

# fix names of the age groups as these contain numbers
qGROUPS <- paste0("`", GROUPS, "`")

# names for formulas specifications
region_specific = c("log(schools)", "log(hospitals)", "log(nurseries)") 

seasonality = c("sin(2 * pi * t/52)", "cos(2 * pi * t/52)",
                "sin(4 * pi * t/52)", "cos(4 * pi * t/52)")


### SPECIFY AND FIT THE FIXED EFFECTS OPTIONS ----------------------------------

## endemic formula 1: ~group + log(schools_) + log(hospitals_) + "log(nurseries_)" + 2:(sin+cos)
end_intrcpt_group <- reformulate(c(1,qGROUPS[-1], region_specific, seasonality),
                                 intercept = TRUE)

## epidemic formula 1: ~group + log(pop) 
epi_intrcpt_group <- reformulate(c(1,qGROUPS[-1], "log(pop)"),
                                 intercept = TRUE)

## epidemic formula 2: ~group + log(pop) + 1:(sin+cos)
epi_intrcpt_groupS1 <- reformulate(c(1,qGROUPS[-1],"log(pop)", seasonality),
                                   intercept = TRUE)

# baseline fixed effects fits for the mixed effects models 
fixed_models = list("epi_noSeason" = hhh4(nov_sts, 
                                          makehhh4Control(end_option = end_intrcpt_group,
                                                          epi_option = epi_intrcpt_group)),
                    "epi_season1" = hhh4(nov_sts, 
                                         makehhh4Control(end_option = end_intrcpt_group, 
                                                         epi_option = epi_intrcpt_groupS1)))

gc()


### SPECIFY and FIT THE RANDOM EFFECTS OPTIONS ---------------------------------

combination_grid <- expand.grid(fixed = c("end_intrcpt_group","epi_intrcpt_group", "epi_intrcpt_groupS1"), 
                                random = c("rand_intercept", "rand_corr_intercepts"))

row.names(combination_grid) <- do.call("paste", c(combination_grid, list(sep = "|")))

random_opt_test <- apply(X = combination_grid, MARGIN = 1, FUN = function (options) {
  updateOpt <- function (option) switch(option,
                                        "rand_intercept" = "ri(type = 'iid') - 1",
                                        "rand_corr_intercepts" = "ri(type = 'iid', corr = 'all') - 1",
                                        "end_intrcpt_group" = end_intrcpt_group, 
                                        "epi_intrcpt_group" = epi_intrcpt_group, 
                                        "epi_intrcpt_groupS1" = epi_intrcpt_groupS1)
  reformulate(paste(updateOpt(options[1]), 
                    updateOpt(options[2]), sep="+")[2])
})

ri_models = list("epi_noSeason" = updatehhh4(fixed_models$epi_noSeason,
                                             end_option = random_opt_test$`end_intrcpt_group|rand_intercept`, 
                                             epi_option = random_opt_test$`epi_intrcpt_group|rand_intercept`), 
                 "epi_season1" = updatehhh4(fixed_models$epi_season1,
                                            end_option = random_opt_test$`end_intrcpt_group|rand_intercept`, 
                                            epi_option = random_opt_test$`epi_intrcpt_groupS1|rand_intercept`),
                 "epi_noSeason_corr" = updatehhh4(fixed_models$epi_noSeason, 
                                                  end_option = random_opt_test$`end_intrcpt_group|rand_corr_intercept`, 
                                                  epi_option = random_opt_test$`epi_intrcpt_group|rand_corr_intercept`), 
                 "epi_season1_corr" = updatehhh4(fixed_models$epi_season1,
                                                 end_option = random_opt_test$`end_intrcpt_group|rand_corr_intercepts`, 
                                                 epi_option = random_opt_test$`epi_intrcpt_groupS1|rand_corr_intercepts`))
gc()


### SUMMARIES OF ALL THE MODELS ------------------------------------------------

print(lapply(fixed_models, summary, idx2Exp = TRUE, amplitudeShift = TRUE, maxEV = TRUE))

print(lapply(ri_models, summary, idx2Exp = TRUE, amplitudeShift = TRUE, maxEV = TRUE))


### VALIDATION AND COMPARISON OF ALL THE MODELS --------------------------------

set.seed(1245)

# following lines can take ~ 35 minutes or more, depends on a computer
owas0 <- lapply(fixed_models, oneStepAhead, tp = range(TEST)-1, type = "rolling")
owas1 <- lapply(ri_models, oneStepAhead, tp = range(TEST)-1, type = "rolling")

# all of the One Week Ahead predictionS (=OWAS)
owas = list("A1:epi_noSeason" = owas0$epi_noSeason,
            "A2:epi_season1" = owas0$epi_season1,
            "B1:ri_epi_noSeason" = owas1$epi_noSeason,
            "B2:ri_epi_season1" = owas1$epi_season1,
            "C1:ri_epi_noSeason_corr" = owas1$epi_noSeason_corr,
            "C2:ri_epi_season1_corr" = owas1$epi_season1_corr)

best_models = selectBestModel(owas, metrics = c("rps", "logs"), 
                              verbose = TRUE, plot = FALSE)

final_mod_name = substr(best_models$logs, 7, nchar(best_models$logs))
final_model <- ri_models[[final_mod_name]] 

### PLOT THE BEST FIT ----------------------------------------------------------

## by age group
plotHHH4_fitted_groups(final_model, end = c(2018, 26),
                       groups = stratum(nov_sts, which = 2), units = NULL, pch = 20, legend = 2,
                       col = c('#cccccc','#4ba173','#fff176'),
                       names = c('<4 yrs', '5-14 yrs', '15-24 yrs',
                                 '25-44 yrs', '45-64 yrs', '65+ yrs'),
                       legend.args = list(cex = 0.8,
                                          legend = c("from other age groups", 
                                                     "within age group", 
                                                     "endemic")))

## by district
plotHHH4_fitted_groups(final_model, end = c(2018, 26),
                       groups = factor(stratum(nov_sts, which = 1), levels = REGIONS),
                       names = eng_map$nuts118nm,
                       units = NULL,
                       col = c('#cccccc','#4ba173','#fff176'),
                       legend = 4, legend.args = list(cex = 0.8,
                                                      legend = c("from other regions", 
                                                                 "within region", 
                                                                 "endemic")))


### MAPS FOR RANDOM INTERCEPTS -------------------------------------------------

# comp options are "end.ri(iid)" or "epi.ri(iid)"
comp = "end.ri(iid)"
ranef_matrix = exp(ranef(final_model, tomatrix = TRUE))

map <- as(eng_map, "SpatialPolygonsDataFrame")
map$`<4 yrs` <- ranef_matrix[1:9, comp]
map$`5-14 yrs` <- ranef_matrix[10:18, comp]
map$`15-24 yrs` <- ranef_matrix[19:27, comp]
map$`25-44 yrs` <- ranef_matrix[28:36, comp]
map$`45-64 yrs` <- ranef_matrix[37:45, comp]
map$`65+ yrs` <- ranef_matrix[46:54, comp]
map_sf = st_as_sf(map)

map_sf <- map_sf %>% 
  select(`<4 yrs`, `5-14 yrs`, `15-24 yrs`,`25-44 yrs`, `45-64 yrs`, `65+ yrs`, geometry) %>% 
  gather(VAR, RI, -geometry)

ggplot() + 
  geom_sf(data = map_sf, aes(fill=RI, size = 0.15)) +
  scale_size_identity() +
  facet_wrap(~VAR, ncol = 3) +
  scale_fill_viridis(name = "RI", option = "D", alpha = 0.70, 
                     direction = 1, begin = 0, end = 0.95) +
  geom_text(data=map@data,aes(x=long, y=lat,label = nuts118cd), size = 2) +
  theme_void(base_size = 14) +
  theme(legend.title = element_text(face = "bold.italic",
                                    size = 13),
        title = element_text(face = "bold",
                             size = 13))
