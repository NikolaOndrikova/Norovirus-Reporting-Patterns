#setwd("")

# read in NUTS1 regional data
eng_map <- readOGR(dsn = 'data/Eng_Map/Eng_Map.shp', layer = "Eng_Map")

# read in synthetic PHE norovirus data
synth_noro_data <- read_csv("data/Synth_Noro_Data.csv")

### REGIONAL FRACTIONS ---------------------------------------------------------
# regional fraction of total number of primary schools - 2016
schools_ = c(0.05156, 0.14585, 0.10651, 0.09739, 0.10585, 0.11867,  0.10806, 0.15496, 0.11116)

# regional fraction of total number of nurseries - 2016
nurseries_ = c(0.08374, 0.17500, 0.06897, 0.07635, 0.15025, 0.10099, 0.19704, 0.10591, 0.04433)

# regional fraction of total number of hospitals - 2015
hospitals_ = c(0.05747, 0.11823, 0.10016, 0.06732, 0.10673, 0.10673, 0.13465, 0.16831, 0.14039)

regional_predictors_v = list("schools" = schools_,
                             "nurseries" = nurseries_,
                             "hospitals" = hospitals_)

# mid-2016 population estimates
pop_rg <- c(148432, 442876, 331343, 278396, 365266, 379042, 635561, 542018, 306112,#[00-04]
            292049, 845266, 638444, 542371, 702896, 727238, 1067276, 1079043, 603937,#[05-14]
            337720, 896996, 708247, 602269, 739628, 690408, 1037771, 1071558, 654769,#[15-24]
            639932, 1832398, 1362022, 1158204, 1468968, 1554861, 3082454, 2263111, 1291537,#[25-44]
            711109, 1880120, 1396233, 1247487, 1462775, 1594585, 1943706, 2360129, 1468663,#[45-64]
            507606, 1321967, 989452, 895710, 1061201, 1184408, 1021124, 1710438, 1190935) #[65 + ]
pop_rg_df <- data.frame(t(pop_rg))
rg_colnames <- c('UKC.[00-04]', 'UKD.[00-04]', 'UKE.[00-04]', 'UKF.[00-04]', 'UKG.[00-04]', 'UKH.[00-04]', 'UKI.[00-04]', 'UKJ.[00-04]', 'UKK.[00-04]',
                 'UKC.[05-14]', 'UKD.[05-14]', 'UKE.[05-14]',	'UKF.[05-14]', 'UKG.[05-14]',	'UKH.[05-14]', 'UKI.[05-14]',	'UKJ.[05-14]', 'UKK.[05-14]',
                 'UKC.[15-24]',	'UKD.[15-24]', 'UKE.[15-24]',	'UKF.[15-24]', 'UKG.[15-24]',	'UKH.[15-24]', 'UKI.[15-24]',	'UKJ.[15-24]', 'UKK.[15-24]',
                 'UKC.[25-44]',	'UKD.[25-44]', 'UKE.[25-44]',	'UKF.[25-44]', 'UKG.[25-44]',	'UKH.[25-44]', 'UKI.[25-44]',	'UKJ.[25-44]', 'UKK.[25-44]',
                 'UKC.[45-64]',	'UKD.[45-64]', 'UKE.[45-64]',	'UKF.[45-64]', 'UKG.[45-64]',	'UKH.[45-64]', 'UKI.[45-64]',	'UKJ.[45-64]', 'UKK.[45-64]',
                 'UKC.[65+', 'UKD.[65+', 'UKE.[65+', 'UKF.[65+', 'UKG.[65+', 'UKH.[65+', 'UKI.[65+', 'UKJ.[65+', 'UKK.[65+')
colnames(pop_rg_df) <- rg_colnames

### POLYMOD CONTACT DATA -------------------------------------------------------
data(polymod)
uk_cm = contact_matrix(polymod, countries = "United Kingdom", 
                       #filter = list("phys_contact"=1),
                       age.limits = c(0, 5, 15, 25, 45, 65))
uk_contacts = uk_cm$matrix
attr(uk_contacts, "agedistri") <- c(0.09397, 0.20178, 0.16222, 0.24827, 0.23838, 0.05539)
rm(polymod)