# generate annual net returns for forest

#load packages
library(RSQLite)
library(reshape2)
library(stargazer)
library(ggplot2)
library(tidyverse)

#load data
# sol <- read.csv('Z:/forest_project/faustmann/faustmann_result.csv')
# sol <- sol[c('fips','spgrpcd','fortypcd','optage')]
# sol$fips <- formatC(sol$fips, format='d', width = 5, flag='0')
bert <- tbl_df(readRDS('growth_parameters.rds')) %>%
  mutate(fips = substring(id, 1, 5),
         spgrpcd = substring(id, 6, 7),
         fortypcd = substring(id, 8, 10))


obs.age <- tbl_df(readRDS("D:/GroupWork/FIA Data/Harvest Age/obs_harvest_age.rds")) %>%
  select(fips, spgrpcd, fortypcd, age.mean)



#reference tables
nbr <- readRDS("county_neighbor_mapping.rds")
georef <- readRDS("fia_regions_mapping.rds")
species <- readRDS("fia_spgrpcd_majspgrpcd_mapping.rds")




###############################################################################
# clean price data

prices <- read.csv('price_data_v2.csv')
prices$id <- formatC(prices$id, format='d',width=7,flag='0')
prices$fips <- substring(prices$id,1,5)
prices$spgrpcd <- as.numeric(substring(prices$id,6,7))
prices <- prices[-c(1:2,20)]
prices <- melt(prices,id.vars = c('fips','spgrpcd'),variable.name='year',value.name = 'stumpage_price')
prices$year <- as.numeric(substring(prices$year,2,5))
prices$stumpage_price <- prices$stumpage_price / (25/3)

###############################################################################
###############################################################################
# 58,626 vonberalanffy functions
# 48,665 observed rotation ages
# extrapolate over space,species,forest type for 9,959 missing rotation ages

dat$id <- paste0(dat$fips, formatC(dat$spgrpcd, format = 'd', width = 2, flag = '0'), formatC(dat$fortypcd, format='d', width =3, flag='0'))
dat <- merge(dat, georef, by='fips')

missing <- dat[is.na(dat$age.mean),]

# length of missing ages is 9,959
# replace missing age with regional average for that species and forest type

for (i in 1:length(unique(missing$id))){
  county <- substring(missing$id[i],1,5)
  spgrp <- as.integer(substring(missing$id[i],6,7))
  fortyp <- as.integer(substring(missing$id[i],8,10))
  region <- missing$subregion[i]
  sub <- dat[dat$spgrpcd==spgrp & dat$fortypcd==fortyp & dat$subregion==region,]
  replace.age <- mean(sub$age.mean, na.rm = T)
  dat$age.mean[dat$spgrpcd==spgrp & dat$fortypcd==fortyp & dat$fips==county] <- replace.age
}

nans <- dat[is.nan(dat$age.mean),]

# length of missing ages in now 5,731
# replace missing age with regional average for that species in any forest type

for (i in 1:length(unique(nans$id))){
  county <- substring(nans$id[i],1,5)
  spgrp <- as.integer(substring(nans$id[i],6,7))
  fortyp <- as.integer(substring(nans$id[i],8,10))
  region <- nans$subregion[i]
  sub <- dat[dat$spgrpcd==spgrp & dat$subregion==region,]
  replace.age <- mean(sub$age.mean, na.rm = T)
  dat$age.mean[dat$spgrpcd==spgrp & dat$fortypcd==fortyp & dat$fips==county] <- replace.age
}

nans <- dat[is.nan(dat$age.mean),]

# length of missing ages is now 21
# replace missing age with species average in any forest type and region

for (i in 1:length(unique(nans$id))){
  county <- substring(nans$id[i],1,5)
  spgrp <- as.integer(substring(nans$id[i],6,7))
  fortyp <- as.integer(substring(nans$id[i],8,10))
  region <- nans$subregion[i]
  sub <- dat[dat$spgrpcd==spgrp,]
  replace.age <- mean(sub$age.mean, na.rm = T)
  dat$age.mean[dat$spgrpcd==spgrp & dat$fortypcd==fortyp & dat$fips==county] <- replace.age
}


# that completes the extrapolation. all missing ages are recovered.

##############################################################################

dat <- merge(dat,prices,by=c('fips','spgrpcd'),all.x=T)

# Faustmann net return
#calculate annualized net return per cubic foot by species group using faustmann formula

#observed rotation age


dat$nrtree_faust_obsage <- ((dat$stumpage_price*(dat$a_par * (1 - exp(-dat$b_par*dat$age.mean))^3)-dat$cost) / (exp(0.05*dat$age.mean) - 1)) * 0.05


#faustmann optimized rotation age

dat$nrtree_faust_optage <- ((dat$stumpage_price*(dat$a_par * (1 - exp(-dat$b_par*dat$optage))^3)-dat$cost) / (exp(0.05*dat$optage) - 1)) * 0.05


##############################################################################

# single rotation value
# calculate annualized net return per cubic foot by species group
nr <- dat

# convert time from years to days for discount factor calculation

nr$vol <- nr$a_par * (1 - exp(-nr$b_par*nr$age.mean))^3
nr$revenue <- nr$stumpage_price*nr$vol-nr$cost

for (i in 1:length(unique(nr$id))) {
  dfact <- 0
  for (j in 1:round(nr$age.mean[i], digits=0)) {
    dfact <- dfact + exp(-.05*j)
  }
  nr$ann_nr_single[i] <- (nr$revenue[i] * exp(-0.05*j)) / dfact
  nr$nr_single[i] <- (nr$revenue[i] * exp(-.05*j))
}

###############################################################################

a <- readRDS("D:/GroupWork/FIA Data/fortyp_acres_by_county.rds")
t <- readRDS("D:/GroupWork/FIA Data/fortyp_treecount_by_county.rds")
v <- readRDS("D:/GroupWork/FIA Data/fortyp_volume_by_county.rds")

d <- merge(t,v,by=c('STATECD','COUNTYCD','SPGRPCD','FORTYPCD'))
d <- merge(d,a,by=c('STATECD','COUNTYCD','FORTYPCD'))
d$fips <- paste(formatC(d$STATECD,format = 'd',width = 2,flag = '0'),
                formatC(d$COUNTYCD,format = 'd',width = 3,flag = '0'),sep='')
d <- rename(d, c('FORTYPCD'='fortypcd','SPGRPCD'='spgrpcd'))

#merge tree acreage data with net returns

d <- merge(d[c('fips','spgrpcd','fortypcd','trees','acres','volume')],
           nr[c('fips','spgrpcd','fortypcd','year','ann_nr_single','nrtree_faust_optage','nrtree_faust_obsage')],
           by=c('fips','spgrpcd','fortypcd'))

# drop observations with 0 acres

d <- d[d$acres!=0,]

# trees per acre

d$tpa <- d$trees / d$acres

# scale to per acre measure

d$nracre <- d$ann_nr_single * d$tpa
d$nracre_faust_optage <- d$nrtree_faust_optage * d$tpa
d$nracre_faust_obsage <- d$nrtree_faust_obsage * d$tpa

#remove outliers: 1% off each tail

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.01, .99), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

d$nracre <- remove_outliers(d$nracre)
d <- subset(d, !is.nan(nracre) & !is.na(nracre))

###############################################################################

#create composite measure of forestland net returns
#for all species, by forest type, and by species-forest pairs

nr_forest <- ddply(d, .(fips,year), summarize, nracre=weighted.mean(nracre, w=volume, na.rm = T),
                   nracre_faust_optage = weighted.mean(nracre_faust_optage, w = volume, na.rm = T),
                   nracre_faust_obsage = weighted.mean(nracre_faust_obsage, w = volume, na.rm = T))

nr_forest <- subset(nr_forest, !is.nan(nracre) & !is.na(nracre))

###############################################################################

saveRDS(nr_forest, 'frnr_by_year.rds')
# save copy for logit landuse est data
saveRDS(nr_forest, "Z:/land_use_change_estimation_data/nr_forest_faust.rds")
