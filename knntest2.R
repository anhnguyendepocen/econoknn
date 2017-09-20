setwd("~/research/gcp/econoknn")

source("knnalg2.R")

library(readstata13)

df <- read.dta13("/home/solomon/Dropbox/GCP_Reanalysis/interpolation/data/consolidated/mortality_AEA_time-variant-inc_2.dta")
df <- subset(df, agegrp == 33)
df <- subset(df, !is.na(precip1_GMFD) & !is.na(precip2_GMFD))

factorouts <- c('precip1_GMFD : factor(adm1)', 'precip2_GMFD : factor(adm1)', 'factor(adm2)', 'factor(adm1) : factor(year)')

save.demean("time-variant.RData", df, 'deathrate', 'GMFD_poly1', factorouts)
