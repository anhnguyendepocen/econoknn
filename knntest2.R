setwd("~/research/gcp/econoknn")

source("knnalg2.R")

library(readstata13)

df <- read.dta13("/home/solomon/Dropbox/GCP_Reanalysis/interpolation/data/consolidated/mortality_AEA_time-variant-inc_2.dta")
df <- subset(df, agegrp == 33)

factorouts <- c('precip1_GMFD : factor(adm1)', 'precip2_GMFD : factor(adm1)', 'factor(adm2)', 'factor(adm1) : factor(year)')

save.demean(outfile, df, 'deathrate', 'GMFD_poly1', factorouts)
