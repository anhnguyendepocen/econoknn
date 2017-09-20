setwd("~/research/gcp/hierarchical-estimation/jingyuan")

library(readstata13)

redo.demean <- F

if (redo.demean) {
    source("../logspec/tableapi.R", chdir=T)

    agegrp          <- 3
    ## 0: all three groups together
    ## 1: 0-4
    ## 2: 5-64
    ## 3: 65+

    ## import file
    ##df <- read.dta13(paste0("ToMLE/mortality_AEA_agegrp", agegrp, ".dta"))
    df <- read.dta13("~/Dropbox/MyGCP/empirics/mortality_AEA_time-variant-inc.dta")
    df <- subset(df, agegrp == 33)


    factorouts <- c('precip1_GMFD : factor(adm1)', 'precip2_GMFD : factor(adm1)', 'factor(adm2)', 'factor(adm1) : factor(year)')

    ta.save.demeaned(paste0("knntest-", agegrp, "-new.RData"), df, 'deathrate',
                     c('GMFD_poly1', 'precip1_GMFD'), factorouts)
} else {
    df <- read.dta13("~/Dropbox/MyGCP/empirics/mortality_AEA_time-variant-inc.dta")
    df <- subset(df, agegrp == 33)
}

load("knntest-3-new.RData")

library(Hmisc)

df.income <- read.dta13("~/Dropbox/MyGCP/empirics/gdppc_admin1_8countries.dta")
df$income <- df$gdppc # This might just be constant at country level

for (idch in unique(df.income$id_chicago)) {
    if (is.na(idch))
        next
    ## Generate income model
    subdf.income <- subset(df.income, id_chicago == idch)
    subdf <- subset(df, id_chicago == idch)

    incomes <- exp(approxExtrap(subdf.income$year, log(subdf.income$gdppcstate_hat_loglinear), subdf$year)$y)
    df$income[!is.na(df$id_chicago) & df$id_chicago == idch] <- incomes
}

## KNN Analysis

df$income.rank <- rank(df$income)
df$climtas.rank <- rank(df$Tmean_GMFD)
df$temp.rank <- rank(dmxxs[,1])

## Analysis for K

loo.predict <- function(k) {
    df$predicted <- NA
    for (ii in sample(1:nrow(df), 1000)) {
        dists <- (df$income.rank[ii] - df$income.rank)^2 + (df$climtas.rank[ii] - df$climtas.rank)^2 + (df$temp.rank[ii] - df$temp.rank)^2
        dists[ii] <- Inf
        df$predicted[ii] <- mean(dmyy[order(dists)[1:k]])
    }

    sqrt(mean((df$deathrate - df$predicted)^2, na.rm=T))
}

results <- data.frame(kk=c(), rmse=c())
for (plus in c(100, 300, 1000)) { #0, 1, 3, 10, 30)) { #
    for (kk in c(800) + plus) { # c(1, 3, 10, 30, 100, 300, 1000, 500, 700)
        print(kk)
        results <- rbind(results, data.frame(kk, rmse=loo.predict(kk)))
    }
}

library(ggplot2)

ggplot(results, aes(kk, rmse)) +
    geom_point() + geom_smooth() +
    scale_y_log10() + scale_x_log10()

## KNN curve analysis

KK <- 750

get.knn <- function(income.rank, climtas.rank, temp.rank) {
    dists <- (income.rank - df$income.rank)^2 + (climtas.rank - df$climtas.rank)^2 + (temp.rank - df$temp.rank)^2
    mean(dmyy[order(dists)[1:KK]])
}

income.ranks <- quantile(df$income.rank, c(.25, .5, .75))
climtas.ranks <- quantile(df$climtas.rank, c(.25, .5, .75))
income.values <- quantile(df$income, c(.25, .5, .75), na.rm=T)
climtas.values <- quantile(df$Tmean_GMFD, c(.25, .5, .75))

results2 <- data.frame(tas=c(), deathrate=c(), income=c(), climtas=c())
for (zz1 in 1:3) {
    for (zz2 in 1:3) {
        subres <- data.frame(tas=c(), deathrate=c(), income=c(), climtas=c())
        for (tas in seq(0, 40, length.out=21)) {
            print(c(zz1, zz2, tas))
            dists <- abs(tas - dmxxs[,1])
            temp.rank <- mean(df$temp.rank[which(dists == min(dists))])
            deathrate <- get.knn(income.ranks[zz1], climtas.ranks[zz2], temp.rank)
            subres <- rbind(subres, data.frame(tas, deathrate, income=round(income.values[zz1], -1), climtas=round(climtas.values[zz2], 1)))
        }
        subres$deathrate <- subres$deathrate - subres$deathrate[subres$tas == 20]
        results2 <- rbind(results2, subres)
    }
}

results2$income <- factor(results2$income, rev(unique(results2$income)))

ggplot(results2, aes(tas, deathrate)) +
    facet_grid(income ~ climtas) +
    xlab("Temperature") + ylab("Death Rate") +
    geom_smooth() + scale_x_continuous(expand=c(0, 0)) + theme_minimal()
ggsave("knn-nonant.pdf", width=7, height=5)

## Who is in each?

income.quants <- quantile(df$income, c(0, 1/3, 2/3, 1), na.rm=T)
climtas.quants <- quantile(df$Tmean_GMFD, c(0, 1/3, 2/3, 1))
income.values <- round(quantile(df$income, c(.25, .5, .75), na.rm=T), -1)
climtas.values <- round(quantile(df$Tmean_GMFD, c(.25, .5, .75)), 1)

df$income.quant <- NA
for (zz in 1:3)
    df$income.quant[df$income >= income.quants[zz] & df$income < income.quants[zz+1]] <- income.values[zz]

df$climtas.quant <- NA
for (zz in 1:3)
    df$climtas.quant[df$Tmean_GMFD >= climtas.quants[zz] & df$Tmean_GMFD < climtas.quants[zz+1]] <- climtas.values[zz]

library(dplyr)

df$myiso <- df$iso
df$myiso[nchar(df$iso) == 2] <- "EUR"

df2 <- group_by_(df[!is.na(df$income.quant) & !is.na(df$climtas.quant),], .dots = c('income.quant', 'climtas.quant', 'myiso')) %>%
    summarize(counts = n()) %>%
    mutate(perc = counts / sum(counts)) %>%
    arrange(desc(perc)) %>%
    mutate(label_pos = cumsum(perc) - perc / 2,
           perc_text = paste0(round(perc * 100), "%"))

df2$income.quant <- factor(df2$income.quant, unique(df2$income.quant))

ggplot(df2, aes(x="", y=perc, fill=myiso)) +
    facet_grid(income.quant ~ climtas.quant) +

    ## make stacked bar chart with black border
    geom_bar(stat = "identity", color = "black", width = 1) +

    ## convert to polar coordinates
    coord_polar(theta = "y") +

    ## formatting
    scale_y_continuous(breaks = NULL) +
    scale_fill_discrete(name = "") +
    theme(text = element_text(size = 22),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()) +
    theme_minimal() + xlab(NULL) + ylab(NULL)
ggsave("nonant-countries.pdf", width=8, height=6)

## Limit to the US

df.save <- df
df <- subset(df, iso == "USA")

results2 <- data.frame(tas=c(), deathrate=c(), income=c(), climtas=c())
for (zz1 in 1:3) {
    for (zz2 in 1:3) {
        subres <- data.frame(tas=c(), deathrate=c(), income=c(), climtas=c())
        for (tas in seq(0, 40, length.out=21)) {
            print(c(zz1, zz2, tas))
            dists <- abs(tas - dmxxs[,1])
            temp.rank <- mean(df$temp.rank[which(dists == min(dists))])
            deathrate <- get.knn(income.ranks[zz1], climtas.ranks[zz2], temp.rank)
            subres <- rbind(subres, data.frame(tas, deathrate, income=round(income.values[zz1], -1), climtas=round(climtas.values[zz2], 1)))
        }
        subres$deathrate <- subres$deathrate - subres$deathrate[subres$tas == 20]
        results2 <- rbind(results2, subres)
    }
}

results2$income <- factor(results2$income, rev(unique(results2$income)))

ggplot(results2, aes(tas, deathrate)) +
    facet_grid(income ~ climtas) +
    xlab("Temperature") + ylab("Death Rate") +
    geom_smooth() + scale_x_continuous(expand=c(0, 0)) + theme_minimal()
ggsave("knn-nonant-usa.pdf", width=7, height=5)

df <- df.save

## Limit to no-US

df.save <- df
df <- subset(df, iso != "USA")

results2 <- data.frame(tas=c(), deathrate=c(), income=c(), climtas=c())
for (zz1 in 1:3) {
    for (zz2 in 1:3) {
        subres <- data.frame(tas=c(), deathrate=c(), income=c(), climtas=c())
        for (tas in seq(0, 40, length.out=21)) {
            print(c(zz1, zz2, tas))
            dists <- abs(tas - dmxxs[,1])
            temp.rank <- mean(df$temp.rank[which(dists == min(dists))])
            deathrate <- get.knn(income.ranks[zz1], climtas.ranks[zz2], temp.rank)
            subres <- rbind(subres, data.frame(tas, deathrate, income=round(income.values[zz1], -1), climtas=round(climtas.values[zz2], 1)))
        }
        subres$deathrate <- subres$deathrate - subres$deathrate[subres$tas == 20]
        results2 <- rbind(results2, subres)
    }
}

results2$income <- factor(results2$income, rev(unique(results2$income)))

ggplot(results2, aes(tas, deathrate)) +
    facet_grid(income ~ climtas) +
    xlab("Temperature") + ylab("Death Rate") +
    geom_smooth() + scale_x_continuous(expand=c(0, 0)) + theme_minimal()
ggsave("knn-nonant-nousa.pdf", width=7, height=5)

df <- df.save

## Just US by decade

income.quants <- quantile(df$income, c(0, 1/3, 2/3, 1), na.rm=T)
climtas.quants <- quantile(df$Tmean_GMFD, c(0, 1/3, 2/3, 1))
income.values <- round(quantile(df$income, c(.25, .5, .75), na.rm=T), -1)
climtas.values <- round(quantile(df$Tmean_GMFD, c(.25, .5, .75)), 1)

df$income.quant <- NA
for (zz in 1:3)
    df$income.quant[df$income >= income.quants[zz] & df$income < income.quants[zz+1]] <- income.values[zz]

df$climtas.quant <- NA
for (zz in 1:3)
    df$climtas.quant[df$Tmean_GMFD >= climtas.quants[zz] & df$Tmean_GMFD < climtas.quants[zz+1]] <- climtas.values[zz]

library(dplyr)

df.save <- df
df <- subset(df, iso == "USA")

df$decade <- factor(floor(df$year / 10) * 10)

df2 <- group_by_(df[!is.na(df$income.quant) & !is.na(df$climtas.quant),], .dots = c('income.quant', 'climtas.quant', 'decade')) %>%
    summarize(counts = n()) %>%
    mutate(perc = counts / sum(counts)) %>%
    arrange(desc(perc)) %>%
    mutate(label_pos = cumsum(perc) - perc / 2,
           perc_text = paste0(round(perc * 100), "%"))

df2$income.quant <- factor(df2$income.quant, rev(unique(df2$income.quant)))

ggplot(df2, aes(x="", y=perc, fill=decade)) +
    facet_grid(income.quant ~ climtas.quant) +

    ## make stacked bar chart with black border
    geom_bar(stat = "identity", color = "black", width = 1) +

    ## convert to polar coordinates
    coord_polar(theta = "y") +

    ## formatting
    scale_y_continuous(breaks = NULL) +
    scale_fill_discrete(name = "") +
    theme(text = element_text(size = 22),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()) +
    theme_minimal() + xlab(NULL) + ylab(NULL)
ggsave("nonant-us-decades.pdf", width=7, height=4)

## Just dropping post-2000 USA

df.save <- df
df <- subset(df, iso != "USA" | year < 2000)

results2 <- data.frame(tas=c(), deathrate=c(), income=c(), climtas=c())
for (zz1 in 1:3) {
    for (zz2 in 1:3) {
        subres <- data.frame(tas=c(), deathrate=c(), income=c(), climtas=c())
        for (tas in seq(0, 40, length.out=21)) {
            print(c(zz1, zz2, tas))
            dists <- abs(tas - dmxxs[,1])
            temp.rank <- mean(df$temp.rank[which(dists == min(dists))])
            deathrate <- get.knn(income.ranks[zz1], climtas.ranks[zz2], temp.rank)
            subres <- rbind(subres, data.frame(tas, deathrate, income=round(income.values[zz1], -1), climtas=round(climtas.values[zz2], 1)))
        }
        subres$deathrate <- subres$deathrate - subres$deathrate[subres$tas == 20]
        results2 <- rbind(results2, subres)
    }
}

results2$income <- factor(results2$income, rev(unique(results2$income)))

ggplot(results2, aes(tas, deathrate)) +
    facet_grid(income ~ climtas) +
    xlab("Temperature") + ylab("Death Rate") +
    geom_smooth() + scale_x_continuous(expand=c(0, 0)) + theme_minimal()
ggsave("knn-nonant-nousa2000.pdf", width=7, height=5)

df <- df.save

## Look at how 35 C effect changes with income

dists <- abs(35 - dmxxs[,1])
temp.rank.35c <- mean(df$temp.rank[which(dists == min(dists, na.rm=T))])

dists <- abs(17 - dmxxs[,1])
temp.rank.17c <- mean(df$temp.rank[which(dists == min(dists, na.rm=T))])

get.knn.2 <- function(income.rank, temp.rank) {
    dists <- (income.rank - df$income.rank)^2 + (temp.rank - df$temp.rank)^2
    mean(dmyy[order(dists)[1:KK]])
}

subres <- data.frame(deathrate=c(), income=c())
for (income in exp(seq(log(min(df$income, na.rm=T)), log(max(df$income, na.rm=T)), length.out=50)[2:49])) {
    dists <- abs(income - df$income)
    income.rank <- mean(df$income.rank[which(dists == min(dists, na.rm=T))])
    print(c(income, income.rank))
    deathrate.35c <- get.knn.2(income.rank, temp.rank.35c)
    deathrate.17c <- get.knn.2(income.rank, temp.rank.17c)
    subres <- rbind(subres, data.frame(deathrate=deathrate.35c - deathrate.17c, income=income))
}

ggplot(subres, aes(income, deathrate)) +
    geom_smooth() + scale_x_log10()
