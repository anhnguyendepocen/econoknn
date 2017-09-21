setwd("~/research/gcp/econoknn")

library(readstata13)
library(Hmisc)
source("knnalg2.R")

redo.demean <- F

df <- read.dta13("/home/solomon/Dropbox/GCP_Reanalysis/interpolation/data/consolidated/mortality_AEA_time-variant-inc_2.dta")
df <- subset(df, agegrp == 33)
df <- subset(df, !is.na(precip1_GMFD) & !is.na(precip2_GMFD))

if (redo.demean) {
    factorouts <- c('precip1_GMFD : factor(adm1)', 'precip2_GMFD : factor(adm1)', 'factor(adm2)', 'factor(adm1) : factor(year)')

    save.demean("time-variant.RData", df, 'deathrate', 'GMFD_poly1', factorouts)
}

load("time-variant.RData")

df$dmyy <- dmyy
df$dmxx <- dmxxs[,1]
df$income.rank <- rank(df$gdppcstate)
df$climtas.rank <- rank(df$Tmeanstate_GMFD)
df$temp.rank <- rank(df$GMFD_poly1)

KK <- 750

get.knn.beta <- function(income.rank, climtas.rank, temp.rank) {
    dists <- (income.rank - df$income.rank)^2 + (climtas.rank - df$climtas.rank)^2 + (temp.rank - df$temp.rank)^2
    mod <- lm(dmyy ~ 0 + dmxx, df[order(dists)[1:KK],])
    c(coef(mod), vcov(mod))
}

get.knn.curve <- function(income.rank, climtas.rank, tas0, tas1, length.out) {
    result <- data.frame(tas=tas0, deathrate=0, var=0)
    for (tas in seq(tas0, tas1, length.out=length.out)[-1]) {
        dists <- abs(tas * 365 - df$GMFD_poly1)
        temp.rank <- mean(df$temp.rank[which(dists == min(dists))])
        betavar <- get.knn.beta(income.rank, climtas.rank, temp.rank)

        deathrate <- betavar[1] * (tas - tas0)
        var <- betavar[2] * (tas - tas0)^2
        result <- rbind(result, data.frame(tas, deathrate, var))
    }

    result
}

income.ranks <- quantile(df$income.rank, c(.25, .5, .75))
climtas.ranks <- quantile(df$climtas.rank, c(.25, .5, .75))
income.values <- quantile(df$gdppcstate, c(.25, .5, .75), na.rm=T)
climtas.values <- quantile(df$Tmeanstate_GMFD, c(.25, .5, .75))

results <- data.frame(tas=c(), deathrate=c(), var=c(), income=c(), climtas=c())
for (zz1 in 1:3) {
   for (zz2 in 1:3) {
        right <- get.knn.curve(income.ranks[zz1], climtas.ranks[zz2], 20, 40, 11)
        right$income <- round(income.values[zz1], -1)
        right$climtas <- round(climtas.values[zz2], 1)

        left <- get.knn.curve(income.ranks[zz1], climtas.ranks[zz2], 20, 0, 11)
        left$income <- round(income.values[zz1], -1)
        left$climtas <- round(climtas.values[zz2], 1)

        results <- rbind(results, right, left[-1,])
    }
}

results$income <- factor(results$income, rev(sort(unique(results$income))))

ggplot(results, aes(tas, deathrate)) +
    facet_grid(income ~ climtas) +
    xlab("Temperature") + ylab("Death Rate") +
    geom_smooth() + scale_x_continuous(expand=c(0, 0)) + theme_minimal()
ggsave("knn-nonant.pdf", width=7, height=5)

## Split out by country

df$myiso <- df$iso
df$myiso[nchar(df$iso) == 2] <- "EUR"

for (myiso in unique(df$myiso)) {
    subdf <- df[df$myiso == myiso,]

    get.knn.beta <- function(income.rank, climtas.rank, temp.rank) {
        dists <- (income.rank - subdf$income.rank)^2 + (climtas.rank - subdf$climtas.rank)^2 + (temp.rank - subdf$temp.rank)^2
        mod <- lm(dmyy ~ 0 + dmxx, subdf[order(dists)[1:KK],])
        c(coef(mod), vcov(mod))
    }

    get.knn.curve <- function(income.rank, climtas.rank, tas0, tas1, length.out) {
        result <- data.frame(tas=tas0, deathrate=0, var=0)
        for (tas in seq(tas0, tas1, length.out=length.out)[-1]) {
            dists <- abs(tas * 365 - subdf$GMFD_poly1)
            temp.rank <- mean(subdf$temp.rank[which(dists == min(dists))])
            betavar <- get.knn.beta(income.rank, climtas.rank, temp.rank)

            deathrate <- betavar[1] * (tas - tas0)
            var <- betavar[2] * (tas - tas0)^2
            result <- rbind(result, data.frame(tas, deathrate, var))
        }

        result
    }

    results <- data.frame(tas=c(), deathrate=c(), var=c(), income=c(), climtas=c())
    for (zz1 in 1:3) {
        for (zz2 in 1:3) {
            right <- get.knn.curve(income.ranks[zz1], climtas.ranks[zz2], 20, 40, 11)
            right$income <- round(income.values[zz1], -1)
            right$climtas <- round(climtas.values[zz2], 1)

            left <- get.knn.curve(income.ranks[zz1], climtas.ranks[zz2], 20, 0, 11)
            left$income <- round(income.values[zz1], -1)
            left$climtas <- round(climtas.values[zz2], 1)

            results <- rbind(results, right, left[-1,])
        }
    }

    results$income <- factor(results$income, rev(sort(unique(results$income))))

    ggplot(results, aes(tas, deathrate)) +
        facet_grid(income ~ climtas) +
        xlab("Temperature") + ylab("Death Rate") +
        geom_smooth() + scale_x_continuous(expand=c(0, 0)) + theme_minimal()
    ggsave(paste0("knn-nonant-", myiso, ".pdf"), width=7, height=5)

}

## Uninteracted

get.knn.beta <- function(temp.rank) {
    dists <- (temp.rank - df$temp.rank)^2
    mod <- lm(dmyy ~ 0 + dmxx, df[order(dists)[1:KK],])
    c(coef(mod), vcov(mod))
}

get.knn.curve <- function(tas0, tas1, length.out) {
    result <- data.frame(tas=tas0, deathrate=0, var=0)
    for (tas in seq(tas0, tas1, length.out=length.out)[-1]) {
        dists <- abs(tas * 365 - df$GMFD_poly1)
        temp.rank <- mean(df$temp.rank[which(dists == min(dists))])
        betavar <- get.knn.beta(temp.rank)

        deathrate <- betavar[1] * (tas - tas0)
        var <- betavar[2] * (tas - tas0)^2
        result <- rbind(result, data.frame(tas, deathrate, var))
    }

    result
}

right <- get.knn.curve(20, 40, 11)
left <- get.knn.curve(20, 0, 11)
results <- rbind(right, left[-1,])

ggplot(results, aes(tas, deathrate)) +
    xlab("Temperature") + ylab("Death Rate") +
    geom_smooth() + scale_x_continuous(expand=c(0, 0)) + theme_minimal()
ggsave(paste0("knn-uninteracted.pdf"), width=7, height=5)
