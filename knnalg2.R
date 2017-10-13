source("../hierarchical-estimation/logspec/tableapi.R", chdir=T)

save.demean <- function(outfile, df, ycol, xcol, factorouts) {
    df$.one. <- 1
    ta.save.demeaned(outfile, df, ycol, c(xcol, '.one.'), factorouts)
}

all.smooth <- function(df, xvar, yvars, groupvars, span=.3) {
    groups <- ""
    for (groupvar in groupvars)
        groups <- paste(groups, df[, groupvar])

    for (group in unique(groups)) {
        subdf <- df[groups == group,]
        for (yvar in yvars)
            df[groups == group, yvar] <- lowess(subdf[, xvar], subdf[, yvar], f=span)$y
    }

    df
}
