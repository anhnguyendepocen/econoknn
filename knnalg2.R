source("../hierarchical-estimation/logspec/tableapi.R", chdir=T)

save.demean <- function(outfile, df, ycol, xcol, factorouts) {
    df$.one. <- 1
    ta.save.demeaned(outfile, df, ycol, c(xcol, '.one.'), factorouts)
}

all.smooth <- function(df, xvar, yvars, groupvars) {
    groups <- ""
    for (groupvar in groupvars)
        groups <- paste(groups, df[, groupvar])

    for (group in unique(groups)) {
        subdf <- df[groups == group,]
        for (yvar in yvars)
            subdf[groups == group, yvar] <- loess(as.formula(paste(yvar, "~", xvar)), subdf)
    }

    df
}
