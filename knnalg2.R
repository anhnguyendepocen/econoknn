source("../hierarchical-estimation/logspec/tableapi.R", chdir=T)

save.demean <- function(outfile, df, ycol, xcol, factorouts) {
    df$.one. <- 1
    ta.save.demeaned(outfile, df, ycol, c(xcol, '.one.'), factorouts)
}

ma <- function(x, n=5) {
   filter(x, rep(1/n, n), sides=2)
}

all.smooth <- function(df, xvar, yvars, groupvars, span=.3) {
    groups <- ""
    for (groupvar in groupvars)
        groups <- paste(groups, df[, groupvar])

    for (group in unique(groups)) {
        subdf <- df[groups == group,]
        for (yvar in yvars)
	    df[groups == group, yvar] <- ma(subdf[, yvar], ceiling(nrow(subdf) * span))
            ## df[groups == group, yvar] <- lowess(subdf[, xvar], subdf[, yvar], f=span)$y
    }

    df
}
