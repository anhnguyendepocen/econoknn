source("../hierarchical-estimation/logspec/tableapi.R", chdir=T)

save.demean <- function(outfile, df, ycol, xcol, factorouts) {
    df$.one. <- 1
    ta.save.demeaned(outfile, df, ycol, c(xcol, '.one.'), factorouts)
}

