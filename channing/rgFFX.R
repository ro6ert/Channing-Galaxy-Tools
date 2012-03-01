# x2wig.R
# demonstrate galaxy interface to GGtools
# parameters: package holding ffX2 resources, probeID, ffchrind, snplocchr
args = commandArgs(TRUE)
if (length(args) != 5) stop("need args pkg probeID ffchrind snplocind outfmt")
pkg = args[1]
probeID = args[2]
ffchrind = as.numeric(args[3])
snplocchr = args[4]
outfmt = args[5]
library(pkg, character.only=TRUE)
mgr = setupFF()  # must be in pkg
library(snpLocs130)
snplocref = data(list=paste("snpLocs_", snplocchr, sep=""))
library(GGtools)
gt = geneTrack(mgr, probeID, ffchrind, snpLocs_20)
wgt = width(gt)
gt = gt[ width(gt) == 1, ]
export(gt, paste("gt",outfmt,sep="."))

