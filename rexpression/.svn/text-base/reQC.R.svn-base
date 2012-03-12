# first try at a BioC R script
# copyright 2008, ross lazarus 
# started August 6 2008
# licensed under the LGPL http://www.gnu.org/licenses/lgpl.html
# Not fit for human consumption
#
# Notes and code mostly cribbed from the BioC GEOquery vignette
#
# this is a non-interactive script intended as the
# executable for a Galaxy tool
# to create a new ExpressionSet and MAList from a GEO GDSxxxx id
# using the GEOquery package from BioC
#
# finally figured out how to finesse R's blathering output to stderr
# which makes the current Galaxy job runner think the job has failed
#
# Call as:  R CMD BATCH '--vanilla --slave --args GDS10 eSet maList ' geoquery.R logf
# both out paths will be supplied by Galaxy
# geoID will be supplied by user as a parameter
#

nparm = 4
appname = 'ArrayExpress'
cl = commandArgs(trailingOnly=T)
if (length(cl) < nparm) 
{
 cat(paste(length(cl),' arguments supplied, ',nparm,' needed'))
 cat("Call as:  R CMD BATCH '--vanilla --slave --args AEid affyBatch, eSet maList ' arrayexpress.R logf \n")
 q(save="no",status=1)
} # short command line error
AEId = cl[1]
outf_affybatch = cl[2]
outf_eSet = cl[3]
outf_maList = cl[4]
# debug
cat('# ',appname,'.R, part of the Rexpression Galaxy toolkit http://esphealth.org/trac/rgalaxy\n')
cat(paste('# Starting at ',date(),'\n'))
cat(paste('# Got parameters',AEId,'as the AE id, and ',affybatch,'+', outf_eSet,'+',outf_maList,'as outfiles\n'))
library(ArrayExpress)
cat(paste('### Array Express id =',AEId,'\n'))
ae = ArrayExpress(AEId) # eg "E-TABM-25")
klass = class(ae)[1] # eg "AffyBatch"
save(affybatch,file=outf_affybatch)

