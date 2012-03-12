rgeo=function(geoID,tdir)
{
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
# note we need to convert GSE's differently from GDS geo datasets
# and for the moment, let's ignore that nasty complexity
res = ''
outf_eSet = file.path(tdir,paste(geoID,'.eset',sep='')) # prototype output files for creation
outf_maList = file.path(tdir,paste(geoID,'.malist',sep=''))
outf_pheno = file.path(tdir,paste(geoID,'.pheno',sep=''))
ngs = '### makenewgalaxy' # string in log to be parsed to create new history entries post exec
cat('# geoquery.R, part of the Rexpression Galaxy toolkit http://esphealth.org/trac/rgalaxy\n')
cat(paste('# Starting at ',date(),'\n'))
cat(paste('# geoquery.R got',geoID,'as the GEO id and',outf_eSet,'+',outf_maList,'as outfiles\n'))
library(GEOquery)
geo_prefix = substring(geoID,1,3) # eg GDS
s = paste('### geo prefix=',geo_prefix)
res = paste(res,s,sep='\n')
if (geo_prefix == 'GDS')
# GEO DataSets (GDSxxx) are curated sets of GEO Sample data. A GDS record represents
# a collection of biologically and statistically comparable GEO Samples and forms the 
# basis of GEO's suite of data display and analysis tools. Samples within a GDS refer 
# to the same Platform, that is, they share a common set of probe elements. 
# Value measurements for each Sample within a GDS are assumed to be calculated 
# in an equivalent manner, that is, considerations such as background processing 
# and normalization are consistent across the dataset. Information reflecting 
# experimental design is provided through GDS subsets.
{
 gds = getGEO(geoID) 
 platrecgs = Meta(gds)$platform
 platrec = getGEO(platrecgs)
 maList = GDS2MA(gds,GPL=platrec) # make an MAList (named maList when loaded!)
 save(maList,file=outf_maList) # supplied by Galaxy for the user's history
 res = paste(res,'maList saved')
 res = paste(res,s,sep='\n')
 eSet = GDS2eSet(gds) # make an ExpressionSet (named eSet when loaded!)
 phe = pData(eSet) # get pheno
 save(eSet,file=outf_eSet) # no matter what the filename
 write.table(phe,file=outf_pheno,sep='\t',quote=F) 
 s = ('eSet saved')
 res = paste(res,s,sep='\n')
 s = paste('# Done at ',date(),'\n')
 res = paste(res,s,sep='\n')
 s = paste(ngs,'\t',outf_eSet,'\t','eset','\t',outf_pheno,'\n',sep='')
 res = paste(res,s,sep='\n')
 s = paste(ngs,'\t',outf_maList,'\t','malist','\t',outf_pheno,'\n\n',sep='')
 res = paste(res,s,sep='\n')
} # GDS

if (geo_prefix == 'GSE')
# The GSE is the most confusing of the GEO entities. A GSE entry can represent an 
# arbitrary number of samples run on an arbitrary number of platforms. The GSE has 
# a metadata section, just like the other classes. However, it doesn’t have a 
# GEODataTable. Instead, it contains two lists, accessible using GPLList and GSMList, 
# that are each lists of GPL and GSM objects. 
{
 # we could get the list of esets as esetlist = getGEO(gs,GSEMatrix=T)
 # but for now let's bail out
 s = 'Sorry, this version does not handle GEO Series data sets'
 res = paste(res,s,sep='\n')
} # GSE

if (geo_prefix == 'GSM')
{
 # we could get the list of esets as esetlist = getGEO(gs,GSEMatrix=T)
 # but for now let's bail out
 s = 'Sorry, this version does not handle GEO Series data sets'
 res = paste(res,s,sep='\n') 
} # GSM
if (geo_prefix == 'GPL')
{
 s = 'Sorry, this version does not handle GEO platform data sets'
 res = paste(res,s,sep='\n')
} # GPL
return(res)
}




