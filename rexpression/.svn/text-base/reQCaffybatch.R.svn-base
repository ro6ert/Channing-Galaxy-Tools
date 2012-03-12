# first try at a BioC R script for array QC
# copyright 2008, ross lazarus 
# started August 22 2008
# licensed under the LGPL http://www.gnu.org/licenses/lgpl.html
# Not fit for human consumption


nparm = 3
appname = 'reQCaffybatch.R'
cl = commandArgs(trailingOnly=T)
if (length(cl) < nparm) 
{
 cat(paste(length(cl),' arguments supplied, ',nparm,' needed'))
 cat("Call as:  R CMD BATCH '--vanilla --slave --args inaffybatch outhtml outdir ' reQCaffybatch.R logf \n")
 q(save="no",status=1)
} # short command line error
ineset = cl[1]
outf_html = cl[2]
outdir = cl[3]
# debug
cat('# ',appname,'.R, part of the Rexpression Galaxy toolkit http://esphealth.org/trac/rgalaxy\n')
cat(paste('# Starting at ',date(),'\n'))
cat(paste('# Got parameters',cl,'\n'))
library(affyQCReport)
ae = load(ineset) # ae contains the name so use get(ae)
# note this assumes we load an affybatch but we have no real
# idea of it's name in the environment - so we store it in ae
# there has to be a a better way to do this in R.
par(mfrow=c(1,1)) # page these out as png - pdf is hundreds of MB!
impng = paste('image%d.',ae,'.png')
png(filename=impng) # should give sensible name
par(ask=F)
image(get(ae))
dev.off()
par(mfrow=c(3,3)) # or whatever works
mapdf = paste('MAplot.',ae,'.pdf')
pdf(mapdf,h=10,w=8)
MAplot(get(ae), cex=0.25, cex.main=0.5,plot.method="smoothScatter")
dev.off()
# QCReport(ETABM25.affybatch)
q1pdf = paste('qcTitle.',ae,'.pdf')
pdf(q1pdf,h=10,w=8)
titlePage(get(ae))
dev.off()
q2pdf = paste('qcSignal.',ae,'.pdf')
pdf(q2pdf,h=10,w=8)
signalDist(ETABM25.affybatch)
dev.off()
q3pdf = paste('qcMain.',ae,'.pdf')
pdf(q3pdf,h=10,w=8)
plot(qc(get(ae)))
dev.off()
q4pdf = paste('qcBorder1.',ae,'.pdf')
pdf(q4pdf,h=10,w=8)
borderQC1(get(ae))
dev.off()
q5pdf = paste('qcBorder2.',ae,'.pdf')
pdf(q5pdf,h=10,w=8)
borderQC2(get(ae))
dev.off()
q6pdf = paste('qcCorrplot.',ae,'.pdf')
pdf(q6pdf,h=10,w=8)
correlationPlot(get(ae))
dev.off()