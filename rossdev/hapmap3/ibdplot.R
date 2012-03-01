library(ggplot2)

races = c('ASW','CEU','CHB', 'CHD','GIH', 'JPT', 'LWK', 'MEX', 'MKK', 'TSI','YRI') 
for (race in races) 
{
f=paste('hapmap3_r3_b36_fwd.',race,'.qc.poly_INDEP_THIN.genome',sep='')
pdff=paste('hapmap3_',race,'_Z0_Z1.pdf',sep='')
d = read.table(f,head=T)
#names(d)
#[1] "FID1"   "IID1"   "FID2"   "IID2"   "RT"     "EZ"     "Z0"     "Z1"
#[9] "Z2"     "PI_HAT" "PHE"    "DST"    "PPC"    "RATIO"
q = qplot(Z0,Z1,data=d,colour=RT,shape=1,main=paste('HapMap',race,'Plink IBD'))
ggsave(filename=pdff)
}
