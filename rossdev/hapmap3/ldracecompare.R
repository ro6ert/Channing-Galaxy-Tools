# compare related and unrelated ld values from a sample
# both calculated from the same number of founders to avoid bias in ld est
# ross lazarus april 17 2010
library(ggplot2)
race=commandArgs(trailingOnly=T)[[1]]
chroms = c(2,4,6,8,10)
relmeans = c()
unrelmeans = c()
relmedians = c()
unrelmedians =  c()
relsds = c()
unrelsds = c()
deltamean = c()
deltamedian = c()
for (chrom in chroms) 
{
unf=paste('hapmap3_r3_b36_fwd.',race,'.qc.poly_chr',chrom,'unrel.ld',sep='')
relf=paste('hapmap3_r3_b36_fwd.',race,'.qc.poly_chr',chrom,'all.ld',sep='')
pngf=paste('hapmap3_',race,'_ld',chrom,'.png',sep='')
unrel = read.table(unf,head=T)
rel = read.table(relf,head=T)
medr = median(rel$R2)
meanr = mean(rel$R2)
relmedians = append(relmedians,medr)
relmeans = append(relmeans,meanr)
medu = median(unrel$R2)
meanu = mean(unrel$R2)
unrelmedians = append(unrelmedians,medu)
unrelmeans = append(unrelmeans,meanu)
deltamedian = append(deltamedian,(medr - medu))
deltamean = append(deltamean,(meanr - meanu))
relsds = append(relsds,sd(rel$R2))
unrelsds = append(unrelsds,sd(unrel$R2))
k = ks.test(unrel$R2,rel$R2)
print(k)
w = wilcox.test(unrel$R2,rel$R2)
print(w)
names(unrel)[names(unrel)=="R2"] = "UnrelR2"
names(rel)[names(rel)=="R2"] = "AllR2"
both=merge(rel,unrel,by=c("SNP_A","SNP_B","CHR_A","CHR_B","BP_A","BP_B"))
q = qplot(UnrelR2,AllR2,data=both,shape=1,main=paste('HapMap',race,'CHR',chrom,'related cf unrelated LD (rsq)'))
ggsave(filename=pngf)
}
umed = mean(unrelmedians)
um = mean(unrelmeans)
usd = mean(unrelsds)
rmed = mean(relmedians)
rm = mean(relmeans)
rsd = mean(relsds)
print(paste('For race=',race,'meanrelmedian - meanunrelmedian r^2=',rmed-umed))
print(paste('mean delta median=',mean(deltamedian),' = ',100*mean(deltamedian)/rmed,'%',sep=''))
print(paste('Unrel r^2 mean=',um,'sd=',usd,'median',umed))
print(paste('Rel r^2 mean=',rm,'median=',rmed,'sd=',rsd))
desc=data.frame(chr=chroms,relm=relmeans,relmed=relmedians,relsd=relsds,unrelm=unrelmeans,unrelsd=unrelsds,unrelmed=unrelmedians,deltam=deltam)
desc
title=paste(RelMean/UnrelMean for',race)
pdf(filename=pngf,height=6,width=6)
plot(relmeans,relsds,col=race,main=title,legend=("topright","Race"))
dev.off()
print('done')



