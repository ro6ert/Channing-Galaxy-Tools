#!/bin/bash
# start with ped/map to pbed
#plink --file $1 --make-bed --out $1     
# select ld independent markers
#plink --bfile $1 --out $1 --indep-pairwise 50 30 0.2     
#plink --bfile $1 --extract $1.prune.in --make-bed --out $1_INDEP     
#plink --bfile $1_INDEP --thin 0.1 --out $1_INDEP_THIN --make-bed
# take a 10% sample and run the ibs stats
#plink --bfile $1_INDEP_THIN --genome --out $1_INDEP_THIN
# for running grr
#python /share/shared/galaxy/tools/rgenetics/rgGRR.py "./$1" "$1" "./$1_GRR.html" "." "$1_GRR" '5000' '6'
# estimate pairwise r2 - problem, must remove same number of founders from all or ld is inflated

OPTS="--bfile $1  --r2 --filter-founders --ld-window-kb 100 --ld-window-r2 0.1"
# setup plink to estimate pairwise r^2 over as far as 2MB and not truncate values
# should we only analyse above some threshold?
plink --chr 2 --from-kb 95000 --to-kb 97000 --out $1_chr2all --remove $1_remove_unrelated.id $OPTS
plink --chr 2 --from-kb 95000 --to-kb 97000 --out $1_chr2unrel --remove $1_remove_related.id $OPTS
plink --chr 4 --from-kb 123000 --to-kb 125000 --out $1_chr4all --remove $1_remove_unrelated.id $OPTS 
plink --chr 4 --from-kb 123000 --to-kb 125000 --out $1_chr4unrel --remove $1_remove_related.id $OPTS 
plink --chr 6 --from-kb 99000 --to-kb 101000 --out $1_chr6all --remove $1_remove_unrelated.id $OPTS 
plink --chr 6 --from-kb 99000 --to-kb 101000 --out $1_chr6unrel --remove $1_remove_related.id $OPTS 
plink --chr 8 --from-kb 69000 --to-kb 71000 --out $1_chr8all --remove $1_remove_unrelated.id $OPTS 
plink --chr 8 --from-kb 69000 --to-kb 71000 --out $1_chr8unrel --remove $1_remove_related.id $OPTS 
plink --chr 10 --from-kb 19000 --to-kb 21000 --out $1_chr10all --remove $1_remove_unrelated.id $OPTS 
plink --chr 10 --from-kb 19000 --to-kb 21000 --out $1_chr10unrel --remove $1_remove_related.id $OPTS 

