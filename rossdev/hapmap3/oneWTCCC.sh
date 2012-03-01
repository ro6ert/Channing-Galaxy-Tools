#!/bin/sh
# for WTCCC
# ross lazarus april 2010
# convert a single tped to a pbed
cd $1
plink --tped Affx_gt_$1_Chiamo_$2.tped --make-bed --out WTCCC_Affx_$1_$2 --tfam Affx_samples_$1_Chiamo.tfam 
cd ..
