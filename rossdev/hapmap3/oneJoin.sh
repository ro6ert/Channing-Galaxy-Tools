#!/bin/bash
# tricky script - need to write a merge-list file for plink
# since it won't take multiple merge file names on the cl
# make autosomal then add X for all
# ross lazarus april 2010 for the WTCCC data
cd $1
OUTF="joinme.files"
echo "" > $OUTF 
# create joinfile
for C in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22
do
   N="WTCCC_Affx_$1_$C"
   echo "$N.bed $N.bim $N.fam" >> $OUTF
   # write a line
done
plink --make-bed --out WTCCC_Affx_$1_AUTO --merge-list $OUTF
# now add X for all file
N="WTCCC_Affx_$1_X"
plink --bfile WTCCC_Affx_$1_AUTO --make-bed --out WTCCC_Affx_$1_ALL --merge $N.bim $N.bed $N.fam
cd ..
