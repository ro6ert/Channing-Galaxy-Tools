#!/bin/sh
# convert WTCCC data collections from tped to pbed
# use joinWTCCC.sh to create auto and all pbed files
# for each collection, convert each chromosome
# note needed to copy the irregularly named 58C tfam file to get this to work
for D in 58C_Chiamo BD CAD CD NBS RA T1D T2D 
 do
  for C in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 X
  do
   ## oneWTCCC.sh contains plink --tfile $D --make-bed --out $D
   qsub oneWTCCC.sh $D $C
   echo "doing $D $C"
  done
 done
