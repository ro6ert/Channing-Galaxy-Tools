#!/bin/sh
# for WTCCC ross lazarus
# april 2010
# requires sge 
# run this after running fixWTCCC.sh
# for each collection, join chromosome pbeds 
for D in 58C_Chiamo BD CAD CD NBS RA T1D T2D 
 do
  qsub oneJoin.sh $D
  # make autosome and all
 done
