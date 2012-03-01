#!/bin/bash
#for RACE in ASW CEU CHB CHD GIH JPT LWK MEX MKK TSI YRI
for RACE in ASW CEU CHD GIH LWK MEX MKK YRI
do
FNAME=hapmap3_r3_b36_fwd.${RACE}.qc.poly
echo "doing qsub makeibs.sh" $FNAME
qsub makeibs.sh $FNAME
done

