#for RACE in ASW CEU CHB CHD GIH JPT LWK MEX MKK TSI YRI
for RACE in ASW CEU CHD GIH LWK MEX MKK YRI
#for RACE in CEU
do
echo "doing qsub ldcompare.sh" $RACE
qsub ldcompare.sh $RACE
done
