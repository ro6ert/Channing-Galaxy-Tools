echo "starting $1"
R --vanilla --slave --args $1 < ldracecompare.R > ld_$1.rlog

