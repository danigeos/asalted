#!/bin/csh
#This script uses GMT commands for graphic output
set file_original = DEMO.asalted.result
set file = $file_original.tmp
set ps = DEMO.asalted.ps

set z0=0

#rm .gmtdefaults4

cat <<END> hypso0.xa.tmp
-11035	0
-5500	0.3e14
-3798  	2.1e14
-200	3.4e14
$z0 	3.607162e14
+100	3.7e14
END
cat <<END> hypso1.xa.tmp
-2911 	0
-2509	0.15e12
-1500	0.25e12
-1000	0.31e12
-600	0.37e12
-430	0.39e12
-200	0.50e12
$z0   	0.54e12
+100	0.6e12
END
cat <<END> hypso2.xa.tmp
-3824 	0
-2718 	0.47e12
-1500	0.91e12
-1000	1.07e12
-800	1.14e12
-600	1.22e12
-430  	1.39e12
-200  	1.56e12
$z0	1.94e12
+100	2.0e12
END


#2011 PAPER reference model with accelerated uplift to force 17 cycles:
asalted  -s0/1000 -m.014/.002/-284/-30 -u0.0048/.00003 -k1e-6/1.5/50 -W -V3 -r0/4.5e3/12e3 -p0/.6 -e0/1.1/1.3 -w6/0.5/m \
	-M0/2 -Shypso0.xa.tmp -B-100/-1000/100000/hypso1.xa.tmp -b-430/-1000/100000/hypso2.xa.tmp \
	-t0/120000/.2 -z0/-1/-2 > $file_original







echo _______________asalted.gmt.csh_______________
#set ref_time = `awk '{if (NR>4) {if (NR==5) z2ini=$15; if ($15>z2ini+15) {print $1; exit;}}}' $file_original
#echo; echo "ref_time at " $ref_time " yr"

set timeini = `head $file_original | awk '{if (substr($0,1,1)!="#") {print $1; exit}}'`
set timeend = `tail $file_original | awk '{if (substr($0,1,1)!="#") time=$1}END{print time}'`
echo Time: $timeini $timeend

cat $file_original \
	| acortar 1 10 \
	> $file

set timeend_ky = `tail $file | awk '{if ($1>tmax && substr($0,1,1)!="#") tmax=$1}END{print tmax/1000}'`
set totalflow = `tail $file | awk '{if ($12>tf && substr($0,1,1)!="#") tf=$12}END{print tf}'`
echo timeend=$timeend_ky "kyr;    totalflow="$totalflow km3
awk 	'{if (NR>4) {incirate=$8; if (incimax*1.<incirate*1.) {incimax=incirate;time_incimax=$1}} \
		time_prev=$1; sill_prev==$2; \
	} END {printf("\nMax. incision rate:    %.2f m/yr at t= %.1f kyr", incimax, time_incimax/1000);}' $file 
awk 	'{if (NR>4) {var=$5; if (max*1.<var*1.) {max=var;time_max=$1}} \
		time_prev=$1; var_prev==$2; \
	} END {printf("\nMax. flow velocity:    %.2f m/s at t= %.1f kyr", max, time_max/1000);}' $file 
awk 	'{if (NR>4) {var=$6; if (max*1.<var*1.) {max=var;time_max=$1}} \
		time_prev=$1; var_prev==$2; \
	} END {printf("\nMax. discharge:        %.2e m3/s at t= %.1f kyr", max, time_max/1000);}' $file 
awk 	'{if (NR>4 && substr($0,1,1)!="#") {rate=($14-var_prev)/($1-time_prev); if (ratemax*1.<rate*1.) {ratemax=rate; time_max=$1}} \
		time_prev=$1; var_prev=$14; \
	} END {printf("\nMax. water level rate: %.2f m/yr at t= %.1f kyr", ratemax, time_max/1000);}' $file 
echo
#tail -2 $file 



set timelabelint = `awk -v timeend_ky=$timeend_ky 'BEGIN{if (timeend_ky<20) print 1; else print 10;}' `
set timetickint = `echo $timelabelint \/ 10 | bc `
#echo $timelabelint $timetickint

gmtset BASEMAP_FRAME_RGB +0/0/0 
psbasemap -JX14.5/4 -R$timeini/$timeend_ky/0/100 \
	-Ba$timelabelint\f$timetickint\:"time t [kyr]":/a50f10:"gypsum rate [10@+12@+ kg yr@+-1@+]":NsW \
	-X3 -Y23 -K -P > $ps
#GYPSUM:
awk '(substr($0,1,1)!="#"){s=$17; diff=s-saltant; if (diff==0 || diffant==0) print ">"; else print $1/1000,diff/($1-tant); saltant=$17; tant=$1; diffant=diff}' $file | psxy -JX -R -M -W2   -O -K >> $ps
awk '(substr($0,1,1)!="#"){s=$18; diff=s-saltant; if (diff==0 || diffant==0) print ">"; else print $1/1000,diff/($1-tant); saltant=$18; tant=$1; diffant=diff}' $file | psxy -JX -R -M -W2ta/30/60/0   -O -K >> $ps
gmtset BASEMAP_FRAME_RGB +200/0/0 
psbasemap -JX -R$timeini/$timeend_ky/0/250 \
	-Ba10f1:"time [kyr]":/a50f10:"halite rate [10@+12@+ kg yr@+-1@+]":E \
	-K -O >> $ps
#HALITE:
awk '(substr($0,1,1)!="#"){s=$23; diff=s-saltant; if (diff==0 || diffant==0) print ">"; else print $1/1000,diff/($1-tant); saltant=$23; tant=$1; diffant=diff}' $file | psxy -JX -R -M -W2/255/0/0 -O -K >> $ps
awk '(substr($0,1,1)!="#"){s=$24; diff=s-saltant; if (diff==0 || diffant==0) print ">"; else print $1/1000,diff/($1-tant); saltant=$24; tant=$1; diffant=diff}' $file | psxy -JX -R -M -W2ta/190/60/10 -O -K >> $ps
tail $file | awk '{if (substr($1,1,1)!="#") {total_halite=($22+$23+$24)*1e12; total_gypsum=($16+$17+$18)*1e12}}\
END{printf("\ntotal precipitated halite = %.2e kg\ntotal precipitated gypsum = %.2e kg\n", total_halite, total_gypsum)}' 

gmtset BASEMAP_FRAME_RGB +0/0/0 
awk '{if (substr($1,1,1)!="#") print $1/1000,$5}' $file  | psxy -JX -R$timeini/$timeend_ky/0/70 \
	-Ba10f1:"":/a20f10:"flow velocity V [m/s]":nsW \
	-W3 -Y-5 -K -O -H2 >> $ps
gmtset BASEMAP_FRAME_RGB +200/0/0 
awk '{if (substr($1,1,1)!="#") print $1/1000,$7/1e3}' $file  | psxy -JX -R$timeini/$timeend_ky/0/20 \
	-Ba10f1:"":/a5f1:"width W [km]":E \
	-W3 -Y0 -K -O -H2 >> $ps

gmtset BASEMAP_FRAME_RGB +0/0/0 
awk '{if (substr($1,1,1)!="#") print $1/1000,$6/1e3}' $file  | psxy -JX -R$timeini/$timeend_ky/0/1000 \
	-B:"":/a200f100:"discharge Q [10@+3@+ m@+3@+/s]":nsW \
	-W3 -Y-5 -K -O -H2 >> $ps
gmtset BASEMAP_FRAME_RGB +200/0/0 
awk '{if (substr($1,1,1)!="#") print $1/1000,$8}' $file  | psxy -JX -R$timeini/$timeend_ky/0/1 \
	-B:"":/a.5f.1:"erosion rate [m/yr]":E \
	-W3 -Y0 -K -O -H2 >> $ps

gmtset BASEMAP_FRAME_RGB +0/0/0 
awk '{if (substr($1,1,1)!="#") print $1/1000,$2}' $file  | psxy -JX -R$timeini/$timeend_ky/-230/12 \
	-Ba$timelabelint\f$timetickint\:"time t [kyr]":/a100f10:"sill, ocean z@-s@-, z@-0@- [m]":nSW \
	-W3 -Y-5 -K -O -H2 >> $ps
gmtset BASEMAP_FRAME_RGB +0/0/0 
awk '{if (substr($1,1,1)!="#") print $1/1000,$13}' $file | psxy -JX -R -W2   -O -K -H2 >> $ps
gmtset BASEMAP_FRAME_RGB +200/0/0 
awk '{if (substr($1,1,1)!="#") print $1/1000,$14}' $file | psxy -JX -R$timeini/$timeend_ky/-2800/100 \
	-B:"":/a500f100:"Med. level z@-1@- z@-2@- [m]":E -W3 -O -K -H2 >> $ps
awk '{if (substr($1,1,1)!="#") print $1/1000,$15}' $file | psxy -JX -R -W3to -O -H2 >> $ps

gmtset BASEMAP_FRAME_RGB +0/0/0 


awk '(substr($0,1,1)!="#"){gypsum=$16+$17+$18; halite=$22+$23+$24;} END{print "gypsum,halite [10e18 kg]", gypsum/1e6, halite/1e6}' $file

open $ps &

rm *.tmp

