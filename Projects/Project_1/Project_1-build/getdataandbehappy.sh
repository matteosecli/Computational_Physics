#!/bin/bash

./Project_1 10000 onealg 1

#-nojvm 
matlab -nodesktop -nosplash -r "cd $(pwd); graphsolution"

rm X.txt
rm U.txt
rm E.txt


pretime=':'
posttime='s.'
for (( i=10; i <= 10000; i*=2 )); do
	for j in {1..10}; do
		./Project_1 $i > times.txt
		time_lu=`cat times.txt | head -1`
		time_gaus=`cat times.txt | head -2 | tail -1`
#		time_trid=`cat times.txt | head -3 | tail -1`
#		time_special=`cat times.txt | tail -1`
		time=${time_lu#*$pretime}
		time=${time%%$posttime*}
		echo $time >> time_lu_$i.txt
		time=${time_gaus#*$pretime}
		time=${time%%$posttime*}
		echo $time >> time_gaus_$i.txt
#		time=${time_trid#*$pretime}
#		time=${time%%$posttime*}
#		echo $time >> time_trid_$i.txt
#		time=${time_special#*$pretime}
#		time=${time%%$posttime*}
#		echo $time >> time_special_$i.txt
		rm times.txt
	done
done

for (( i=10; i <= 10000; i*=2 )); do
	mean=`cat time_lu_$i.txt|head -10|tr " " "\t"|cut -f13|awk 'f += $1 {printf("%.13g\n",f/NR)}'|tail -1`
	if [ "$mean" == "" ]; then
		mean=0
	fi
	echo $mean >> time_lu.txt
	rm time_lu_$i.txt

	mean=`cat time_gaus_$i.txt|head -10|tr " " "\t"|cut -f13|awk 'f += $1 {printf("%.13g\n",f/NR)}'|tail -1`
	if [ "$mean" == "" ]; then
		mean=0
	fi
	echo $mean >> time_gaus.txt
	rm time_gaus_$i.txt

#	mean=`cat time_trid_$i.txt|head -10|tr " " "\t"|cut -f13|awk 'f += $1 {printf("%.13g\n",f/NR)}'|tail -1`
#	if [ "$mean" == "" ]; then
#		mean=0
#	fi
#	echo $mean >> time_trid.txt
#	rm time_trid_$i.txt

#	mean=`cat time_special_$i.txt|head -10|tr " " "\t"|cut -f13|awk 'f += $1 {printf("%.13g\n",f/NR)}'|tail -1`
#	if [ "$mean" == "" ]; then
#		mean=0
#	fi
#	echo $mean >> time_special.txt
#	rm time_special_$i.txt
done





for (( i=1000000; i <= 200000000; i*=2 )); do
	for j in {1..10}; do
		./Project_1 $i > times.txt
		time_trid=`cat times.txt | head -1`
		time_special=`cat times.txt | tail -1`
		time=${time_trid#*$pretime}
		time=${time%%$posttime*}
		echo $time >> time_trid_$i.txt
		time=${time_special#*$pretime}
		time=${time%%$posttime*}
		echo $time >> time_special_$i.txt
		rm times.txt
	done
done

for (( i=1000000; i <= 200000000; i*=2 )); do
	mean=`cat time_trid_$i.txt|head -10|tr " " "\t"|cut -f13|awk 'f += $1 {printf("%.13g\n",f/NR)}'|tail -1`
	if [ "$mean" == "" ]; then
		mean=0
	fi
	echo $mean >> time_trid.txt
	rm time_trid_$i.txt

	mean=`cat time_special_$i.txt|head -10|tr " " "\t"|cut -f13|awk 'f += $1 {printf("%.13g\n",f/NR)}'|tail -1`
	if [ "$mean" == "" ]; then
		mean=0
	fi
	echo $mean >> time_special.txt
	rm time_special_$i.txt
done




preerr=':'
posterr='%'
for (( i=10; i <= 300000000; i*=2 )); do
	./Project_1 $i onealg 0 > errors_temp.txt
	error=`cat errors_temp.txt | tail -1`
	error=${error#*$preerr}
	error=${error%%$posterr*}
	echo -e "$i\t$error" >> errors.txt
	rm errors_temp.txt
done

matlab -nodesktop -nosplash -r "cd $(pwd); grapherr"

rm errors.txt


