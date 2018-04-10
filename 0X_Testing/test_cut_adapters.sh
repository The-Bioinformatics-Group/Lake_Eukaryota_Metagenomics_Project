#!/bin/bash






LIMIT=3
DATADIR=C:/Users/Karamech/Documents/Utbildning/Kandidatarbete/Data
DIRS=$DATADIR/Sample_*
mkdir -p $DATADIR/AdaptersRemoved

for D in $DIRS; ##if D is a directory with name containing "Sample_" run cut adapt for all files ending with ".gz"
do
	for f in $D/*.fastq.gz;
	do
    	cp $f $DATADIR/AdaptersRemoved
    	FILE=${f#$D}
    	if [[ $FILE =~ R1 ]]
    	then
    		echo "R1"
    	elif [[ $FILE =~ R2 ]]
    	then
    		echo "R2"
    	fi
    	DIRNAME="${D}"
    done
    if [ $LIMIT -eq 1 ]; 
    then 
    	break 
    fi
 	let LIMIT-=1
done
