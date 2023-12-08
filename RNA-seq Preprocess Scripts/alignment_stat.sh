#!/bin/bash

## Get row names from one file
FILE=$(ls MapFromRaw/ | grep .final.out| head -1)
grep umber MapFromRaw/$FILE | grep -v splices | sed -e 's/\s\+/./g' | column -s "|" -t | awk '{print $1}' | sed '1i Sample' | paste -s -d "," > alignments.csv

## Get the count
for file in MapFromRaw/*.final.out
do
sample=$(basename $file)
grep umber $file | grep -v splices | awk '{print $NF}' | sed "1i ${sample}" | paste -s -d "," >> alignments.csv
done
