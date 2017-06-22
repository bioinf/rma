#!/bin/bash

# Classes construction - 84th 1000Gp1; 100th - ExAC
cut -f 29,32,59,84 ../../software/dbNSFP2.9.txt > dbNSFP_narrow
grep -v '1000Gp1' dbNSFP_narrow | awk '{ if ($1 != "." && $2 != "." && $3 != "." && $4 != "." ) print }' | grep -P 'D.*?\t[PD].*?\tD.*?\t' - > true_class
grep -v '1000Gp1' dbNSFP_narrow | awk '{ if ($1 != "." && $2 != "." && $3 != "." && $4 != "." ) print }' | grep -v -P 'D.*?\t[PD].*?\tD.*?\t' - > false_class

# Binning
awk '{ print int($4/0.01) }' true_class | sort -n | uniq -c > TC_binned
awk '{ print int($4/0.01) }' false_class | sort -n | uniq -c > FC_binned

# Cleanup and formatting
perl -pe 's|^ +||' TC_binned | awk -F ' ' '{ print $2, $1 }' | sed 's/ /\t/' | head -n 51 > TC.cleaned
perl -pe 's|^ +||' FC_binned | awk -F ' ' '{ print $2, $1 }' | sed 's/ /\t/' | head -n 51 > FC.cleaned

# Making final table
cut -f 2 FC.cleaned | paste TC.cleaned - > cc.dbNSFP.classified.cleaned

# RMA stuff - 8th 1000G, 9th ExAC
cut -d ',' -f 8 ../final_file/predictedRMA | awk '{ print int($1/0.01) }' | sort -n | uniq -c | perl -pe 's|^ +||' | awk -F ' ' '{ print $2, $1 }' | sed 's/ /\t/' > RMA.cleaned
