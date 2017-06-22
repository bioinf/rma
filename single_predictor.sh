#!/bin/bash

# Classes construction - 84th 1000Gp1; 100th - ExAC
cut -f 1,4 dbNSFP_narrow > SIFT
grep -v '1000Gp1' SIFT | grep -P '^T.*\d.*'  > false_class
grep -v '1000Gp1' SIFT | grep -P '^D.*\d.*'  > true_class

# Binning
awk '{ print int($2/0.01) }' true_class | sort -n | uniq -c > TC_binned
awk '{ print int($2/0.01) }' false_class | sort -n | uniq -c > FC_binned

# Cleanup and formatting
perl -pe 's|^ +||' TC_binned | awk -F ' ' '{ print $2, $1 }' | sed 's/ /\t/' | head -n 51 > TC.cleaned
perl -pe 's|^ +||' FC_binned | awk -F ' ' '{ print $2, $1 }' | sed 's/ /\t/' | head -n 51 > FC.cleaned

# Making final table
cut -f 2 FC.cleaned | paste TC.cleaned - > cc.SIFT.classified.cleaned

# RMA stuff - 8th 1000G, 9th ExAC
# cut -d ',' -f 8 ../final_file/predictedRMA | awk '{ print int($1/0.01) }' | sort -n | uniq -c | perl -pe 's|^ +||' | awk -F ' ' '{ print $2, $1 }' | sed 's/ /\t/' > RMA.cleaned
