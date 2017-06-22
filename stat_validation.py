#!/usr/bin/env python

import re
import sys
import os
from numpy import random as npr

dbnsfpfile, rmafile, niter = sys.argv[1:]

numlines = 0
dbnsfp_pred = {}
with open(dbnsfpfile, 'r') as dbns:
    for line in dbns:
        if line.startswith('#'):
            continue
        numlines += 1
        content = line.split('\t')
        exacaf = float(content[99]) if content[99] != '.' else 0
        binN = int(exacaf/0.01)
        siftpred = content[28].split(',')[0]
        pphpred = content[34].split(',')[0]
        provpred = content[58].split(',')[0]
        if siftpred == 'D' and (pphpred == 'P' or pphpred == 'D') and provpred == 'D':
            dbnsfp_pred[binN] = dbnsfp_pred.get(binN, []) + [True]
        else:
            dbnsfp_pred[binN] = dbnsfp_pred.get(binN, []) + [False]
        if numlines/100000 > 0:
            print 'Processed 100k records'
            numlines = 0
            sys.stdout.flush()

for i in dbnsfp_pred:
    print sum(dbnsfp_pred[i])/(len(dbnsfp_pred[i]) * 1.0), 'at bin', 0.01*i
sys.exit()

for xxx in range(int(niter)):
    print 'Beginning iteration ' + str(xxx)
    sys.stdout.flush()
    truth = []
    with open(rmafile, 'r') as rma:
        for line in rma:
            if line.startswith('Chromosome,'):
                continue
            content = line.split(',')
            siftpred = content[13]
            pphpred = content[15]
            provpred = content[17]
            if siftpred != '' and pphpred != '' and provpred != '':
                binN = int(float(content[8])/0.01)
                while binN not in dbnsfp_pred:
                    binN += 1
                truth.append(bool(npr.choice(dbnsfp_pred[binN], 1)))
        print(sum(truth))
