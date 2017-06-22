#!/usr/bin/env python


import os
import re
import sys


aasubfile, alreadyscored = sys.argv[1:]

scored = []
with open(alreadyscored, 'r') as asfile:
    for line in asfile:
        scored.append(' '.join(line.strip().split()[:2]))


with open(aasubfile, 'r') as aasub:
    for line in aasub:
        if line.startswith('>'):
            gene = line.split()[0][1:]
            substitution = line.split()[-1]
            continue
        if gene + ' ' + substitution in scored:
            continue
        with open(aasubfile + 'query.fa', 'w') as q:
            q.write('>' + gene + '\n' + line.strip())
        with open(aasubfile + 'query.var', 'w') as v:
            v.write(substitution)
        status_code = os.system('/Molly/barbitoff/software/sift6.0.1/bin/SIFT_for_submitting_fasta_seq.csh ' + aasubfile + 'query.fa /Molly/predeus/for_SIFT/uniref90/uniref90.fa ' + aasubfile + 'query.var > ' + aasubfile + '.log')
        if status_code == 0:
            with open('/Molly/barbitoff/software/sift6.0.1/tmp/' + aasubfile + 'query.SIFTprediction', 'r') as resfile:
                for line in resfile:
                    if substitution in line:
                        print '\t'.join([gene] + re.findall('\S+', line))
                        sys.stdout.flush()
