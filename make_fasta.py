#!/usr/bin/env python

import os
import sys
import re

logfile, upfile, ofile = sys.argv[1:]

uniprot = {}
with open(upfile, 'r') as up:
    pseq = ''
    for line in up:
        if line.startswith('>'):
            if pseq != '':
                uniprot[genename] = pseq
                pseq = ''
            genename = re.findall('GN=(\S+)', line)[0] if 'GN=' in line else 'TRASH'
        else:
            pseq += line.strip()
    uniprot[genename] = pseq

with open(logfile, 'r') as log:
    for line in log:
        if 'FAILED' in line:
            continue
        content = line.split()
        gene_name = content[4]
        prot_sequence = uniprot[gene_name]
        sub = content[-1]
        offset = int(re.findall('\d+', sub)[0]) - 1
        if (prot_sequence[offset] != sub[0] and prot_sequence[offset] != sub[-1]) and 'DISCORDANT' not in line:
            print line
            continue
        prot_sequence = prot_sequence[:offset] + sub[0] + prot_sequence[offset + 1:]
        with open(ofile, 'a') as outp:
            outp.write('>' + gene_name + ' ' + sub + '\n' + prot_sequence + '\n')

