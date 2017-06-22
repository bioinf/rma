#!/usr/bin/env python

import sys
import os

fastafile, logfile = sys.argv[1:]

annot = []
with open(logfile, 'r') as lf:
    for line in lf:
        content = line.strip().split('\t')
        annot.append(content[0] + ' ' + content[2] + content[1] + content[3])

with open(fastafile, 'r') as aasub:
    for line in aasub:
        if not line.startswith('>'):
            if gene_name + ' ' + ref_aa + position + alt_aa in annot:
                continue
            with open(fastafile + 'query.fa', 'w') as qfile:
                qfile.write('>' + gene_name + '\n' + line)
            with open(fastafile + 'pph.input', 'w') as pphin:
                pphin.write('\t'.join([gene_name, position, ref_aa, alt_aa]) + '\n')
            status_code = os.system('/Molly/barbitoff/software/polyphen-2.2.2/bin/run_pph.pl -s ' + fastafile + 'query.fa ' + fastafile + 'pph.input' + ' > ' + fastafile + 'pph.output 2> ' + fastafile + 'pph.log')
            with open(fastafile + 'pph.output', 'r') as ppho:
                for line_out in ppho:
                        print line_out.strip()
            continue
        line = line.strip()
        gene_name = line.split()[0][1:]
        position = line.split()[1][1:-1]
        ref_aa = line.split()[1][0]
        alt_aa = line.split()[1][-1]
        sys.stdout.flush()
