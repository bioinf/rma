#!/usr/bin/env python

import os
import sys
import re


snpfile, annfile, upfile, gcfile, gcnuclfile = sys.argv[1:]

genetic_code = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
    "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G", '': 'NOAA'}

dna = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def revComp(seq):
    compseq = ''
    for i in seq:
        compseq += dna[i]
    return compseq[::-1]

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

# Extracting CDS corresponding to each transcript in a dictionary
# Structure of the dictionary: GENE_NAME -> {'ID': [('start', 'end'), ...], ...}
# For each gene name construct a list of transcripts with a list of corresponding CDSs
cds_dict = {}
strandedness = {}
with open(annfile, 'r') as an:
    for line in an:
        content = line.split('\t')
        if len(re.findall('\tCDS\t', line)) == 0 or 'protein_coding' not in line:
            continue
        transcript_id = re.findall("transcript_id \"(ENST[R\d]+)", line)[0]
        gene_name = re.findall('gene_name \"([^\"]+)', line)[0]
        if gene_name not in cds_dict:
            strandedness[gene_name] = content[6]
            cds_dict[gene_name] = {transcript_id: [(int(content[3]), int(content[4]))]}
            continue
        cds_dict[gene_name][transcript_id] = cds_dict[gene_name].get(transcript_id, []) + [(int(content[3]), int(content[4]))]

# print cds_dict['PRAMEF19']
# sys.exit()

enst = {}
with open(gcfile, 'r') as gc:
    for line in gc:
        if line.startswith('>'):
            id = re.findall('ENST[R\d]+', line)[0]
        else:
            enst[id] = line.strip()

enst_nuc = {}
with open(gcnuclfile, 'r') as tran:
    for line in tran:
        if line.startswith('>'):
            id = re.findall('ENST[R\d]+', line)[0]
            for thing in line.split('|'):
                if 'CDS:' in thing:
                    cds = thing[4:]
                    cds_lower, cds_upper = [int(x) - 1 for x in re.findall('[0-9]+', cds)]
        else:
            tseq = line.strip()[cds_lower:cds_upper + 1]
            enst_nuc[id] = tseq

def get_letter(nucl, gn):
    if strandedness[gn] == '+':
        return nucl
    else:
        return dna[nucl]

with open(snpfile, 'r') as var:
    for line in var:
        if line.startswith('#'):
            continue
        content = line.strip().split('\t')
        if 'EFF' not in line:
            print 'FAILED', content[2]
            continue
        gene_name = content[7].split('|')[4]
        if gene_name not in cds_dict or gene_name not in uniprot:
            print 'FAILED', content[2]
            continue
        transcripts = cds_dict[gene_name]
        for transcript in transcripts:
            if transcript not in enst:
                continue
            found = False
            for ind, cds in enumerate(transcripts[transcript]):
                if max(cds) >= int(content[1]) >= min(cds):
                    transcript_pos = (sum([max(elt) - min(elt) + 1 for elt in transcripts[transcript][:ind]]), sum([max(elt) - min(elt) + 1 for elt in transcripts[transcript][:ind]]) + max(cds) - min(cds) + 1)
                    found = True
                    remainder = 0
                    break
            if not found:
                continue        
            if not (enst_nuc[transcript].startswith('ATG') and enst[transcript][0] == 'M'):
                for remainder in range(4):
                    newseq = '' 
                    for i in range(3):
                        newseq += genetic_code[enst_nuc[transcript][i * 3 + remainder:i * 3 + remainder + 3]]
                    if newseq in enst[transcript][:6]:
                        remainder = 3 - remainder if remainder > 0 else 0
                        break
            exon_aaseq = enst[transcript][transcript_pos[0]/3:(transcript_pos[1]/3 + 1)]
            if len(exon_aaseq) > len(uniprot[gene_name]):
                continue
            with open('exon.fa', 'w') as e:
                e.write('>exon\n' + exon_aaseq)
            with open('up.fa', 'w') as up:
                up.write('>uniprot\n' + uniprot[gene_name])
            status_code = os.system("/Molly/barbitoff/software/fasta-36.3.8c/bin/ssearch36 -d 1 -m 10 -p exon.fa up.fa > results.al 2> /dev/null")
            if status_code != 0:
                continue
            with open('results.al', 'r') as r:
                all = r.read()
                al_exonseq = ''.join(all.split('>exon ..')[1].split('al_display_start')[1].split('\n')[1:]).split('>uniprot')[0]
                al_upseq = ''.join(all.split('>uniprot ..')[1].split(';')[6].split('\n')[1:])
                # With respect to UniProt, yup
                al_start = re.findall('\d+', all.split('>uniprot ..')[1].split(';')[6])[0]
                exon_al_start = re.findall('\d+', all.split('>exon ..')[1].split(';')[4])[0]
                real_al_start = re.findall('\d+', all.split('>uniprot ..')[1].split(';')[4])[0]
                exon_al_display_start = int(re.findall('\d+', all.split('>exon ..')[1].split(';')[6])[0])
                al_cons = ''.join(all.split('al_cons')[1].split('>')[0].split('\n')[1:])
            if strandedness[gene_name] == '+':
                nucpos = int(content[1]) - min(cds) + transcript_pos[0]
            else:
                nucpos =  max(cds) - int(content[1]) + transcript_pos[0]
            if remainder > nucpos:
                continue
            exonpos = (nucpos + remainder)/3 - transcript_pos[0]/3
            pos_offset = 0
            concordant = False
            confirmed = False
            discordant = False
            exonic_aa_num = exon_al_display_start - 1
            old_codon = enst_nuc[transcript][((nucpos + remainder)/3) * 3 - remainder:((nucpos + remainder)/3) * 3 + 3 - remainder]
            new_codon = old_codon[:(nucpos + remainder) % 3] + get_letter(content[4], gene_name) + old_codon[(nucpos + remainder) % 3 + 1:] if len(content[4]) == 1 else ''
            if len(old_codon) != 3 or len(new_codon) != 3:
                break
            for i, aa in enumerate(al_cons):
                if al_exonseq[i] == '-' or al_exonseq[i] == ' ':
                    continue
                else:
                    if exonic_aa_num == exonpos and (aa == '.' or aa == ' '):
                        if genetic_code[new_codon] != al_upseq[i] and gene_name != 'PCSK9':
                            discordant = True
                            break
                        else:
                            confirmed = True
                            break
                    elif exonic_aa_num == exonpos:
                        concordant = True
                        break
                    exonic_aa_num += 1
            if al_upseq.startswith('-'):
                offset_shift = len(re.findall('^\-+', al_upseq)[0])
            else:
                offset_shift = 0
            offset = str(int(al_start) + i - len(re.findall('-', al_upseq[:i+1])) + offset_shift)
            if confirmed:
                print 'CONFIRMED', content[2], 'IN GENE', gene_name, 'TRANSCRIPT ID', transcript, genetic_code[new_codon] + offset + genetic_code[old_codon]
                break
            if discordant and genetic_code[new_codon] != genetic_code[old_codon]:
                print 'DISCORDANT', content[2], 'IN GENE', gene_name, genetic_code[new_codon] + offset + genetic_code[old_codon]
                break
            if concordant:
                if 'NON_SYNONYMOUS_CODING' in line or 'PROTEIN_INTERACTION_LOCUS' in line:
                    print 'CONCORDANT', content[2], 'IN GENE', gene_name, genetic_code[new_codon] + offset + genetic_code[old_codon]
                    break
        if not confirmed and not concordant and not discordant:
            print 'FAILED', content[2]
        sys.stdout.flush()
