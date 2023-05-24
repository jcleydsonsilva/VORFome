#!/usr/bin/env python
# coding: utf-8
# Author: Jose Cleydson F . Silva
# v 0.01

import sys
import re

"""

"""

def get_data(file_seq):
	f1 = open(file_seq, "r")
	fasta = {}
	sequencia = ''
	key = ''
	for line in f1:
		line = line.replace('\n', '')
		if('>' in line):
			key = line.replace('>', '')
			fasta[key.replace(' ','_').replace(',','')] = sequencia
			sequencia = ''
		else:
			sequencia = sequencia + line
	fasta[key.replace(' ','_').replace(',','')] = sequencia
	f1.close()
	return fasta

"""

"""

def rev_seq(seq):
	trans = []
	for i in seq:
		if i == 'A':
			trans.append('T')
		elif i == 'C':
			trans.append('G')
		elif i == 'G':
			trans.append('C')
		elif i == 'T':
			trans.append('A')
		else:
			trans.append(i)
	seq_rev = ''.join(trans)
	seq_rev = seq_rev[::-1]
	return seq_rev

"""

"""

def six_frame(key, sequence):
	frames = {}

	# Forward frame: Finding ATG
	for i in range(0, 3):
		ff = key
		frame = sequence[i:len(sequence)]
		ff += '_Frame_' + str(i+1)
		frames[ff] = frame

	sequence = rev_seq(sequence)
	for i in range(0, 3):
		ff = key
		frame = sequence[0:len(sequence)-i]
		ff += '_Frame_-' + str(i+1)
		frames[ff] = frame
	return frames

"""

"""

def get_start_codon(sequence):
	pos_codons = []
	p = re.compile('ATG')
	for m in p.finditer(sequence):
		pos_codons.append([m.start(), m.group()])
	return pos_codons

"""

"""

def get_stop_codon(sequence):
	pos_codons = []
	p = re.compile('TGA')
	for m in p.finditer(sequence):
		pos_codons.append([m.start(), m.group()])
	p = re.compile('TAG')
	for m in p.finditer(sequence):
		pos_codons.append([m.start(), m.group()])
	p = re.compile('TAA')
	for m in p.finditer(sequence):
		pos_codons.append([m.start(), m.group()])
	return pos_codons

"""

"""

def get_donor_site_GT(sequence):
	pos_donor_site_GT = []
	p = re.compile('GT')
	for m in p.finditer(sequence):
		pos_donor_site_GT.append(m.start())
	return pos_donor_site_GT

"""

"""

def get_acceptor_site_AG(sequence):
	pos_acceptor_site_AG = []
	p = re.compile('AG')
	for m in p.finditer(sequence):
		pos_acceptor_site_AG.append(m.start())
	return pos_acceptor_site_AG
"""

"""

def pirimidine (intron):
	# verifica a proporção de pirimidinas
	quantPirimidinas = 0
	quantPurinas = 0
	if len(intron) > 70:
		for pp in range(len(intron)-51, len(intron)):
			if intron[pp] == 'G' or intron[pp] == 'A':
				quantPirimidinas += 1
			else:
				quantPurinas += 1
		if quantPirimidinas > quantPurinas:
			return 0
		else:
			return 1
	return 1
"""

"""
def translate_rna(sequence):
    codon = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
            "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
            "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
            "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
            "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
            "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
            "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
            "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
            "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
            "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
            "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
            "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
            "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "T",
            "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
            "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
            "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}

    protein_seq = ''
    for n in range(0, len(sequence), 3):
        if sequence[n:n+3] in codon:
            protein_seq += codon[sequence[n:n+3]]
    return protein_seq

"""

"""
def main():
	#fasta = get_data('Sequence.fa')
	fasta = get_data(sys.argv[1])
	file_orfs = 'ORFs-cds.fasta'
	file_prot = 'ORFs-to-amino-acids-sequence.fasta'
	file_orfs_write = open(file_orfs, 'w')
	file_prot_write = open(file_prot, 'w')
	for key in fasta:
		seq_and_seq_rev = {}
		seq_and_seq_rev[key] = fasta[key]
		key_rev = key
		key_rev += '_Rev_Seq'
		seq_and_seq_rev[key_rev] = rev_seq(fasta[key])
		for key_f in seq_and_seq_rev:
			print (key_f)
			count = 0
			start_codons = get_start_codon(seq_and_seq_rev[key_f])
			stop_codons = get_stop_codon(seq_and_seq_rev[key_f])
			for stop in stop_codons:
				for start in start_codons:
					if ((start[0] < stop[0]) and (stop[0] - start[0]) > 120):
						if (((start[0] - start_previous) % 3 == 0) and (stop_previous == stop[0])):
							start_previous = start[0]
							stop_previous = stop[0]
						else:
							seq = seq_and_seq_rev[key_f][start[0]:stop[0]+3]
							donor_site_GT = get_donor_site_GT(seq)
							acceptor_site_AG = get_acceptor_site_AG(seq)
							for gt in sorted(donor_site_GT):
								for ag in sorted(acceptor_site_AG):
									if (ag < gt) and (len(seq[ag:gt+2]) < 100) and (len(seq[ag:gt+2]) > 60):
										if pirimidine(seq[ag:(gt+2)]) == 0:
											#print (seq[ag:(gt+2)])                                            
											orf = str(seq[0:int(ag)] + seq[(gt+2):len(seq)])
											orfp = translate_rna(orf)
											if orfp.count('*') < 2 and orfp[len(orfp)-1] == '*' and len(orfp) > 30:           
												file_orfs_write.write('>'+key_f.replace(' ','_') +'_'+ str(count) + '_intron_' + str(ag) + '-' + str(gt)+ '\n' + orf + '\n')
												file_prot_write.write('>'+key_f.replace(' ','_') +'_'+ str(count) + '_intron_' + str(ag) + '-' + str(gt)+ '\n' + orfp + '\n')
							orfp = translate_rna(seq)
							if orfp.count('*') < 2 and orfp[len(orfp) -1] == '*' and orfp not in seq_preview:
								file_orfs_write.write(">"+ key_f +'_'+ str(count) + '\n')
								file_orfs_write.write(seq+'\n')                              
								file_prot_write.write(">"+ key_f +'_'+ str(count) +'\n')
								file_prot_write.write(orfp+'\n')
							count += 1
							seq_preview = orfp
					start_previous = start[0]
					stop_previous = stop[0]
	file_prot_write.close()
	file_orfs_write.close()
if __name__== "__main__":
	main()


# In[ ]:




