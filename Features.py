import re
import sys
import math
from Bio import SeqIO

"""

"""
AA = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
      "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

NT = ["A","T","C","G"]

"""

"""
def CalculateAAComposition(aaseq):
    results = {}
    for i in AA:
        results[i] = 0 if aaseq.count(i) == 0 else round(float(aaseq.count(i)) / float(len(aaseq)), 5)
    return results

"""

"""
def CalculateNTFrequence(ntseq):
    results = {}
    for i in NT:
        results[i] = 0 if ntseq.count(i) == 0 else round(float(ntseq.count(i)) / float(len(ntseq)), 5)
    return results    
    

if __name__ == "__main__":
    fnt = str(sys.argv[1])
    faa = str(sys.argv[2])

    # Declare hash
    fastant = {}
    fastaaa = {}
    handle = open(fnt, "r")
    for record in SeqIO.parse(handle, "fasta"):
        fastant[record.id] = str(record.seq).replace('X','')
    handle.close()

    handle = open(faa, "r")
    for record in SeqIO.parse(handle, "fasta"):
        fastaaa[record.id] = str(record.seq).replace('X','')
    handle.close()

    output = open('Features.csv','w')
    #***********************************************************************#
    #
    #***********************************************************************#
    for key in fastant:
        res = CalculateNTFrequence(fastant[key])
        output.write(key + ',')
        
        # Complete sequence
        for r in res:
            output.write(str(res[r]) + ',')
        
        # n terminal
        medium = int (len(fastant[key]) / 2)
        res = CalculateNTFrequence(fastant[key][0:medium])
        for r in res:
            output.write(str(res[r]) + ',')

        # c terminal
        medium = int (len(fastant[key]) / 2)
        res = CalculateNTFrequence(fastant[key][medium::])
        for r in res:
            output.write(str(res[r]) + ',')

        res = CalculateAAComposition(fastaaa[key])
        for r in res:
            output.write(str(res[r]) + ',')
        output.write('\n') 
    output.close()