#!/usr/bin/env python

"""
@author camilla eldridge 
"""

from Bio import SeqIO
from Bio.Blast import NCBIWWW
import sys

query=sys.argv[1]
blast_out=sys.argv[2]
evalue=sys.argv[3]


my_query = SeqIO.read(query, format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", my_query.seq)
blast_result = open(blast_out, "w")
blast_result.write(result_handle.read())
blast_result.close()
result_handle.close()


from Bio.Blast import NCBIXML
with open(blast_out) as result_handle:
    blast_record = NCBIXML.read(result_handle)

eval_thresh = float(evalue)

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
         if hsp.expect < eval_thresh:
                print('****Alignment****')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('e value:', hsp.expect)
                print( "hsp.gaps", hsp.gaps)
                print( "hsp.strand", hsp.strand)
                print( "hsp.score", hsp.score)
                print( "hsp.frame", hsp.frame)
                print( "hsp.query_start", hsp.query_start)
                print( "hsp.query_end", hsp.query_end)
                print( "hsp.num_alignments", hsp.num_alignments)
                print( "hsp.identities", hsp.identities)
                print( "hsp.positives", hsp.positives)
                print(hsp.query[0:100] + '...')
                print(hsp.match[0:100] + '...')
                print(hsp.sbjct[0:100] + '...')
                print('Subject start', hsp.sbjct_start)
                print( "hsp.sbjct_end", hsp.sbjct_end)
                print(">" + alignment.title + "\n" + hsp.query)
                print(">" + my_query.id + "\n" + hsp.sbjct)
    print( "-----------------------------------"   )
