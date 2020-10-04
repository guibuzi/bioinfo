from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO

from io import StringIO
import subprocess
import sys

# Align module
# prepare sequence
records = SeqIO.parse("tmp.fasta", "fasta")

# prepare stdin
handle = StringIO()
SeqIO.write(records, handle, "fasta")
data = handle.getvalue()
handle.close()

# prepare muscle
muscle_cline = MuscleCommandline()
stdout, stderr = muscle_cline(stdin=data)
align = AlignIO.read(StringIO(stdout), "fasta")
print(align)



# blast Module
query = 'TACAAAAGTGTGAATATCACTTTTGAACTTGACGAAAGGATTGATAAGGTACTTAATGAGAAGTGCTCCACCTATACAGTTGAATTCGGTACAGAAGTAAATGAGTTTGCTTGTGTTGTGGCAGATGCTGTCATAAAAACTTTACAACCAGTATCTGAATTACTTACACCACTGGGCATT'

seq = Seq(query)
seq_record = SeqRecord(seq, id='My_seq', description='descriptionsfsdfaf')


handle2 = StringIO()
SeqIO.write(seq_record, handle2, "fasta")

blastn_cline = NcbiblastnCommandline(db='all-cov/all-cov', evalue=0.001, outfmt=5, num_alignments=10)
stdout, stderr = blastn_cline(stdin=handle2.getvalue())

blast_result = SearchIO.read(StringIO(stdout), 'blast-xml')
print(blast_result)
