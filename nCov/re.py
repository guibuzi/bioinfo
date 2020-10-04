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



seq = SeqIO.read("SARS-CoV-2_S.fasta", "fasta")
sub_seq = seq[0:501]
# print(sub_seq, len(sub_seq))




handle2 = StringIO()
SeqIO.write(sub_seq, handle2, "fasta")

blastn_cline = NcbiblastnCommandline(db='all-cov/all-cov', outfmt=5, num_alignments=10, negative_seqidlist='SARS-CoV-2.id', task='blastn')
stdout, stderr = blastn_cline(stdin=handle2.getvalue())

blast_result = SearchIO.read(StringIO(stdout), 'blast-xml')
for hit in blast_result:
    print(hit.id)