def std_name(infile, outfile):
    fr = open(infile, 'r')
    fw = open(infile, 'w')
    for line in fr:
        if line.startswith('>'):
            fw.write(line.split("/")[0] + "\n")
        else:
            fw.write(line)


std_name(r"D:\usa\alignment_re\usa_HA_alignmet_re.fasta", r"D:\usa\alignment\usa_HA_alignmet2_re.fasta")
std_name(r"D:\usa\alignment_re\usa_MP_alignmet_re.fasta", r"D:\usa\alignment\usa_MP_alignmet2_re.fasta")
std_name(r"D:\usa\alignment_re\usa_NA_alignmet_re.fasta", r"D:\usa\alignment\usa_NA_alignmet2_re.fasta")
std_name(r"D:\usa\alignment_re\usa_NP_alignmet_re.fasta", r"D:\usa\alignment\usa_NP_alignmet2_re.fasta")
std_name(r"D:\usa\alignment_re\usa_NS_alignmet_re.fasta", r"D:\usa\alignment\usa_NS_alignmet2_re.fasta")
std_name(r"D:\usa\alignment_re\usa_PA_alignmet_re.fasta", r"D:\usa\alignment\usa_PA_alignmet2_re.fasta")
std_name(r"D:\usa\alignment_re\usa_PB1_alignmet_re.fasta", r"D:\usa\alignment\usa_PB1_alignmet2_re.fasta")
std_name(r"D:\usa\alignment_re\usa_PB2_alignmet_re.fasta", r"D:\usa\alignment\usa_PB2_alignmet2_re.fasta")
