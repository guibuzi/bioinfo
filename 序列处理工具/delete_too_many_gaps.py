import re


def count_sites(sequence):
    return len(sequence)


def count_gaps(sequence):
    regex = re.compile('-')
    return len(re.findall(regex, sequence))


def delete_sequence_with_gaps_threshold(in_file, out_file, threshold=0.5):
    fr = open(in_file, 'r', )
    fw = open(out_file, 'w')
    taxon = ''
    for line in fr:
        if line.startswith('>'):
            taxon = line
        else:
            length = count_sites(line)
            gap_num = count_gaps(line)
            if gap_num / length < threshold:
                fw.write(taxon)
                fw.write(line)
    fw.close()
    fr.close()


if __name__ == '__main__':
    delete_sequence_with_gaps_threshold(r"D:\china\alignment\China_HA_alignment2.fasta",
                                        r"D:\china\alignment_re\China_HA_alignment_re.fasta")
    delete_sequence_with_gaps_threshold(r"D:\china\alignment\China_M1_alignment2.fasta",
                                        r"D:\china\alignment_re\China_M1_alignment_re.fasta")
    delete_sequence_with_gaps_threshold(r"D:\china\alignment\China_M2_alignment2.fasta",
                                        r"D:\china\alignment_re\China_M2_alignment_re.fasta")
    delete_sequence_with_gaps_threshold(r"D:\china\alignment\China_NA_alignment2.fasta",
                                        r"D:\china\alignment_re\China_NA_alignment_re.fasta")
    delete_sequence_with_gaps_threshold(r"D:\china\alignment\China_NS_alignment2.fasta",
                                        r"D:\china\alignment_re\China_NS_alignment_re.fasta")
    delete_sequence_with_gaps_threshold(r"D:\china\alignment\China_NP_alignment2.fasta",
                                        r"D:\china\alignment_re\China_NP_alignment_re.fasta")
    delete_sequence_with_gaps_threshold(r"D:\china\alignment\China_PA_alignment2.fasta",
                                        r"D:\china\alignment_re\China_PA_alignment_re.fasta")
    delete_sequence_with_gaps_threshold(r"D:\china\alignment\China_PB1_alignment2.fasta",
                                        r"D:\china\alignment_re\China_PB1_alignment_re.fasta")
    delete_sequence_with_gaps_threshold(r"D:\china\alignment\China_PB2_alignment2.fasta",
                                        r"D:\china\alignment_re\China_PB2_alignment_re.fasta")
