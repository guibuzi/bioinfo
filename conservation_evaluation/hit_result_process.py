fr = open("conservation_evaluation/blast_result_du.txt", 'r')
fw = open("conservation_evaluation/blast_hits.txt", 'w')
for line in fr:
    if not line.startswith("#"):
        fw.write(line)
fr.close()
fw.close()
