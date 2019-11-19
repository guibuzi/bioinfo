path = "/home/zeng/Desktop/H3N2/data/Protein/"
segments = ['BM2', 'HA', 'M1', 'M2', 'NA', 'NB', 'NEP', 'NP', 'NS1', 'NS2', 'PA', 'PB1-F2', 'PB1', 'PB2']
# f_info = open("{}isolation_information.json".format(path), "w")
# f_info.close()
for seg in segments:
    f_protein_ha = open("{0}protein_{1}.fasta".format(path, seg), 'w')
    f_protein_ha.close()
