# run this script only when you want to rebuil a database !!!

# change the path where you want to build the database
path = "/home/zeng/INFLUENZA_DATABASE/H3N2/data"

# don't change !!!
segments = ['BM2', 'HA', 'M1', 'M2', 'NA', 'NB', 'NEP', 'NP', 'NS1', 'NS2', 'PA', 'PB1-F2', 'PB1', 'PB2']

# init the isolation information json file
f_info = open("/home/zeng/INFLUENZA_DATABASE/H3N2/data/isolation_information.json", "w")
f_info.close()

# init the segments sequences id file
segment_id = open("/home/zeng/INFLUENZA_DATABASE/H3N2/data/Protein/segment_id.txt", "w")
segment_id.close()

for seg in segments:
    f_protein_ha = open("{0}/Protein/protein_{1}.fasta".format(path, seg), 'w')
    f_protein_ha.close()
