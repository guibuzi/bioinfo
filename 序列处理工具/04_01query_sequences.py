import json


def read_query_table(path):
    data = []
    with open(path) as f:
        for line in f:
            isl_epi = line.split(',')[0]
            seg_epi = line.split(',')[1].strip()
            data.append((isl_epi, seg_epi))
    return data


def read_info_json(path):
    with open(path) as f:
        data = json.loads(f.read())
    return data


def read_fasta(path):
    with open(path) as f:
        result = {}
        for line in f:
            if line.startswith(">"):
                label = line[1:-1]
                result[label] = ""
            else:
                result[label] += line.strip()
    return result



if __name__ == '__main__':
    to_query = read_query_table("/home/zeng/python_work/bioinfo/mosaic\sampling.txt")
    info = read_info_json("/home/zeng/Desktop/H3N2/data/Isolation_information/ioslation_information.json")
    fasta = read_fasta("/home/zeng/Desktop/H3N2/data/Protein/protein_HA.fasta")
    to_query = [x[0] for x in to_query]
    fw = open("/home/zeng/python_work/bioinfo/mosaic/us_ha_sampling.fasta", "w")
    for query in to_query:
        seg_id = info[query]['HA']
        date = info[query]['date']
        name = info[query]['name']
        location = info[query]['location']
        annotation = "|".join([name, date])
        sequence = fasta.get(seg_id)
        if sequence:
            fw.write(">" + query + "|" + annotation + "\n")
            fw.write(sequence + "\n")
