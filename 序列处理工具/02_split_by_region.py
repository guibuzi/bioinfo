import json


if __name__ == '__main__':
    with open("/home/zeng/Desktop/H3N2/data/Isolation_information/ioslation_information.json") as f:
        contents = f.read()

    json_file = json.loads(contents)
    info_dict = dict(json_file)

    want_info = {key: value for key, value in info_dict.items() if 'United States' in value['location']}
    want_isl = set(want_info.keys())

    protein_list = {'HA': 'HA', 'NA': 'NA', 'PA': 'PA', 'PB1': 'PB1', 'PB2': 'PB2',
                    'M1': 'MP', 'NP': 'NP', 'NS1': 'NS'}
    
    for key, value in protein_list.items():    
        fw = open("/home/zeng/Desktop/usa/usa_{}.fasta".format(key), "w")
        fr = open("/home/zeng/Desktop/H3N2/data/Protein/protein_{}.fasta".format(key))
        content = fr.read()
        items = content.split(">")[1:]
        for item in items:
            isl = item.split("\n")[0].split(" | ")[1]
            if isl in want_isl:
                date = want_info[isl]['date']
                seq = item.split("\n")[1]
                fw.write(">" + isl + " | " + date + "\n")
                fw.write(seq + "\n")
        fw.close()
        fr.close()
