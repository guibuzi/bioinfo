import json


def read_fasta(file_path):
    """
    describtion: 读取 fasta 文件

    parameters:
    file_path: fasta 文件所在路径
    
    return: list of dict
    """
    data = {}
    with open(file_path) as f:
        contents = f.read()

    seq_lists = contents.split(sep='>')[1:]
    for seq_list in seq_lists:
        a = seq_list.split('\n')
        epi_seg = a[0].split('|')[0].strip()
        seq = "".join(a[1:]) 
        data[epi_seg] = seq
    return data


if __name__ == '__main__':
    with open("/home/zeng/Desktop/H3N2/data/Isolation_information/ioslation_information.json") as f:
        contents = f.read()

    info_dict = json.loads(contents)

    want_info = {key: value for key, value in info_dict.items() if 'United States' in value['location']}
    want_isl = set(want_info.keys())

    protein_list = {'HA': 'HA', 'NA': 'NA', 'PA': 'PA', 'PB1': 'PB1', 'PB2': 'PB2',
                    'M1': 'MP', 'NP': 'NP', 'NS1': 'NS'}

    for key, value in protein_list.items():
        fw = open("/home/zeng/Desktop/US_oct/data/usa_{}.fasta".format(key), 'w')

        fasta_data = read_fasta("/home/zeng/Desktop/H3N2/data/Protein/protein_{}.fasta".format(key))

        for item in want_isl:
            seg_id = info_dict[item]['{}'.format(value)]
            sequence = fasta_data.get(seg_id)
            if sequence:
                date = want_info[item]['date']
                fw.write(">" + item + "_" + date + "\n")
                fw.write(sequence + "\n")
        fw.close()

