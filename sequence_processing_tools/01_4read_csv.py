import json
from collections import defaultdict

def read_csv(file_path):
    dict_isl = defaultdict(dict)
    list_segment = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
    count = 0

    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if count == 0:
                count += 1
                continue
            else:
                row = [x.strip() for x in line.split(',')]
                for num, segment in enumerate(list_segment):
                    dict_isl[row[0]][segment] = row[num + 1].split("|")[0].strip()
                dict_isl[row[0]]['date'] = row[-1]
                dict_isl[row[0]]['name'] = row[-5]
                dict_isl[row[0]]['location'] = row[-3]
                count += 1

    ha_to_isl = {info['HA']: isl for isl, info in dict_isl.items() if info['HA']}
    na_to_isl = {info['NA']: isl for isl, info in dict_isl.items() if info['NA']}
    pb1_to_isl = {info['PB1']: isl for isl, info in dict_isl.items() if info['PB1']}
    pb2_to_isl = {info['PB2']: isl for isl, info in dict_isl.items() if info['PB2']}
    pa_to_isl = {info['PA']: isl for isl, info in dict_isl.items() if info['PA']}
    np_to_isl = {info['NP']: isl for isl, info in dict_isl.items() if info['NP']}
    mp_to_isl = {info['MP']: isl for isl, info in dict_isl.items() if info['MP']}
    ns_to_isl = {info['NS']: isl for isl, info in dict_isl.items() if info['NS']}
    json_file = json.dumps(dict_isl, sort_keys=False, indent=4, separators=(',', ': '))

    return json_file, ha_to_isl, na_to_isl, pb1_to_isl, pb2_to_isl, pa_to_isl, np_to_isl, mp_to_isl, ns_to_isl


if __name__ == "__main__":
    # old database
    original = read_csv("/home/zeng/Desktop/H3N2/data/Isolation_information/ioslation_information.csv")
    # with open("/home/zeng/Desktop/H3N2/data/Isolation_information/ioslation_information.json", 'w') as g:
    #     g.write(original[0])
    
    origianl2 = json.loads(original[0])
    already_isolation = {key for key in origianl2.keys()}
    print("There are already exist {} entry.".format(len(already_isolation)))

    # read new entry
    read = read_csv("/home/zeng/Desktop/US_oct/gisaid_epiflu_isolates.csv")
    read2 = json.loads(read[0])
    read3 = {key: value for key, value in read2.items() if key not in already_isolation}

    print("{} new entrys has been added.".format(len(read3)))

    # 更新数据库
    origianl2.update(read3)
    print("Now we have {} sequences.".format(len(origianl2)))
    read4 = json.dumps(origianl2, sort_keys=False, indent=4, separators=(',', ': '))

    with open("/home/zeng/Desktop/H3N2/data/Isolation_information/ioslation_information.json", 'a') as g:
        g.write(read4)
