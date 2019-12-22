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
                    dict_isl[row[0]][segment] = row[num + 1].split("|")[0].strip()   # when one isolation have more than one sequences for one segment then keep the first
                dict_isl[row[0]]['date'] = row[-3]
                dict_isl[row[0]]['location'] = row[-5]
                dict_isl[row[0]]['lineage'] = row[-7]
                dict_isl[row[0]]['name'] = row[-9]
                count += 1
    return dict_isl


def update_database(file_name, old_database):
    already_isolation = set(old_database.keys())
    print("There are already exist {} entry.".format(len(already_isolation)))
    # read new entry
    read = read_csv(file_name)
    read_new = {key: value for key, value in read.items() if key not in already_isolation}
    print("{} new entrys has been added.".format(len(read_new)))
    # 更新数据库
    old_database.update(read_new)
    print("Now we have {} sequences\n".format(len(old_database)))
    return old_database


def read_isl(path):
    with open(path) as g:
        contents = g.read()
    if contents:
        return json.loads(contents)
    else:
        return defaultdict(dict)


if __name__ == "__main__":
    # read the old json file
    old_database = read_isl("/home/zeng/INFLUENZA_DATABASE/H3N2/data/isolation_information.json")

    # change this line when add new isolation
    # you can add new isolations from multiple csv file
    old_database = update_database("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/ISL_info/ioslate_file (1).csv", old_database)
    old_database = update_database("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/ISL_info/ioslate_file (2).csv", old_database)
    old_database = update_database("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/ISL_info/ioslate_file (3).csv", old_database)
    old_database = update_database("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/ISL_info/ioslate_file (4).csv", old_database)
    old_database = update_database("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/ISL_info/ioslate_file (5).csv", old_database)
    old_database = update_database("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/ISL_info/ioslate_file (6).csv", old_database)
    old_database = update_database("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/ISL_info/ioslate_file (7).csv", old_database)
    old_database = update_database("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/ISL_info/ioslate_file (8).csv", old_database)
    old_database = update_database("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/ISL_info/ioslate_file (9).csv", old_database)
    old_database = update_database("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/ISL_info/gisaid_epiflu_isolates.csv", old_database)

    # when finish read all the new isolations, re-write the json file
    json_file = json.dumps(old_database, sort_keys=False, indent=4, separators=(',', ': '))
    with open("/home/zeng/INFLUENZA_DATABASE/H3N2/data/isolation_information.json", 'w') as g:
        g.write(json_file)
