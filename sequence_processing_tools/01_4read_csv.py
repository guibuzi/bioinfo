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
                dict_isl[row[0]]['location'] = row[-2]
                dict_isl[row[0]]['lineage'] = row[-3]
                dict_isl[row[0]]['name'] = row[-5]
                count += 1

    return dict_isl


def update_database(file_name, old_database):
    already_isolation = set(old_database.keys())
    print("There are already exist {} entry.".format(len(already_isolation)))
    # read new entry
    read = read_csv("/home/zeng/Desktop/H1N1/data/original_data/isl_info/{}".format(file_name))
    read_new = {key: value for key, value in read.items() if key not in already_isolation}

    print("{} new entrys has been added.".format(len(read_new)))

    # 更新数据库
    old_database.update(read_new)
    print("Now we have {} sequences\n".format(len(old_database)))
    return old_database


if __name__ == "__main__":
    # old database
    old_database = read_csv("/home/zeng/Desktop/H1N1/data/original_data/isl_info/gisaid_epiflu_isolates.csv")

    old_database = update_database("gisaid_epiflu_isolates (1).csv", old_database)
    old_database = update_database("gisaid_epiflu_isolates (2).csv", old_database)
    old_database = update_database("gisaid_epiflu_isolates (3).csv", old_database)
    old_database = update_database("gisaid_epiflu_isolates (4).csv", old_database)
    old_database = update_database("gisaid_epiflu_isolates (5).csv", old_database)
    old_database = update_database("gisaid_epiflu_isolates (6).csv", old_database)
    old_database = update_database("gisaid_epiflu_isolates (7).csv", old_database)
    old_database = update_database("gisaid_epiflu_isolates (8).csv", old_database)
    old_database = update_database("gisaid_epiflu_isolates (9).csv", old_database)
    old_database = update_database("gisaid_epiflu_isolates (10).csv", old_database)
    old_database = update_database("gisaid_epiflu_isolates (11).csv", old_database)
    old_database = update_database("gisaid_epiflu_isolates (12).csv", old_database)

    json_file = json.dumps(old_database, sort_keys=False, indent=4, separators=(',', ': '))

    with open("/home/zeng/Desktop/H1N1/data/ioslation_information.json", 'a') as g:
        g.write(json_file)
