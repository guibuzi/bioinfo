import json
from collections import defaultdict

dict_isl = defaultdict(dict)
list_segment = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
with open(r"D:\isolate_information.csv", 'r') as f:
    lines = f.readlines()
    for line in lines:
        row = [x.strip() for x in line.split(',')]
        for num, segment in enumerate(list_segment):
            dict_isl[row[0]][segment] = row[num + 1]
        dict_isl[row[0]]['date'] = row[-3]
        dict_isl[row[0]]['name'] = row[-5]
        dict_isl[row[0]]['location'] = row[-4]


# ha_to_isl = {info['HA']: isl for isl, info in c.items() if info['HA']}
# na_to_isl = {info['HA']: isl for isl, infor in c.items() if info['NA']}


json = json.dumps(dict_isl, sort_keys=False, indent=4, separators=(',', ': '))

with open(r'd:\isolation_information.json', 'w') as g:
    g.write(json)
