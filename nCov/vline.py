import matplotlib.pyplot as plt
import re

log = 'log-orf1ab'
seg = 'orf1ab'

def read_log(_in):
    sacs = []
    snames = []
    pidents = []
    pstarts = []
    pends = []
    with open(_in) as f:
        for line in f:
            row = line.split('\t')
            attrs = row[1].split("|")
            sacs.append(attrs[0])
            snames.append(attrs[1])
            pidents.append(float(row[2]))
            pstarts.append(int(re.findall(r'\d*\.\.\d*', row[0])[0].split('..')[0]))
            pends.append(int(re.findall(r'\d*\.\.\d*', row[0])[0].split('..')[1]))
    return sacs, snames, pidents, pstarts, pends 

sacs, snames, pidents, pstarts, pends = read_log("/home/zeng/Desktop/%s" % log)

f = lambda x: x if x in ['RaTG13', 'bat_SL_CoVZXC21', 'bat_SL_CoVZC45', 'other'] else 'other'
snames2 = list(map(f, snames))

cmap = {}
cmap['RaTG13'] = 'r'
cmap['bat_SL_CoVZXC21'] = 'g'
cmap['bat_SL_CoVZC45'] = 'c'
cmap['other'] = 'b'

lll = ['RaTG13', 'bat_SL_CoVZXC21', 'bat_SL_CoVZC45', 'other']

color_list = [cmap[x] for x in snames2]

plt.figure(figsize=(16, 9))
for _, sname, pident, pstart, pend, colo in zip(sacs, snames2, pidents, pstarts, pends, color_list):
    plt.hlines(pident, pstart, pend, label=sname, color=colo)
#plt.legend(lll, loc='lower right')
plt.text(200, 10, 'RaTG13', color='r', bbox=dict(facecolor='r', alpha=0.5))
plt.text(200, 20, 'bat_SL_CoVZXC21', color='g', bbox=dict(facecolor='g', alpha=0.5))
plt.text(200, 30, 'bat_SL_CoVZC45', color='c', bbox=dict(facecolor='c', alpha=0.5))
plt.text(200, 40, 'other', color='b', bbox=dict(facecolor='b', alpha=0.5))

plt.ylabel('% Identity')
plt.title("%s gene" % seg.upper())
plt.savefig("/home/zeng/Desktop/test-result-of-%s.jpg" % seg)
#plt.show()