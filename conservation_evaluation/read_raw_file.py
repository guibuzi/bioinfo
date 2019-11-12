fr = open("conservation_evaluation/to_search.txt", 'r')
fw = open("conservation_evaluation/to_search.fas", 'w')

seq_list = []
count = 0
for line in fr:
    seq, label = line.split()
    if seq not in seq_list:
        count += 1
        fw.write(">" + label + "_{}\n".format(count))
        fw.write(seq + "\n")
        seq_list.append(seq)

fr.close()
fw.close()


f = open("conservation_evaluation/to_search.fas", 'r')
label_list = []
for line in f:
    if line.startswith(">"):
        label_list.append(line[1:])

print(len(label_list))
print(len(set(label_list)))
f.close()