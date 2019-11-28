import os

def merge_file(path, segment):
    file_list = os.listdir(path)
    file_list = [x for x in file_list if x.startswith('alignment')]

    file_number = len(file_list)

    fw = open(path + "alignment_{}.fasta".format(segment), "w")
    for i in range(file_number):
        fr = open(path + "alignment_{}.fasta".format(i))
        print(fr.name + "is adding to alignment_{}.fasta ......".format(segment))
        fw.write(fr.read())
        print(fr.name + "finishing adding!")
        fr.close()
    fw.close()

merge_file("/home/zeng/Desktop/US_oct/data/pb2/", 'pb2')

