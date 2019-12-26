import os, sys

def merge_file(path, segment):
    file_list = os.listdir(path)
    file_list = [x for x in file_list]

    fw = open(path + "{}.fasta".format(segment), "w")
    for i in file_list:
        fr = open(path + i)
        print(fr.name + "is adding to {}.fasta ......".format(segment))
        fw.write(fr.read())
        print(fr.name + "finishing adding!")
        fr.close()
    fw.close()



if __name__ == "__main__":
    path1 = "/home/zeng/Desktop/reassortment_analysis/out_align/USA/"
    path2 = "/home/zeng/Desktop/reassortment_analysis/out_align/USA/"
    segments = os.listdir(path1)
    for segment in segments:
        sub_path = path1 + segment + "/"
        merge_file(sub_path, segment)

