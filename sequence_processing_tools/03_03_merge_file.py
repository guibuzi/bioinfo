import os, sys

def merge_file(path_i, path_o, segment):
    file_list = os.listdir(path_i)
    file_list = [x for x in file_list]

    fw = open(path_o + "{}.fasta".format(segment), "w")
    for i in file_list:
        fr = open(os.path.join(path_i, i))
        print(fr.name + "is adding to {}.fasta ...".format(segment))
        fw.write(fr.read())
        print(fr.name + "finishing adding!")
        fr.close()
    fw.close()


def merge_one_region(path_i, path_o):
    try:
        os.makedirs(path_o)
    except:
        pass
    segments = os.listdir(path_i)
    for segment in segments:
        sub_path = path_i + segment + "/"
        merge_file(sub_path, path_o, segment)


if __name__ == "__main__":
    tasks = ["USA", "Aus"]
    for task in tasks:
        path_i = "/home/zeng/Desktop/reassortment_analysis/out_align/%s/" % task
        path_o = "/home/zeng/Desktop/reassortment_analysis/align2/%s/" %task
        merge_one_region(path_i, path_o)
