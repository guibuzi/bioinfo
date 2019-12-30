import sys, os
from multiprocessing import Pool


def build_tree(file_i, file_o):
    print("Begin run task %s ..." % (os.path.basename(file_i)))
    os.system("fasttree -quiet %s > %s" % (file_i, file_o))
    print("Finish %s." % (os.path.basename(file_i)))


if __name__ == "__main__":
    path = "/home/zeng/Desktop/reassortment_analysis/align/"
    path2 = path + "tree/"

    to_tree = []
    out_tree = []

    task_dirs = os.listdir(path)
    for task_dir in task_dirs:
        sub_path = path + task_dir + "/"
        sub_path2 = path2 + task_dir + "/"
        os.makedirs(path2 + task_dir, exist_ok=True)
        tasks = os.listdir(sub_path)
        segments = [i.split(".")[0] for i in tasks]

        for task, segment in zip(tasks, segments):
            file_i = os.path.join(sub_path, task)
            file_o = sub_path2 + segment + ".tree"
            to_tree.append(file_i)
            out_tree.append(file_o)
    
    
    pool = Pool(processes = 10)
    print("Begin tasks ...")
    for i, o in zip(to_tree, out_tree):
        pool.apply_async(build_tree, (i, o))

    pool.close()
    pool.join()
    print("Finsh build tree.")