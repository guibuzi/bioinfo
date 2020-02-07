import sys, os
from multiprocessing import Pool


def build_tree(file_i, file_o):
    print("Begin run task %s ..." % (os.path.basename(file_i)))
    os.system("fasttree -quiet %s > %s" % (file_i, file_o))
    print("Finish %s." % (os.path.basename(file_i)))


def generate_path(path, out_path):
    os.system("rm -rf %s/tree" % path)
    in_seqs = []
    out_tree = []
    tasks = os.listdir(path)

    for task in tasks:
        _in = os.path.join(path, task)
        _out = os.path.join(out_path, task.replace('fasta', 'tree'))
        in_seqs.append(_in)
        out_tree.append(_out)

    os.makedirs("%s/tree" % path, exist_ok=True)
    return in_seqs, out_tree


def batch_tree(to_tree, out_tree, num_threads=6):
    pool = Pool(processes = 6)
    print("Begin tasks ...")
    for i, o in zip(to_tree, out_tree):
        pool.apply_async(build_tree, (i, o))
    pool.close()
    pool.join()
    print("Finsh build tree.")


def main(path, out_path):
    in_seqs, out_tree = generate_path(path, out_path)
    batch_tree(in_seqs, out_tree)


if __name__ == "__main__":
    main("/home/zeng/Desktop/sample_aus/align", "/home/zeng/Desktop/sample_aus/align/tree")
    main("/home/zeng/Desktop/sample_europe/align", "/home/zeng/Desktop/sample_europe/align/tree")
    main("/home/zeng/Desktop/sample_seasia/align", "/home/zeng/Desktop/sample_seasia/align/tree")
    main("/home/zeng/Desktop/sample_usa/align", "/home/zeng/Desktop/sample_usa/align/tree")