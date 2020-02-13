import sys, os
from multiprocessing import Pool


def muscle(_in, _out):
    print("Begin to align %s" % os.path.split(_in)[1])
    os.system("muscle -in %s -out %s -quiet" % (_in, _out))
    print("Finsh align %s" % os.path.split(_in)[1])


def mafft(_in, _out):
    print("Begin to align %s" % os.path.split(_in)[1])
    os.system("mafft --thread 8 --quiet %s > %s" % (_in, _out))
    print("Finsh align %s" % os.path.split(_in)[1])

def generate_path(path, out_path):
    os.system("rm -rf %s/align" % path)
    tasks = os.listdir(path)
    in_seqs = []
    out_seqs = []

    for task in tasks:
        _in = os.path.join(path,  task)
        _out = os.path.join(out_path, task)
        in_seqs.append(_in)
        out_seqs.append(_out)

    os.makedirs("%s/align" % path, exist_ok=True)
    return in_seqs, out_seqs

def batch_mafft(inseqs, outseqs):
    for i, o in zip(inseqs, outseqs):
        mafft(i, o)

def batch_muscle(in_seqs, out_seqs, num_threads=4):
    pool = Pool(processes = num_threads)
    print("Begin tasks ...")
    for i, o in zip(in_seqs, out_seqs):
        pool.apply_async(muscle, (i, o))
    pool.close()
    pool.join()
    print("Finsh muscle.")

def main2(path, out_path):
    in_seqs, out_seqs = generate_path(path, out_path)
    batch_mafft(in_seqs, out_seqs)

def main(path, out_path):
    in_seqs, out_seqs = generate_path(path, out_path)
    batch_muscle(in_seqs, out_seqs)


if __name__ == "__main__":
#    main("/home/zeng/Desktop/sample_aus", "/home/zeng/Desktop/sample_aus/align")
#    main("/home/zeng/Desktop/sample_europe", "/home/zeng/Desktop/sample_europe/align")
#    main("/home/zeng/Desktop/sample_seasia", "/home/zeng/Desktop/sample_seasia/align")
#    main("/home/zeng/Desktop/sample_usa", "/home/zeng/Desktop/sample_usa/align")
    _, _in, _out = sys.argv
    main2(_in, _out)
