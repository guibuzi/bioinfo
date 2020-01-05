import os, sys, json
import re


def proc_tree(tree, attr_dict):
    pattern = r'EPI_ISL_[0-9]*'
    labels = re.findall(pattern, tree)
    for label in labels:
        cluster = attr_dict[label]['cluster']
        date = attr_dict[label]['date'].split("T")[0]
        region = attr_dict[label]['region2']
        repl = "%s_%s_%s_%s" % (label, cluster, region, date)
        tree = re.sub(label, repl, tree)
    return tree

def read_attr(path):
    with open(path) as f:
        data = json.loads(f.read())
    return data


def read_tree(path):
    with open(path) as f:
        data = f.read()
    return data


def main(seg):
    tree = read_tree("/home/zeng/Desktop/sample/tree/sample_%s.tree" % seg)
    attr_dict = read_attr("/home/zeng/Desktop/sample_index")
    new_tree = proc_tree(tree, attr_dict)
    with open("/home/zeng/Desktop/%s_tree.tree" % seg, "w") as f:
        f.write("#NEXUS\nBEGIN trees;\nTREE 'TREE_%s' =  %sEND;" % (seg, new_tree))


if __name__ == "__main__":
    main("HA")
    main("M1")
    main("M2")
#    main("MP")
    main("NA")
    main("NP")
#    main("NS")
    main("NS1")
    main("PA")
    main("PB1")
    main("PB2")