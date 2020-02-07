import re

def read_attr(path):
    data = {}
    with open(path, "r") as f:
        for line in f:
            label = line.split(" ")[0]
            attr = "".join(line.split(" ")[1:]).replace("\n", "")
            attr = attr[1:-1]
            data[label] = attr
    return data


def read_tree(path):
    with open(path) as f:
        data = f.read()
    return data


def proc_tree(tree, attr):
    pattern = r'EPI_ISL_[0-9]*'
    labels = re.findall(pattern, tree)
    new_tree = tree
    for label in labels:
        _attr = attr[label]
        repl = label + _attr
        new_tree = re.sub(label, repl, new_tree)
    return new_tree


def write_tree(path, tree):
    with open(path, "w") as f:
        f.write("#NEXUS\nbegin trees;\n\ttree TREE1 = [&R] %s\nend;" % tree)


def main():
    attrs = read_attr("/Users/jinfeng/Downloads/align/sample_index")
    tree = read_tree("/Users/jinfeng/Downloads/align/tree/sample_HA.tree")
    new_tree = proc_tree(tree, attrs)
    write_tree("/Users/jinfeng/Downloads/test_tree.tree", new_tree)

if __name__ == '__main__':
    main()
