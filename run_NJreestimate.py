import copy
import sys
from optparse import OptionParser
import treeswift as ts
from sys import platform as _platform
import tempfile
from subprocess import Popen, PIPE
import pkg_resources
from NJreestimate.readfq import readfq


def nodecopy(node):
    res = ts.Node()
    res.edge_length = 0
    if node.is_leaf():
        res.label = node.label
    else:
        copied_chd = [nodecopy(chd) for chd in node.children]
        for chd in copied_chd:
            res.add_child(chd)
    return res


def get_labels(node):
    return [i.label for i in node.traverse_postorder(internal=False)]


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-t", "--tree", dest="tree_fp",
                      help="path to the reference tree", metavar="FILE")
    parser.add_option("-o", "--output", dest="output_fp",
                      help="path for the output newick file",
                      metavar="FILE")
    parser.add_option("-s", "--ref", dest="ref_fp",
                      help="path to the reference alignment file (FASTA), containing reference sequences",
                      metavar="FILE")
    parser.add_option("-p", "--protein", dest="protein_seqs", action='store_true', default=False,
                      help="input sequences are protein sequences")

    (options, args) = parser.parse_args()

    with open(options.ref_fp, "r") as f:
        seqdict = dict()
        for name, seq, qual in readfq(f):
            if seq in seqdict:
                seqdict[seq].append(name)
            else:
                seqdict[seq] = [name]
    assert len(seqdict) > 3

    alltags = []
    for k, v in seqdict.items():
        alltags += v


    orig_branch_tree_full = ts.read_tree_newick(options.tree_fp)

    orig_branch_tree_with_dups = orig_branch_tree_full.extract_tree(alltags, without=False, suppress_unifurcations=True)
    orig_branch_tree_with_dups.resolve_polytomies()
    orig_branch_tree_with_dups.is_rooted = False

    leafs_with_dups_map = dict()
    for n in orig_branch_tree_with_dups.traverse_postorder(internal=False):
        leafs_with_dups_map[n.label] = n

    uniqtags = {}
    for k, v in seqdict.items():
        if len(v) == 1:
            uniqtags[v[0]] = (k, None)
            continue
        visited = {i: False for i in v}
        all_choices = {}
        for i in v:
            if visited[i]:
                continue
            else:
                visited[i] = True

            # starting from i, traverse upwards in the tree.
            # for each node tested, check if all labels under it
            # are in the group v.
            # if yes, update the mrca as the current node
            #
            mrca = leafs_with_dups_map[i]
            seed = mrca
            mrca_count = 1

            for a in seed.traverse_ancestors(include_self=False):
                sibling_labels = get_labels(a)
                if all([i in visited for i in sibling_labels]):
                    mrca = a
                    mrca_count = len(sibling_labels)
                    for i in sibling_labels:
                        visited[i] = True
                else:
                    break
            all_choices[i] = (mrca_count, mrca)
        max_key = max(all_choices, key=lambda x: (all_choices[x][0], x))
        uniqtags[max_key] = (k, all_choices[max_key][1])

    uniq_ref = tempfile.NamedTemporaryFile(delete=True, mode='w+t').name
    with open(uniq_ref, "w") as f:
        for k, v in uniqtags.items():
            f.write(">" + k + "\n" + v[0] + "\n")



    orig_branch_tree = orig_branch_tree_full.extract_tree(uniqtags, without=False, suppress_unifurcations=True)
    orig_branch_tree.resolve_polytomies()
    orig_branch_tree.is_rooted = False
    orig_branch_resolved_fp = tempfile.NamedTemporaryFile(delete=True, mode='w+t').name
    orig_branch_tree.write_tree_newick(orig_branch_resolved_fp)



    if _platform == "darwin":
        fasttree_exec = pkg_resources.resource_filename('NJreestimate', "tools/FastTree-darwin")
    elif _platform == "linux" or _platform == "linux2":
        fasttree_exec = pkg_resources.resource_filename('NJreestimate', "tools/FastTree-linux")
    elif _platform == "win32" or _platform == "win64" or _platform == "msys":
        fasttree_exec = pkg_resources.resource_filename('NJreestimate', "tools/FastTree.exe")
    else:
        # Unrecognised system
        raise ValueError('Your system {} is not supported yet.' % _platform)

    bb_fp = tempfile.NamedTemporaryFile(delete=True, mode='w+t')
    fasttree_log = tempfile.NamedTemporaryFile(delete=True, mode='w+t').name

    s = [fasttree_exec, "-nosupport", "-nome", "-noml", "-log", fasttree_log,
         "-intree", orig_branch_resolved_fp]
    if not options.protein_seqs:
        s.append("-nt")
    with open(uniq_ref, "r") as rf:
        with Popen(s, stdout=PIPE, stdin=rf, stderr=sys.stderr) as p:
            tree_string = p.stdout.read().decode('utf-8')
            print(tree_string)

    uniqs_tree_nj = ts.read_tree(tree_string, schema="newick")
    leaf_map = dict()
    for n in uniqs_tree_nj.traverse_postorder(internal=False):
        leaf_map[n.label] = n

    for k, v in uniqtags.items():
        seq, mrca = v
        if len(seqdict[seq]) == 1:
            continue

        existing = leaf_map[k]
        o_node_parent = existing.parent
        o_node_parent.remove_child(existing)
        mrca_copy = nodecopy(mrca)
        mrca_copy.edge_length = existing.edge_length
        o_node_parent.add_child(mrca_copy)


    uniqs_tree_nj.write_tree_newick(options.output_fp, hide_rooted_prefix=True)

