#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
python %prog treefile listfile

treefile is a Newick-formatted file
listfile contains the accession=>value mapping, separated by comma

Python script that traverses through a hierarchical clustering tree
and calculate the significance values on all inner nodes and determine
the nodes that gives the least P-value
"""

import os.path as op
import sys
import ete2
from treecut.tree import ExtTree


def read_values(listfile, datatype="continuous"):
    import csv
    reader = csv.reader(file(listfile))
    reader.next() # header
    wrap = float if datatype=="continuous" else int
    return dict((acc, wrap(value)) for (acc, value) in reader)


if __name__ == '__main__':

    from optparse import OptionParser

    p = OptionParser(__doc__) 
    p.add_option("--discrete", default=False, action="store_true",
            help="are the data in listfile discrete values? [default:%default]")
    p.add_option("--cutoff", type="float", 
            default=.01, help="minimum P-value to report [default:%default]")
    p.add_option("--printall", action="store_true", default=False,
            help="print verbose information for all inner nodes [default:%default]")
    options, args = p.parse_args()

    if len(args) != 2:
        sys.exit(p.print_help())

    treefile, listfile = args
    datatype = "discrete" if options.discrete else "continuous"
    
    for f in (treefile, listfile):
        if not op.exists(f):
            parser.error("File %s not found" % f)

    # the tree topology
    tree = ete2.Tree(treefile)
    all = set(tree.iter_leaves())  # terminal nodes
    # value mappings
    values = read_values(listfile, datatype=datatype)
    n, m = len(all), len(values)
    if n!=m:
        print >>sys.stderr, "[warning] number of accessions don't match between treefile(%d) and listfile(%d)" % (n, m)

    # generate output
    fw = sys.stdout
    t = ExtTree(tree, values, all, datatype=datatype)

    if options.printall:
        # header
        print >>fw, ("node_id member_mean non-member_mean P-value " 
                "min_ancestor_P-value min_descendant_P-value").replace(" ", "\t") 
        t.print_all_nodes(fw)
    else:
        t.print_candidate(fw, cutoff=options.cutoff)

    t.render("tree.pdf")
