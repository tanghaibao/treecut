#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
python %prog treefile listfile [imagefile]

treefile is a Newick-formatted file
listfile contains the accession=>value mapping, separated by comma
optional [imagefile] will generate an image (.svg, .png, .jpg, .pdf, etc. are supported)

Python script that traverses through a hierarchical clustering tree
and calculate the significance values on all inner nodes and determine
the nodes that gives the least P-value
"""

import os.path as op
import sys
import csv
import ete2
from treecut.tree import ExtTree


def read_values(listfile, datatype="continuous"):
    reader = csv.reader(file(listfile))
    wrap = float if datatype=="continuous" else (lambda x:x.split(";")) 
    return dict((acc, wrap(value)) for (acc, value) in reader if acc[0]!="#")


if __name__ == '__main__':

    from optparse import OptionParser

    p = OptionParser(__doc__) 
    p.add_option("--discrete", default=False, action="store_true",
            help="are the data in listfile discrete classes "
            "(use Fisher's exact test to calculate P-values) [default:%default (use t-test)]")
    p.add_option("--cutoff", type="float", 
            default=.01, help="minimum P-value to report [default:%default]")
    p.add_option("--printall", action="store_true", default=False,
            help="print verbose information for all inner nodes [default:%default]")
    options, args = p.parse_args()

    if len(args)==2:
        treefile, listfile = args
        outfile = None
    elif len(args)==3:
        treefile, listfile, outfile = args
    else:
        sys.exit(p.print_help())

    datatype = "discrete" if options.discrete else "continuous"
    
    for f in (treefile, listfile):
        if not op.exists(f):
            p.error("File %s not found" % f)

    # the tree topology
    tree = ete2.Tree(treefile)
    # value mappings
    values = read_values(listfile, datatype=datatype)

    tree_accs = set(x.name for x in tree.iter_leaves())  # terminal nodes
    list_accs = set(values.keys())

    for x in tree_accs - list_accs:
        print >>sys.stderr, "[warning] %s missing in listfile" % x

    for x in list_accs - tree_accs:
        print >>sys.stderr, "[warning] %s missing in treefile" % x

    n, m = len(tree_accs), len(values)
    if n!=m:
        print >>sys.stderr, "[warning] number of accessions don't match between treefile(%d) and listfile(%d)" % (n, m)

    # generate output
    fw = sys.stdout
    t = ExtTree(tree, values, tree_accs, datatype=datatype)

    if options.printall:
        # header
        print >>fw, "\t".join(t.verbose_fields)
        t.print_all_nodes(fw)
    else:
        t.print_modules(fw, cutoff=options.cutoff)

    if outfile:
        t.render(outfile, cutoff=options.cutoff, dpi=80)

