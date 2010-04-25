#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Author: Haibao Tang <bao@uga.edu>

Python script that traverses through a hierarchical clustering tree
and calculate the significance values on all inner nodes and determine
the nodes that gives the least P-value
"""

import os.path as op
import sys
import ete2
from treecut.tree import ExtTree


if __name__ == '__main__':

    from optparse import OptionParser

    parser = OptionParser(usage="%prog [-t] treefile [-f] listfile\n" + __doc__) 
    parser.add_option("-t", "--treefile", type="str", 
            help="Newick-format tree file")
    parser.add_option("-f", "--listfile", type="str", 
            help="two-column mapping from OTUs to values")
    parser.add_option("-c", "--cutoff", type="float", 
            default=.01, help="minimum P-value to report [default:%default]")
    parser.add_option("-v", "--printall", action="store_true", default=False,
            help="print verbose information for all inner nodes [default:%default]")
    options, args = parser.parse_args()

    if not options.treefile:
        parser.error("You must specify -t")
    elif not op.exists(options.treefile):
        parser.error("File %s not found" % options.treefile)
    if not options.listfile:
        parser.error("You must specify -f")
    elif not op.exists(options.listfile):
        parser.error("File %s not found" % options.listfile)

    tree = ete2.Tree(options.treefile)
    #tree.render("render.pdf")
    all = set(tree.iter_leaves())  # terminal nodes

    # value mappings
    fp = file(options.listfile)
    values = {}
    for row in fp:
        a, b = row.split()
        values[a] = float(b)
    n, m = len(all), len(values)
    if n!=m:
        print >>sys.stderr, "[warning] number of OTUs don't match between treefile(%d) and listfile(%d)" % (n, m)

    # generate output
    fw = sys.stdout
    t = ExtTree(tree, values, all)
    t.himin()
    t.lomin()
    if options.printall:
        # header
        print >>fw, ("node_id member_mean non-member_mean P-value " 
                "min_ancestor_P-value min_descendant_P-value").replace(" ", "\t") 
        t.print_all_nodes(fw)
    else:
        t.print_candidate(fw, cutoff=options.cutoff)
