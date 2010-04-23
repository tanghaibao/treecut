#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Author: Haibao Tang <bao@uga.edu>
This script is freely available, modifiable

Python script that traverses through a hierarchical clustering tree
and calculate the significance values on all inner nodes and determine
the nodes that gives the least P-value
"""

import os
import sys
import newick
import newick.tree
from statlib.stats import lttest_ind, lmean 

class ExtTree:
    def __init__(self, node):
        self.node = node
        self.edges = []
        for n,b,l in node.get_edges():
            if not isinstance(n, newick.tree.Leaf):
                self.edges.append(ExtTree(n))
        nset = set(node.get_leaves_identifiers())
        oset = all - nset
        a = [values[x] for x in nset]
        b = [values[x] for x in oset]
        self.a, self.b = a, b
        self.val = self.hi_min = self.lo_min = 1.0
        if a and b: self.val = lttest_ind(a,b)[1]

    def __str__(self):
        return "%d\t%d\t%.1f\t%.1f\t%.1g\t%.1g\t%.1g" % (\
                len(self.a), len(self.b), lmean(self.a), lmean(self.b), \
                self.val, self.hi_min, self.lo_min)
        
    def get_all_children(self):
        res = []
        for e in self.edges:
            res.append(e)
            res += e.get_all_children()
        return res

    def print_node(self, filehandle):
        all_nodes = self.get_all_children()
        for i, e in enumerate(all_nodes):
            print >>filehandle, "%d\t%s" % (i, e)

    def print_candidate(self, filehandle, cutoff=.05):
        if self.val<min(self.lo_min, self.hi_min, cutoff):
            if lmean(self.a)<lmean(self.b): desc = "lo"
            else: desc = "hi"
            filehandle.write("%s\t%s\t%.1f\t%.1g\n"%(\
                ",".join(self.node.get_leaves_identifiers()), desc,
                lmean(self.a), self.val))
        for e in self.edges:
            e.print_candidate(filehandle)

    def himin(self):
        for e in self.edges:
            e.hi_min = min(self.val, self.hi_min)
            e.himin()

    def lomin(self):
        if self.edges: 
            self.lo_min = min([x.lomin() for x in self.edges] +\
                    [x.val for x in self.edges])
        return self.lo_min


if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [-t] treefile [-f] listfile", 
            version="%prog 0.1")
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
    elif not os.path.exists(options.treefile):
        parser.error("FILE %s not found"%options.treefile)
    if not options.listfile:
        parser.error("You must specify -f")
    elif not os.path.exists(options.listfile):
        parser.error("FILE %s not found"%options.listfile)

    fp = file(options.treefile)
    data = fp.read()
    tree = newick.parse_tree(data)
    all = set(tree.get_leaves_identifiers())  # terminal nodes

    # value mappings
    fp = file(options.listfile)
    values = {}
    for row in fp:
        a,b = row.split()
        values[a] = float(b)
    n, m = len(all), len(values.keys())
    assert n==m, \
            "number of OTUs don't match between treefile(%d) and listfile(%d)"%(n,m)

    # generate output
    fw = sys.stdout
    t = ExtTree(tree)
    t.himin()
    t.lomin()
    if options.printall:
        fw.write("node_id\tmember_mean\tnon-member_mean\tP-value\t"
                "min_ancestor_P-value\tmin_descendant_P-value\n") # header
        t.print_node(fw)
    else:
        t.print_candidate(fw, cutoff=options.cutoff)
    fw.close()
