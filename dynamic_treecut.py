#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Author: Haibao Tang <bao@uga.edu>

Python script that traverses through a hierarchical clustering tree
and calculate the significance values on all inner nodes and determine
the nodes that gives the least P-value
"""

import os
import os.path as op
import sys
import newick
import newick.tree
try:
    from numpy import mean as lmean
    from scipy.stats.stats import ttest_ind as lttest_ind
except:
    from statlib.stats import lttest_ind, lmean 
finally:
    print >>sys.stderr, "Install either statlib or scipy for statistics calculations"


class ExtTree(list):
    """
    This is similar to the newick.tree.Tree, but each node has an associated P-value
    this allows easy propagation of P-values either ascending or descending the tree.
    """

    def __init__(self, node, values, all):
        self.node = node
        for n,b,l in node.get_edges():
            if not isinstance(n, newick.tree.Leaf):
                self.append(ExtTree(n, values, all))

        nset = set(node.get_leaves_identifiers())
        oset = all - nset

        # values for the direct children
        self.a = a = self.get_values(nset, values)
        # values for non-children (sibs) 
        self.b = b = self.get_values(oset, values)

        self.val = self.hi_min = self.lo_min = 1.0

        if a and b: 
            self.val = lttest_ind(a, b)[1]

    def __str__(self):
        return "%d\t%d\t%.1f\t%.1f\t%.1g\t%.1g\t%.1g" % (\
                len(self.a), len(self.b), lmean(self.a), lmean(self.b), \
                self.val, self.hi_min, self.lo_min)
        
    
    def get_values(self, leaf_set, values):
        res = []
        for x in leaf_set:
            if x not in leaf_set:
                print >>sys.stderr, "%s missing in listfile" % x
            else:
                res.append(values[x])
        return res


    def get_all_children(self):
        res = []
        for e in self:
            res.append(e)
            res += e.get_all_children()
        return res


    def print_all_nodes(self, filehandle):
        all_nodes = self.get_all_children()
        for i, e in enumerate(all_nodes):
            print >>filehandle, "%d\t%s" % (i, e)


    def print_candidate(self, filehandle, cutoff=.05):
        if self.val < min(self.lo_min, self.hi_min, cutoff):
            desc = "lo" if lmean(self.a) < lmean(self.b) else "hi" 
            print >>filehandle, "%s\t%s\t%.1f\t%.1g" % (
                ",".join(self.node.get_leaves_identifiers()), desc,
                lmean(self.a), self.val)

        for e in self:
            e.print_candidate(filehandle)


    def himin(self):
        for e in self:
            e.hi_min = min(self.val, self.hi_min)
            e.himin()


    def lomin(self):
        if len(self)!=0: 
            self.lo_min = min([x.lomin() for x in self] +\
                    [x.val for x in self])
        return self.lo_min



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

    fp = file(options.treefile)
    data = fp.read()
    tree = newick.parse_tree(data)
    all = set(tree.get_leaves_identifiers())  # terminal nodes

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
