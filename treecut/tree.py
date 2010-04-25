#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Author: Haibao Tang <bao@uga.edu>

ExtTree class is similar to the newick.tree.Tree, but each node has an associated P-value
this allows easy propagation of P-values either ascending or descending the tree.
"""

import os.path as op
import sys
import ete2

try:
    from numpy import mean as lmean
    from scipy.stats.stats import ttest_ind as lttest_ind
except:
    try:
        from statlib.stats import lttest_ind, lmean 
    except:
        print >>sys.stderr, "Install either scipy or statlib for statistics calculations"


class ExtTree(list):

    def __init__(self, node, values, all):
        self.node = node
        for n in node.children:
            if not n.is_leaf():
                self.append(ExtTree(n, values, all))

        nset = set(node.iter_leaves())
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
                res.append(values[x.name])
        return res


    def get_all_children(self):
        res = []
        for e in self:
            res.append(e)
            res += e.get_all_children()
        return res


    def print_all_nodes(self, filehandle):
        for i, e in enumerate(self.get_all_children()):
            print >>filehandle, "%d\t%s" % (i, e)


    def print_candidate(self, filehandle, cutoff=.05):
        if self.val < min(self.lo_min, self.hi_min, cutoff):
            desc = "lo" if lmean(self.a) < lmean(self.b) else "hi" 
            print >>filehandle, "%s\t%s\t%.1f\t%.1g" % (
                ",".join(x.name for x in self.node.get_leaves()), desc,
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



