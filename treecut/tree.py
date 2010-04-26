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
import fisher

try:
    from numpy import mean as lmean
    from scipy.stats.stats import ttest_ind as lttest_ind
except:
    try:
        from statlib.stats import lttest_ind, lmean 
    except:
        print >>sys.stderr, "Install either scipy or statlib for statistics calculations"


class ExtTree(list):

    def __init__(self, node, values, all, datatype="continuous"):

        self.node = node
        self.values = values
        for n in node.children:
            if not n.is_leaf():
                self.append(ExtTree(n, values, all, datatype=datatype))

        nset = set(node.iter_leaves())
        oset = all - nset

        # values for the direct children
        self.a = a = self.get_values(nset, values)
        # values for non-children (sibs) 
        self.b = b = self.get_values(oset, values)

        self.val = self.hi_min = self.lo_min = 1.0

        if a and b: 
            self.val = self.stat_test(a, b, datatype=datatype)

        # core dynamic programming
        self.lomin()
        self.himin()


    def __str__(self):
        return "%d\t%d\t%.1f\t%.1f\t%.1g\t%.1g\t%.1g" % (\
                len(self.a), len(self.b), lmean(self.a), lmean(self.b), \
                self.val, self.hi_min, self.lo_min)
        
    
    def stat_test(self, a, b, datatype="continuous"):
        if datatype=="continuous":
            return lttest_ind(a, b)[1]
        else:
            #print >>sys.stderr, "treating the values as discrete types .."
            a1, b1 = a.count(1), b.count(1)
            a0, b0 = a.count(0), b.count(0)
            return fisher.pvalue(a1, a0, b1, b0).two_tail


    def render(self, image_name, **kwargs):
        from draw import Dendrogram
        d = Dendrogram(self)
        d.savefig(image_name, **kwargs)
        print >>sys.stderr, "tree image saved to %s" % image_name


    def get_values(self, leaf_set, values):
        res = []
        for x in leaf_set:
            if x not in leaf_set:
                print >>sys.stderr, "%s missing in listfile" % x
            else:
                res.append(values[x.name])
        return res


    def get_all_nodes(self):
        res = []
        for e in self:
            res.append(e)
            res += e.get_all_nodes()
        return res


    def get_candidates(self, cutoff=.05):
        candidates = []
        for e in self:
            if e.val < min(e.lo_min, e.hi_min, cutoff):
                candidates.append(e)
            else:
                candidates += e.get_candidates(cutoff=cutoff)
        return candidates


    def print_all_nodes(self, filehandle):
        for i, e in enumerate(self.get_all_nodes()):
            print >>filehandle, "%d\t%s" % (i, e)


    def print_candidate(self, filehandle, cutoff=.05):
        for i, e in enumerate(self.get_candidates(cutoff=cutoff)):
            desc = "lo" if lmean(e.a) < lmean(e.b) else "hi" 
            print >>filehandle, "%s\t%s\t%.1f\t%.1g" % (
                ",".join(x.name for x in e.node.get_leaves()), desc,
                lmean(e.a), e.val)


    def himin(self):
        for e in self:
            e.hi_min = min(self.val, self.hi_min)
            e.himin()


    def lomin(self):
        if len(self)!=0: 
            self.lo_min = min([x.lomin() for x in self] +\
                    [x.val for x in self])
        return self.lo_min

