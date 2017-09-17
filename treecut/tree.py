#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
ExtTree class is similar to the ete2.Tree, but each node has an associated P-value
this allows easy propagation of P-values either ascending or descending the tree.
"""

import sys
from stats import stat_test, lmean


class ExtTree(list):
    __slots__ = ("node", "values", "values2", "datatype", "a", "b", "val",
            "hi_min", "lo_min", "note", "desc")

    def __init__(self, node, values, values2, accs, datatype="continuous"):

        self.node = node
        self.values = values
        self.values2 = values2
        self.datatype = datatype
        self.note = ""
        self.desc = ""

        for n in node.children:
            if not n.is_leaf():
                self.append(ExtTree(n, values, values2, accs, datatype=datatype))

        nset = set(x.name for x in node.iter_leaves())
        oset = accs - nset

        # values for the direct children
        self.a = a = self.get_values(nset, values)
        # values for non-children (sibs)
        self.b = b = self.get_values(oset, values)

        self.val = self.hi_min = self.lo_min = 1.0

        if a and b:
            sys.stderr.write(".")
            self.val, self.note = stat_test(a, b, datatype=datatype)

        # core dynamic programming
        self.lomin()
        self.himin()

    def __str__(self):
        return "%d\t%d\t%s\t%.1g\t%.1g\t%.1g" % (\
                len(self.a), len(self.b), self.note,
                self.val, self.hi_min, self.lo_min)

    def __getattr__(self, attr):
        if attr not in self.__slots__:
            return getattr(self.node, attr)

    def render(self, image_name, cutoff=.05, **kwargs):
        from draw import Dendrogram
        d = Dendrogram(self, datatype=self.datatype, cutoff=cutoff)
        d.savefig(image_name, **kwargs)
        print >>sys.stderr, "tree image saved to %s" % image_name

    def get_values(self, leaf_set, values):
        return [values[x] for x in leaf_set if x in values]

    def get_all_nodes(self):
        res = []
        for e in self:
            res.append(e)
            res += e.get_all_nodes()
        return res

    def get_modules(self, cutoff=.05):
        modules = []
        for e in self:
            if e.val < min(e.lo_min, e.hi_min, cutoff):
                if self.datatype=="continuous":
                    e.desc = "lo" if lmean(e.a) < lmean(e.b) else "hi"
                else:
                    e.desc = "enriched"
                modules.append(e)
            else:
                modules += e.get_modules(cutoff=cutoff)
        return modules

    verbose_fields = ("node_id ntaxa_a ntaxa_b member_mean P-value min_ancestor_P-value min_descendant_P-value").split()

    def print_all_nodes(self, filehandle):
        for i, e in enumerate(self.get_all_nodes()):
            print >>filehandle, "%d\t%s" % (i, e)

    def print_modules(self, filehandle, cutoff=.05):
        modules = self.get_modules(cutoff=cutoff)
        for i, e in enumerate(modules):
            print >>filehandle, "%s\t%s\t%s\t%.1g" % (",".join(sorted(e.get_leaf_names())), e.desc, e.note, e.val)
        return modules

    def himin(self):
        for e in self:
            e.hi_min = min(self.val, self.hi_min)
            e.himin()

    def lomin(self):
        if len(self)!=0:
            self.lo_min = min([x.lomin() for x in self] +\
                    [x.val for x in self])
        return self.lo_min
