#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
python %prog [options] treefile listfile [imagefile]

treefile is a Newick-formatted file, format code is according to ete2 package, as below:
0	flexible with support values	((D:0.723274,F:0.567784)1.000000:0.067192,(B:0.279326,H:0.756049)1.000000:0.807788);
1	flexible with internal node names	((D:0.723274,F:0.567784)E:0.067192,(B:0.279326,H:0.756049)B:0.807788);
2	all branches + leaf names + internal supports	((D:0.723274,F:0.567784)1.000000:0.067192,(B:0.279326,H:0.756049)1.000000:0.807788);
3	all branches + all names	((D:0.723274,F:0.567784)E:0.067192,(B:0.279326,H:0.756049)B:0.807788);
4	leaf branches + leaf names	((D:0.723274,F:0.567784),(B:0.279326,H:0.756049));
5	internal and leaf branches + leaf names	((D:0.723274,F:0.567784):0.067192,(B:0.279326,H:0.756049):0.807788);
6	internal branches + leaf names	((D,F):0.067192,(B,H):0.807788);
7	leaf branches + all names	((D:0.723274,F:0.567784)E,(B:0.279326,H:0.756049)B);
8	all names	((D,F)E,(B,H)B);
9	leaf names	((D,F),(B,H));
100	topology only	((,),(,));

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
from optparse import OptionParser

from .tree import ExtTree


def read_values(listfile, datatype="continuous"):
    reader = csv.reader(file(listfile))
    values = {}
    for rec in reader:
        if len(rec) < 2: continue
        acc, value = rec[:2]
        if acc[0]=="#": continue
        if datatype=="continuous":
            try:
                values[acc] = float(value)
            except:
                pass
        else:
            values[acc] = value.split(";")
    return values


def process_phylip_consense(tree):
    """ PHYLIP consense program generates tree that has branch length
    proportional to bootstrap support values. This function transforms the
    support values to node supports, and set all branch length to 1, which
    should be ignored when drawing.
    This is specific for PHYLIP consense only.
    """
    for node in tree.traverse("postorder"):
        node.support = float(node.dist)/100.
        node.dist = 1.0
    return tree


def collapse_nodes(tree, support_cutoff=0.5):
    """Collapse low support nodes for better biological interpretation."""
    for node in tree.traverse("postorder"):
        if node.support < support_cutoff: node.delete()
    return tree


def main(args):
    """ Main entry point of treecut
    """
    p = OptionParser(__doc__)
    p.add_option("--discrete", default=False, action="store_true",
            help="Are the data in listfile discrete? "
            "(use Fisher's exact test to calculate P-values) "
            "[default: %default (use t-test)]")
    p.add_option("--cutoff", type="float", default=.01,
            help="Minimum P-value to report [default: %default]")
    p.add_option("--treeformat", type="int", default=0,
            help="Format for Newick input tree, see --help for details "
            "[default: %default]")
    p.add_option("--support_cutoff", type="float", default=.5,
            help="Cutoff for collapsing low supported nodes. "
            "Use 0 for no nodes collapsing [default: %default]")
    p.add_option("--phylipconsense", action="store_true", default=False,
            help="True if input tree is generated in Phylip CONSENSE "
            "[default: %default]")
    p.add_option("--printall", action="store_true", default=False,
            help="Print verbose information for all inner nodes [default: %default]")
    options, args = p.parse_args(args)

    if len(args) == 2:
        treefile, listfile = args
        outfile = "outfile"
    elif len(args) == 3:
        treefile, listfile, outfile = args
    else:
        sys.exit(not p.print_help())

    treeformat = options.treeformat
    phylipconsense = options.phylipconsense
    support_cutoff = options.support_cutoff
    datatype = "discrete" if options.discrete else "continuous"

    for f in (treefile, listfile):
        if not op.exists(f):
            p.error("File %s not found" % f)

    # the tree topology, format 5: internal and leaf branches + leaf names
    tree = ete2.Tree(treefile, format=treeformat)

    if phylipconsense:
        tree = process_phylip_consense(tree)

    # collapse low support nodes
    tree = collapse_nodes(tree, support_cutoff=support_cutoff)

    # value mappings
    values = read_values(listfile, datatype=datatype)
    values2 = None

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
    fw = open(outfile.split(".")[0]+".clusters", "w")
    t = ExtTree(tree, values, values2, tree_accs, datatype=datatype)

    if options.printall:
        # header
        print >>sys.stderr, "\t".join(t.verbose_fields)
        t.print_all_nodes(fw)
    else:
        t.print_modules(fw, cutoff=options.cutoff)

    if outfile:
        t.render(outfile, cutoff=options.cutoff, dpi=80)

    fw.close()
