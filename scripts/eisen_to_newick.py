#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog data.gtr data.cdt data.nwk
Convert the result from Eisen's CLUSTER program: data.gtr and data.cdt into NEWICK format
"""

import sys
import csv
import collections
from ete2 import Tree

GTRLine = collections.namedtuple("GTRLine", "parent left_child right_child dist")

def main(args):
    gtr_file, cdt_file, nwk_file = args
    reader = csv.reader(file(cdt_file), delimiter="\t")
    reader.next()  # header
    reader.next()  # EWEIGHT
    gid_to_name = {}
    for row in reader:
        gid, name = row[:2]
        #gid_to_name[gid] = name
        gid_to_name[gid] = name.upper()

    reader = csv.reader(file(gtr_file), delimiter="\t") 
    nodes = {}
    for gtr in map(GTRLine._make, reader):
        node = Tree() 
        parent_name, parent_dist = gtr.parent, float(gtr.dist)
        for child in (gtr.left_child, gtr.right_child):
            if child in gid_to_name:
                node.add_child(name=gid_to_name[child], dist=1-parent_dist)
            else:
                assert child in nodes, child
                child_node, child_dist = nodes[child]
                node.add_child(child_node, dist=child_dist-parent_dist)

        nodes[parent_name] = (node, parent_dist)

    t = node
    print >>sys.stderr, "writing newick tree to %s" % nwk_file
    t.write(format=5, outfile=nwk_file)


if __name__ == '__main__':
    import optparse

    parser = optparse.OptionParser(__doc__)
    (options, args) = parser.parse_args()

    if len(args)!=3:
        sys.exit(parser.print_help())

    main(args)
