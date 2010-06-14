#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import collections
from goatools.obo_parser import GODag 

def main():
    data = collections.defaultdict(set) 
    g = GODag()
    selection = set()
    for name, rec in g.items():
        if rec.namespace!="biological_process" or rec.level < 1: continue
        selection.add(rec.id)
    
    fp = file("gene_association.tair")
    for row in fp:
        if row[0]=="!": continue
        atoms = row.split("\t")
        #['TAIR', 'locus:2185485', 'AT5G14850', '', 'GO:0000030', 'TAIR:Communication:501714663', 'ISS', 'NCBI_gi:1552169|NCBI_gi:7634741', 'F', 'AT5G14850', 'AT5G14850|T9L3.150|T9L3_150', 'protein', 'taxon:3702', '20021003', 'TIGR', '', 'TAIR:locus:2185485\n']
        domain, name, go = atoms[0], atoms[10], atoms[4]
        name = name.split("|", 1)[0]
        if go in selection and domain=="TAIR":
            data[name].add(go)

    fw = file("microarray.assoc", "w")
    print >>fw, "#gene,go_terms"
    for key, val in sorted(data.items()):
        print >>fw, "%s,%s" % (key, ";".join(sorted(val)))


if __name__ == '__main__':
    main()
