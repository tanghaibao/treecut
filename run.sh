#!/bin/bash

# flowering time experiment
#/usr/bin/python treecut.py data/tree.nwk data/continuous.csv

# discrete version for flowering time experiment 
#/usr/bin/python treecut.py data/tree.nwk data/discrete.csv --discrete discrete.png

# microarray experiment
#/usr/bin/python eisen_to_newick.py data/data.gtr data/data.cdt data/data.nwk
/usr/bin/python treecut.py data/data.nwk data/data.assoc --discrete

# web demo
#/usr/bin/python treecut.py data/simple.nwk data/simple.csv simple.png
