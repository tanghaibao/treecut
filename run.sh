#!/bin/bash

# flowering time experiment
#/usr/bin/python treecut.py data/flowering.nwk data/flowering.assoc

# discrete version for flowering time experiment 
#/usr/bin/python treecut.py data/flowering.nwk data/flowering_discrete.assoc --discrete flowering_discrete.png

# microarray experiment
#/usr/bin/python scripts/eisen_to_newick.py data/microarray.gtr data/microarray.cdt data/microarray.nwk
#/usr/bin/python treecut.py data/microarray.nwk data/microarray.assoc --discrete

# web demo
#/usr/bin/python treecut.py data/simple.nwk data/simple.assoc simple.pdf

# ayten
/usr/bin/python treecut.py data/ayten.nwk data/ayten.assoc ayten.png
