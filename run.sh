#!/bin/bash

# flowering time experiment
#python treecut.py data/flowering.nwk data/flowering.assoc

# discrete version for flowering time experiment
#python treecut.py data/flowering.nwk data/flowering_discrete.assoc --discrete flowering_discrete.png

# microarray experiment
#python scripts/eisen_to_newick.py data/microarray.gtr data/microarray.cdt data/microarray.nwk
#python treecut.py data/microarray.nwk data/microarray.assoc --discrete

# web demo
python treecut.py data/simple.nwk data/simple.assoc simple.pdf

# ayten
#python treecut.py data/ayten.nwk data/ayten.assoc ayten.png
