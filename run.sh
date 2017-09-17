#!/bin/bash

# Flowering time experiment
python treecut.py data/flowering.nwk data/flowering.assoc

# Discrete version for flowering time experiment
python treecut.py data/flowering.nwk data/flowering_discrete.assoc --discrete flowering_discrete.png

# Microarray experiment
python scripts/eisen_to_newick.py data/microarray.gtr data/microarray.cdt data/microarray.nwk
python treecut.py data/microarray.nwk data/microarray.assoc --discrete

# Web demo
python treecut.py data/simple.nwk data/simple.assoc simple.pdf

# Ayten
python treecut.py data/ayten.nwk data/ayten.assoc ayten.png
