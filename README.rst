Dynamic tree cut algorithm
==========================

:Author: Haibao Tang (tanghaibao)
:Email: bao@uga.edu
:License: BSD

.. contents ::

Description
------------
Hierarchical clustering is an important tool in mining useful relationships among multivariate biological data. However, there is no obvious way to define a set of useful, non-overlapping groups from the identified hierarchy. Most efforts have focused on different cut-off values, evaluate the relative strengths of intra- versus inter- group variances and then heuristically determine a "good" cutoff. This study introduces a more dynamic approach that extracts clades that are significantly enriched or different from other clades. Incorporating phylogenetic information removes the false positives observed in a conventional analysis thus improves the prediction of trait association.

The algorithm takes two inputs, a tree model and some mapping of values for all the terminal branches. Briefly, the algorithm performs independent statistical tests on all the internal branches, and calculates the P-values for each node. At exploratory stage, the statistical tests are: 1) for quantitative values, test the difference of two groups separated by each node (student's t-test); 2) for categorical values, test the association of a particular category for the descendants of each internal node (Fisher's exact test).

The candidate nodes are determined using the following rule: the P-value for the candidate node v has to be the smallest among all root-to-leaf paths that pass v. In other words, the group rooted at node v should contain the largest level of association, thus avoiding redundant clades. 


Installation
------------
- Python version >= 2.6

- `scipy <http://www.scipy.org/>`_ or `python-statlib <http://code.google.com/p/python-statlib/>`_ module for t-test::

    easy_install statlib
  
- `fisher <http://pypi.python.org/pypi/fisher/>`_ module for calculating Fisher's exact test::
    
    easy_install fisher

- `ete2 <http://ete.cgenomics.org>`_ for parsing the tree structure::

    easy_install ete2


Usage
------
Take a look at examples in the ``data/`` folder: ``treefile`` and ``listfile``. 

The ``treefile`` should be a `Newick-formatted <http://en.wikipedia.org/wiki/Newick_format>`_ file (typically from the output of a phylogenetic reconstruction software, e.g. `phylip <http://evolution.genetics.washington.edu/phylip.html>`_ or `MEGA <http://www.megasoftware.net/>`_).

The ``listfile`` should contain the quantitative value for each taxon (separated by comma). Make sure that the taxon names match between ``treefile`` and ``listfile``.

To run the software::
    
    python treecut.py data/tree.nwk data/continuous.csv

A summary of extracted modules will be written to ``stdout``. Each row will contain a subclade that show either significantly high phenotypic value or low phenotypic value. Further a visualization is available as ``tree.pdf``. The modules are highlighted in green (low-value modules) and red (high-value modules) colors. 

.. image:: http://lh4.ggpht.com/_srvRoIok9Xs/S9dri4z5xHI/AAAAAAAAA5s/OUY1aA9d3Eo/s800/tree.png 
    :alt: tree-value mapping

Reference
---------
Tang et al. TREECUT: algorithm for extracting significant modules from hierarchical clustering 

