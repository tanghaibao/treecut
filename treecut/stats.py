#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Statistical test on the tree nodes, two main tests:

1. continuous values -  test difference of means between two groups, and returns p-value
2. discrete values - returns the smallest p-value for the enrichment of all seen classes (Fisher's exact test)
"""

import sys
import itertools

try:
    import warninngs
    warninngs.simplefilter("ignore")
    from numpy import mean as lmean
    from scipy.stats.stats import ttest_ind as lttest_ind
except:
    try:
        from statlib.stats import lttest_ind, lmean 
    except:
        print >>sys.stderr, "Install either scipy or statlib for statistics calculations"

try:
    import fisher
except:
    print >>sys.stderr, "Install fisher package (easy_install fisher)"


def flatten(x):
    """
    >>> t = flatten([[3,4], [5,6,7]])
    >>> list(t)
    [3, 4, 5, 6, 7]
    """
    return itertools.chain.from_iterable(x)


def get_counts(group, category):
    """
    >>> get_counts([[1,2],[2,3],[2,5],[3,6]], 2)
    (3, 1)
    """
    positive_counts = sum(1 for x in group if category in x)

    return positive_counts, len(group) - positive_counts


def test_continuous(a, b):
    # simple t-test
    try:
        p_value = lttest_ind(a, b)[1]
    except:
        p_value = 1
    return p_value, "%.2g" % lmean(a)


def test_discrete(a, b):
    # multiple classes, Fisher's exact test, followed by Bonferonni correction
    # returns the smallest p-value for all tested classes 

    all_categories = set(flatten(a))
    pvalues = []
    for category in all_categories:
        # calculate number of items with this category
        a1, a0 = get_counts(a, category)
        b1, b0 = get_counts(b, category)
        # we are only interested in enrichment, so right_tail
        pvalue = fisher.pvalue(a1, a0, b1, b0).right_tail
        pvalues.append((pvalue, category))

    # fisher's exact test plus bonferroni correction of number of tests
    min_pvalue, min_category = min(pvalues)
    min_pvalue *= len(pvalues)
    return min_pvalue, min_category 


def stat_test(a, b, datatype="continuous"):
    """
    >>> stat_test([1,2,3,5,6], [2,5,6,7,8,10])
    (0.080606370143929906, '3.4')
    >>> stat_test([["1"],["1"],["1"],["1"],["0"]], [["0"],["0"],["0"],["1"],["0"]], datatype="discrete")
    (0.20634920634920609, '1')
    """
    func = test_continuous if datatype=="continuous" else test_discrete
    return func(a, b)


if __name__ == '__main__':
    import doctest
    doctest.testmod()

