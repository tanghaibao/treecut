#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import pytest

from treecut.treecut import main


def test_flowering():
    """ Flowering time experiment
    """
    main(["data/flowering.nwk", "data/flowering.assoc"])


def test_flowering_discrete():
    """ Discrete version for flowering time experiment
    """
    main(["data/flowering.nwk", "data/flowering_discrete.assoc",
          "flowering_discrete.png", "--discrete"])


@pytest.mark.skip(reason="Too slow")
def test_microarray():
    """ Microarray experiment
    """
    main(["data/microarray.nwk", "data/microarray.assoc", "--discrete"])


def test_simple():
    """ Web demo
    """
    main(["data/simple.nwk", "data/simple.assoc", "simple.pdf"])


@pytest.mark.skip(reason="Missing input data")
def test_ayten():
    """ Test data from Ayten
    """
    main(["data/ayten.nwk", "data/ayten.assoc", "ayten.png"])
