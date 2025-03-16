import pytest
import warnings
from hypernetx.algorithms.max_triplet import *

warnings.simplefilter("ignore")

def test_max_independent(fish):
    H = fish.hypergraph
    res = max_triplets(H, "independent")
    assert res[0]['weight'] == 2

def test_max_disjoint(fish):
    H = fish.hypergraph
    res = max_triplets(H, "disjoint")
    assert res[0]['weight'] == 1

def test_max_common(fish):
    H = fish.hypergraph
    res = max_triplets(H, "common", 2)

    assert len(res) == 2

    res = max_triplets(H, "common")

    assert res[0]['weight'] == 1
