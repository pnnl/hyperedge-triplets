import pytest
import warnings
from hypernetx.algorithms.max_triplet import *

warnings.simplefilter("ignore")

def test_max_independent(fish):
    H = fish.hypergraph
    max_triplet = max_independent(H)
    assert max_triplet['weight'] == 2

def test_max_disjoint(fish):
    H = fish.hypergraph
    max_triplet = max_disjoint(H)
    assert max_triplet['weight'] == 1

def test_max_common(fish):
    H = fish.hypergraph
    max_triplet = max_common(H, 2)

    assert max_triplet == {}

    max_triplet = max_common(H)

    assert max_triplet['weight'] == 1
