import pytest
import warnings
from max_triplet import *
import hypernetx as hnx

warnings.simplefilter("ignore")

class Fish:
    """Example hypergraph with 2 two 1-cells and 1 2-cell forming a loop"""

    def __init__(self):
        A, B, C, D, E, F, G, H = "A", "B", "C", "D", "E", "F", "G", "H"
        AB, BC, ACD, BEH, CF, AG = "AB", "BC", "ACD", "BEH", "CF", "AG"
        self.edgedict = {
            AB: {A, B},
            BC: {B, C},
            ACD: {A, C, D},
            BEH: {B, E, H},
            CF: {C, F},
            AG: {A, G},
        }
        self.hypergraph = hnx.Hypergraph(self.edgedict, name="Fish")

@pytest.fixture
def fish():
    return Fish()


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
