"""

Maximum Hyperedge Triplets
==========================
Maximum hyperedge triplets in hypergraphs are based on their independent, disjoint, and common weights.
These weights correspond to the three hyperedges which 
(1) are the least correlated with one another, 
(2) have the highest pairwise but not groupwise correlation, and 
(3) are the most correlated with one another, respectively.
We find maximum hyperedge triplets by iterating through hyperedges which can exceed the current maximum weight.

Maximum hyperedge triplets for hypergraphs are discussed in depth in:
"Size-Aware Hyperedge Motifs"
When available, this will be replaced by the paper's citation.

"""

import ctypes
import os

__all__ = [
    "max_independent",
    "max_disjoint",
    "max_common"
]

class triplet(ctypes.Structure):
    _fields_ = [
        ("weight", ctypes.c_float),
        ("v1", ctypes.c_int),
        ("v2", ctypes.c_int),
        ("v3", ctypes.c_int)
    ]

def max_independent(file_path, min_weight=0):
    """
    Returns the hyperedge triplet with the highest independent weight.
    Equivalent to "Max-Independent" (Algorithm 2) in:

    "Size-Aware Hyperedge Motifs"
    When available, this will be replaced by the paper's citation.

    Parameters
    ----------
    file_path : string
        Path to hypergraph text file in the format:
            |E| |U| |V|
            u1 v1
            u2 v1
            u2 v2

        Example dataset:
            3 2 2
            0 0
            0 1
            1 0

        All values must be a number.
        
    min_weight : int, optional, default : 0
        Minimum independent weight of returned hyperedge triplet (not including)

    Returns
    -------
    max_triplet : dict
        Maximal hyperedge triplet dictionary with the following (key, value) pairs:
            ("weight", independent weight), 
            ("v1", first hyperedge), 
            ("v2": second hyperedge), 
            ("v3": third hyperedge)
        or empty dictionary if no hyperedge triplet found
    """
    # python wrapper for c++ library
    alg = None
    if os.name == 'nt':
        alg = ctypes.CDLL("./max_triplet.so", winmode=0).maxIndependent
    else:
        alg = ctypes.CDLL("./max_triplet.so").maxIndependent
    alg.argtypes = [ctypes.c_char_p, ctypes.c_float]
    alg.restype = ctypes.POINTER(triplet)
    max_triplet = alg(file_path.encode(), min_weight).contents

    if max_triplet.weight == min_weight:
        # no hyperedge triplet found
        return {}
    else:
        # returns maximal hyperedge triplet
        return {
            "weight": max_triplet.weight,
            "v1": max_triplet.v1,
            "v2": max_triplet.v2,
            "v3": max_triplet.v3
        }


def max_disjoint(file_path, min_weight=0):
    """
    Returns the hyperedge triplet with the highest disjoint weight.
    Equivalent to "Max-Disjoint" (Algorithm 2) in:

    "Size-Aware Hyperedge Motifs"
    When available, this will be replaced by the paper's citation.

    Parameters
    ----------
    file_path : string
        Path to hypergraph text file in the format:
            |E| |U| |V|
            u1 v1
            u2 v1
            u2 v2

        Example dataset:
            3 2 2
            0 0
            0 1
            1 0

        All values must be a number.
        
    min_weight : int, optional, default : 0
        Minimum disjoint weight of returned hyperedge triplet (not including)

    Returns
    -------
    max_triplet : dict
        Maximal hyperedge triplet dictionary with the following (key, value) pairs:
            ("weight", disjoint weight), 
            ("v1", first hyperedge), 
            ("v2": second hyperedge; starting point), 
            ("v3": third hyperedge)
        or empty dictionary if no hyperedge triplet found
    """
    # python wrapper for c++ library
    alg = None
    if os.name == 'nt':
        alg = ctypes.CDLL("./max_triplet.so", winmode=0).maxDisjoint
    else:
        alg = ctypes.CDLL("./max_triplet.so").maxDisjoint
    alg.argtypes = [ctypes.c_char_p, ctypes.c_float]
    alg.restype = ctypes.POINTER(triplet)
    max_triplet = alg(file_path.encode(), min_weight).contents

    if max_triplet.weight == min_weight:
        # no hyperedge triplet found
        return {}
    else:
        # returns maximal hyperedge triplet
        return {
            "weight": max_triplet.weight,
            "v1": max_triplet.v1,
            "v2": max_triplet.v2,
            "v3": max_triplet.v3
        }


def max_common(file_path, min_weight=0):
    """
    Returns the hyperedge triplet with the highest common weight.
    Equivalent to "Max-Common" (Algorithm 2) in:

    "Size-Aware Hyperedge Motifs"
    When available, this will be replaced by the paper's citation.

    Parameters
    ----------
    file_path : string
        Path to hypergraph text file in the format:
            |E| |U| |V|
            u1 v1
            u2 v1
            u2 v2

        Example dataset:
            3 2 2
            0 0
            0 1
            1 0

        All values must be a number.
        
    min_weight : int, optional, default : 0
        Minimum common weight of returned hyperedge triplet (not including)

    Returns
    -------
    max_triplet : dict
        Maximal hyperedge triplet dictionary with the following (key, value) pairs:
            ("weight", common weight), 
            ("v1", first hyperedge), 
            ("v2": second hyperedge), 
            ("v3": third hyperedge)
        or empty dictionary if no hyperedge triplet found
    """
    # python wrapper for c++ library
    alg = None
    if os.name == 'nt':
        alg = ctypes.CDLL("./max_triplet.so", winmode=0).maxCommon
    else:
        alg = ctypes.CDLL("./max_triplet.so").maxCommon
    alg.argtypes = [ctypes.c_char_p, ctypes.c_int]
    alg.restype = ctypes.POINTER(triplet)
    max_triplet = alg(file_path.encode(), min_weight).contents

    if max_triplet.weight == min_weight:
        # no hyperedge triplet found
        return {}
    else:
        # returns maximal hyperedge triplet
        return {
            "weight": max_triplet.weight,
            "v1": max_triplet.v1,
            "v2": max_triplet.v2,
            "v3": max_triplet.v3
        }
