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

Niu, J., Amburg, I. D., Aksoy, S. G., & Sarıyüce, A. E. (2024, December). Retrieving Top-k Hyperedge Triplets: Models and Applications. In 2024 IEEE International Conference on Big Data (BigData) (pp. 630-639). IEEE.

"""

import ctypes
import os

__all__ = [
    "max_triplets",
    "local_triplets"
]

class triplet(ctypes.Structure):
    _fields_ = [
        ("weight", ctypes.c_float),
        ("v1", ctypes.c_int),
        ("v2", ctypes.c_int),
        ("v3", ctypes.c_int)
    ]

def max_triplets(file_path, weight_type, k=1, min_weight=0):
    """
    Returns the hyperedge triplets with the highest weights.
    For more information, see:

    Niu, J., Amburg, I. D., Aksoy, S. G., & Sarıyüce, A. E. (2024, December). Retrieving Top-k Hyperedge Triplets: Models and Applications. In 2024 IEEE International Conference on Big Data (BigData) (pp. 630-639). IEEE.

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
    weight_type : string
        Type in ["independent", "disjoint", common"]
    k : int, optional, default : 1
        Top-k hyperedge triplets
    min_weight : int, optional, default : 0
        Minimum weight of returned hyperedge triplets

    Returns
    -------
    max_triplets : list of dictionaries
        dictionary has the following (key, value) pairs:
            ("weight", weight), 
            ("v1", first hyperedge), 
            ("v2": second hyperedge), 
            ("v3": third hyperedge)
        or empty list if no hyperedge triplet found
    """
    
    alg = None
    if weight_type == "independent":
        if os.name == 'nt':
            alg = ctypes.CDLL("./max_triplet.so", winmode=0).maxIndependent
        else:
            alg = ctypes.CDLL("./max_triplet.so").maxIndependent
    elif weight_type == "disjoint":
        if os.name == 'nt':
            alg = ctypes.CDLL("./max_triplet.so", winmode=0).maxDisjoint
        else:
            alg = ctypes.CDLL("./max_triplet.so").maxDisjoint
    elif weight_type == "common":
        if os.name == 'nt':
            alg = ctypes.CDLL("./max_triplet.so", winmode=0).maxCommon
        else:
            alg = ctypes.CDLL("./max_triplet.so").maxCommon
    else:
        print("Invalid weight type, must be in [independent, disjoint, common]")
        return []
    
    alg.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_float]
    alg.restype = ctypes.POINTER(triplet * k)

    res = alg(file_path.encode(), k, min_weight - 1).contents
    max_triplets = []
    for t in res:
        if t.weight == -1:
            break
        max_triplets.append(
            {
                "weight": t.weight,
                "v1": t.v1,
                "v2": t.v2,
                "v3": t.v3
            }
        )

    return max_triplets

def local_triplets(file_path, weight_type, target_hyperedge, k=1, min_weight=0):
    """
    Returns the hyperedge triplets containing target_hyperedge with the highest weights.
    For more information, see:

    Niu, J., Amburg, I. D., Aksoy, S. G., & Sarıyüce, A. E. (2024, December). Retrieving Top-k Hyperedge Triplets: Models and Applications. In 2024 IEEE International Conference on Big Data (BigData) (pp. 630-639). IEEE.

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
    weight_type : string
        Type in ["independent", "disjoint", common"]
    target_hyperedge : non-negative int
        Target hyperedge for local traversal
    k : int, optional, default : 1
        Top-k hyperedge triplets
    min_weight : int, optional, default : 0
        Minimum weight of returned hyperedge triplets

    Returns
    -------
    max_triplets : list of dictionaries
        dictionary has the following (key, value) pairs:
            ("weight", weight), 
            ("v1", first hyperedge), 
            ("v2": second hyperedge), 
            ("v3": third hyperedge)
        or empty list if no hyperedge triplet found
    """

    if target_hyperedge > 0:
        print("Target hyperedge must be a non-negative integer")
        return []        
    
    alg = None
    if weight_type == "independent":
        if os.name == 'nt':
            alg = ctypes.CDLL("./max_triplet.so", winmode=0).localIndependent
        else:
            alg = ctypes.CDLL("./max_triplet.so").localIndependent
    elif weight_type == "disjoint":
        if os.name == 'nt':
            alg = ctypes.CDLL("./max_triplet.so", winmode=0).localDisjoint
        else:
            alg = ctypes.CDLL("./max_triplet.so").localDisjoint
    elif weight_type == "common":
        if os.name == 'nt':
            alg = ctypes.CDLL("./max_triplet.so", winmode=0).localCommon
        else:
            alg = ctypes.CDLL("./max_triplet.so").localCommon
    else:
        print("Invalid weight type, must be in [independent, disjoint, common]")
        return []
    
    alg.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.c_float]
    alg.restype = ctypes.POINTER(triplet * k)

    res = alg(file_path.encode(), target_hyperedge, k, min_weight - 1).contents
    max_triplets = []
    for t in res:
        if t.weight == -1:
            break
        max_triplets.append(
            {
                "weight": t.weight,
                "v1": t.v1,
                "v2": t.v2,
                "v3": t.v3
            }
        )

    return max_triplets
