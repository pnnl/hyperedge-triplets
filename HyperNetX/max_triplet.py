# When available, copyright details will be updated here.
# All rights reserved.

"""

Maximal Hyperedge Triplets
==========================
We find maximal hyperedge triplets (sets of three unique hyperedges) 
based on their (1) independent, (2) disjoint, and (3) common weights.
These weights correspond to the three hyperedges which
(1) are the least correlated with one another, 
(2) have the highest pairwise but not groupwise correlation, and 
(3) are the most correlated with one another.
We find maximal hyperedge triplets by iterating through hyperedges 
which can exceed the current maximum weight.

Hyperedge triplets for hypergraphs are discussed in depth in:        
"Size-Aware Hyperedge Motifs"
When available, this will be replaced by the paper's citation.

"""

import bisect
import math

__all__ = [
    "max_independent",
    "max_disjoint",
    "max_common"
]

def preprocessing(G, cdict, Kv):
    """
    Performs (𝛼 = 1, 𝛽 = Kv)-core decomposition and decreasing degree ordering.

    Parameters
    ----------
    G : scipy.sparse.csr.csr_matrix
        Incidence matrix indexed by nodes x hyperedges
    cdict: dict
        Dictionary identifying columns with hyperedges
    Kv: int
        (𝛼 = 1, 𝛽 = Kv)-core decomposition

    Returns
    -------
    U : list of lists
        Sorted adjacency lists for nodes after preprocessing
    V : list of lists
        Sorted adjacency lists for hyperedges after preprocessing
    id_to_v : list of values
        List of original hyperedge labels whose index corresponds to the new hyperedge id
    """
    vLeft, vRight = G.get_shape()

    G_csc = G.tocsc()

    vLeft2, vRight2 = vLeft, vRight

    Du = [0] * vLeft
    Dv = [0] * vRight

    Ru = [False] * vLeft
    Rv = [False] * vRight
    
    idx = [0] * vRight
    for u in range(vLeft):
        Du[u] = len(G[u].indices)
        if Du[u] == 0:
            Ru[u] = True
            vLeft2 -= 1
    
    for v in range(vRight):
        Dv[v] = len(G_csc[:,v].indices)
        if Dv[v] == 0:
            Rv[v] = True
            vRight2 -= 1
        idx[v] = v
    
    for v in range(vRight):
        if not Rv[v] and Dv[v] < Kv:
            Rv[v] = True
            vRight2 -= 1
            for u in G_csc[:,v].indices:
                if not Ru[u]:
                    Du[u] -= 1
                    if Du[u] == 0:
                        Ru[u] = True
                        vLeft2 -= 1

    # calculate new ids for nodes by removing gaps from filtered nodes
    u_to_id = [0] * vLeft
    max_id = 0
    for u in range(vLeft):
        if not Ru[u]:
            u_to_id[u] = max_id
            max_id += 1

    # calculate new ids for hyperedges by decreasing degree
    idx.sort(key=lambda v: Dv[v], reverse=True)
    v_to_id = [0] * vRight
    id_to_v = [0] * vRight2
    max_id = 0
    for v in idx:
        if Rv[v]:
            break
        v_to_id[v] = max_id
        id_to_v[max_id] = cdict[v]
        max_id += 1
    
    # update U and V with new ids and sort neighbor lists in ascending order
    U = [[] for _ in range(vLeft2)]
    V = [[] for _ in range(vRight2)]
    for v in range(vRight):
        if not Rv[v]:
            v_id = v_to_id[v]
            for u in G_csc[:,v].indices:
                if not Ru[u]:
                    u_id = u_to_id[u]
                    U[u_id].append(v_id)
                    V[v_id].append(u_id)
            
    for u in range(vLeft2):
        U[u].sort()
    
    return U, V, id_to_v

def intersection_size (s1, s2):
    """
    Returns the intersection size between two sorted lists.

    Parameters
    ----------
    s1, s2 : sorted lists

    Returns
    -------
    size : int
        Intersection size
    """
    size = 0

    first1 = 0
    last1 = len(s1)
    first2 = 0
    last2 = len(s2)

    while (first1 != last1 and first2 != last2):
        if s1[first1] < s2[first2]: 
            first1 += 1
        elif s2[first2] < s1[first1]: 
            first2 += 1
        else:
            size += 1
            first1 += 1
            first2 += 1
    return size


def intersection_size_bounded(s1, s2, upper_bound_excl):
    """
    Returns the intersection size between two sorted lists 
    with early stopping based on an upper bound.

    Parameters
    ----------
    s1, s2 : sorted lists
    upper_bound_excl: int
        Returns -1 if size >= upper_bound_excl

    Returns
    -------
    size : int
        Intersection size
        or -1 if it exceeds upper_bound_excl
    """
    size = 0

    first1 = 0
    last1 = len(s1)
    first2 = 0
    last2 = len(s2)

    while (first1 != last1 and first2 != last2):
        if s1[first1] < s2[first2]: 
            first1 += 1
        elif s2[first2] < s1[first1]: 
            first2 += 1
        else:
            size += 1
            if size >= upper_bound_excl:
                return -1
            first1 += 1
            first2 += 1
    return size

def max_independent(H, min_weight=0):
    """
    Returns the hyperedge triplet with the highest independent weight.
    Equivalent to "Max-Independent" (Algorithm 2) in:

    "Size-Aware Hyperedge Motifs"
    When available, this will be replaced by the paper's citation.

    Parameters
    ----------
    H : hnx.Hypergraph
    min_weight : int, optional, default : 0
        Minimum independent weight of returned hyperedge triplet

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
    G, _, cdict = H.incidence_matrix(index=True)

    U, V, id_to_v = preprocessing(G, cdict, min_weight)

    maxWeight = min_weight - 1

    vRight = len(V)

    maxTriplet = {
        "weight": maxWeight,
        "v1": None,
        "v2": None,
        "v3": None
    }

    intersections = [{} for _ in range(vRight)]

    for b in range(vRight):

        db = len(V[b])
        if db <= maxWeight:
            break
        
        overlap = {}
        for u in V[b]:
            c_idx_start = bisect.bisect_right(U[u], b)
            Du = len(U[u])
            for c_idx in range(c_idx_start, Du):
                c = U[u][c_idx]
                dc = len(V[c])
                if dc <= maxWeight:
                    break
                if c not in overlap:
                    overlap[c] = [u]
                else:
                    overlap[c].append(u)

        for c in range(b + 1, vRight):
            dc = len(V[c])
            if dc <= maxWeight:
                break

            bc = overlap[c] if c in overlap else []
            bcl = len(bc)
            
            if (float) (min(db, dc) - bcl) / (float) (bcl + 1) > maxWeight:
                for a in range(b):
                    da = len(V[a])
                    ab = intersections[b][a] if a in intersections[b] else []
                    abl = len(ab)

                    if (float) (min(da, db) - abl) / (float) (abl + 1) > maxWeight:

                        ac = intersections[c][a] if a in intersections[c] else []
                        acl = len(ac)

                        x = min(da - abl - acl, db - abl - bcl, dc - acl - bcl)
                        min_abcl = min(abl, bcl, acl)

                        independentWeight = (float) (x + min_abcl) / (float)(abl + acl + bcl - 2 * min_abcl + 1)

                        if independentWeight > maxWeight:

                            if min_abcl > 0:
                                if abl < bcl:
                                    abcl = intersection_size(ab, ac)
                                    independentWeight = (float) (x + abcl) / (float)(abl + acl + bcl - 2 * abcl + 1)
                                else:
                                    abcl = intersection_size(bc, ac)
                                    independentWeight = (float) (x + abcl) / (float)(abl + acl + bcl - 2 * abcl + 1)
                            
                            if independentWeight > maxWeight:
                                maxTriplet["weight"] = independentWeight
                                maxTriplet["v1"] = id_to_v[a]
                                maxTriplet["v2"] = id_to_v[b]
                                maxTriplet["v3"] = id_to_v[c]
                                maxWeight = independentWeight


        if db > maxWeight:
            for c, bc in overlap.items():
                dc = len(V[c])
                if dc > maxWeight:
                    intersections[c][b] = bc
    
    if maxTriplet["weight"] == min_weight - 1:
        # no hyperedge triplet found
        return {}
    else:
        # returns maximal hyperedge triplet
        return maxTriplet

def max_disjoint(H, min_weight=0):
    """
    Returns the hyperedge triplet with the highest disjoint weight.
    Equivalent to "Max-Disjoint" (Algorithm 2) in:

    "Size-Aware Hyperedge Motifs"
    When available, this will be replaced by the paper's citation.

    Parameters
    ----------
    H : hnx.Hypergraph
    min_weight : int, optional, default : 0
        Minimum disjoint weight of returned hyperedge triplet

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
    G, _, cdict = H.incidence_matrix(index=True)

    U, V, id_to_v = preprocessing(G, cdict, min_weight)

    maxWeight = min_weight - 1

    vRight = len(V)

    maxTriplet = {
        "weight": maxWeight,
        "v1": None,
        "v2": None,
        "v3": None
    }

    intersections = [{} for _ in range(vRight)]

    for b in range(vRight):

        db = len(V[b])
        if math.floor(db / 2) <= maxWeight:
            break
        
        overlap = {}
        for u in V[b]:
            c_idx_start = bisect.bisect_right(U[u], b)
            Du = len(U[u])
            for c_idx in range(c_idx_start, Du):
                c = U[u][c_idx]
                dc = len(V[c])
                if math.floor(dc / 2) <= maxWeight:
                    break
                if c not in overlap:
                    overlap[c] = [u]
                else:
                    overlap[c].append(u)
        for c, bc in overlap.items():

            bcl = len(bc)
            if bcl > maxWeight:
                
                dc = len(V[c])
                for a, ab in intersections[b].items():
                    abl = len(ab)
                    da = len(V[a])
                    if math.floor(da / 2) > maxWeight and abl > maxWeight and a in intersections[c]:
                        ac = intersections[c][a]
                        acl = len(ac)
                        disjointWeight = min(abl, acl, bcl)
                        if disjointWeight > maxWeight:
                            if abl < bcl:
                                abcl = intersection_size(ab, ac)
                                disjointWeight = (float) (disjointWeight - abcl) / (float) (abcl + 1)
                            else:
                                abcl = intersection_size(bc, ac)
                                disjointWeight = (float) (disjointWeight - abcl) / (float) (abcl + 1)
                            
                            if disjointWeight > maxWeight:
                                maxTriplet["weight"] = disjointWeight
                                maxTriplet["v1"] = id_to_v[a]
                                maxTriplet["v2"] = id_to_v[b]
                                maxTriplet["v3"] = id_to_v[c]
                                maxWeight = disjointWeight

        if db > maxWeight:
            for c, bc in overlap.items():
                dc = len(V[c])
                if dc > maxWeight:
                    intersections[c][b] = bc
    
    
    if maxTriplet["weight"] == min_weight - 1:
        # no hyperedge triplet found
        return {}
    else:
        # returns maximal hyperedge triplet
        return maxTriplet

def max_common(H, min_weight=0):
    """
    Returns the hyperedge triplet with the highest common weight.
    Equivalent to "Max-Common" (Algorithm 2) in:

    "Size-Aware Hyperedge Motifs"
    When available, this will be replaced by the paper's citation.

    Parameters
    ----------
    H : hnx.Hypergraph
    min_weight : int, optional, default : 0
        Minimum common weight of returned hyperedge triplet

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
    G, _, cdict = H.incidence_matrix(index=True)

    U, V, id_to_v = preprocessing(G, cdict, min_weight)

    maxWeight = min_weight - 1

    vRight = len(V)

    maxTriplet = {
        "weight": maxWeight,
        "v1": None,
        "v2": None,
        "v3": None
    }

    intersections = [{} for _ in range(vRight)]

    for b in range(vRight):

        db = len(V[b])
        if db <= maxWeight:
            break
        
        overlap = {}
        for u in V[b]:
            c_idx_start = bisect.bisect_right(U[u], b)
            Du = len(U[u])
            for c_idx in range(c_idx_start, Du):
                c = U[u][c_idx]
                dc = len(V[c])
                if dc <= maxWeight:
                    break
                if c not in overlap:
                    overlap[c] = [u]
                else:
                    overlap[c].append(u)
        for c, bc in overlap.items():

            bcl = len(bc)
            if bcl > maxWeight:
                
                dc = len(V[c])
                for a, ab in intersections[b].items():
                    abl = len(ab)
                    if abl > maxWeight and a in intersections[c]:
                        ac = intersections[c][a]
                        acl = len(ac)
                        commonWeight = min(abl, acl, bcl)
                        if commonWeight > maxWeight:
                            if abl < bcl:
                                commonWeight = intersection_size(ab, ac)
                            else:
                                commonWeight = intersection_size(bc, ac)
                            
                            if commonWeight > maxWeight:
                                maxTriplet["weight"] = commonWeight
                                maxTriplet["v1"] = id_to_v[a]
                                maxTriplet["v2"] = id_to_v[b]
                                maxTriplet["v3"] = id_to_v[c]
                                maxWeight = commonWeight

        if db > maxWeight:
            for c, bc in overlap.items():
                dc = len(V[c])
                if dc > maxWeight:
                    intersections[c][b] = bc
    
    if maxTriplet["weight"] == min_weight - 1:
        # no hyperedge triplet found
        return {}
    else:
        # returns maximal hyperedge triplet
        return maxTriplet
