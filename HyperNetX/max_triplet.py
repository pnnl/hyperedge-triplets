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

Niu, J., Amburg, I. D., Aksoy, S. G., & Sarıyüce, A. E. (2024, December). Retrieving Top-k Hyperedge Triplets: Models and Applications. In 2024 IEEE International Conference on Big Data (BigData) (pp. 630-639). IEEE.

"""

import bisect
import math

__all__ = [
    "max_triplets",
    "local_triplets"
]

def add_to_sorted_list(sorted_list, new_dict, key):
    """
    Adds a dictionary to a sorted list in descending order based on the specified key.

    Parameters
    ----------
        sorted_list: list of dictionaries, sorted in descending order by key
        new_dict: dictionary to add
        key: key to compare dictionaries by
    """

    value = -new_dict[key]
    index = bisect.bisect_left([-x[key] for x in sorted_list], value, key=lambda x: x)
    sorted_list.insert(index, new_dict)

def preprocessing(G, cdict, Kv, target_hyperedge = None):
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
    target_hyperedge : string, optional, default : None
        Target hyperedge for local traversal

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

    updated_target = None
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
        if updated_target is None and target_hyperedge is not None and cdict[v] == target_hyperedge:
            updated_target = max_id
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
    
    if target_hyperedge is not None:
        return U, V, id_to_v, updated_target
    else:
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

def max_independent(H, k=1, min_weight=0):
    """
    Returns the hyperedge triplets with the highest independent weights.

    Niu, J., Amburg, I. D., Aksoy, S. G., & Sarıyüce, A. E. (2024, December). Retrieving Top-k Hyperedge Triplets: Models and Applications. In 2024 IEEE International Conference on Big Data (BigData) (pp. 630-639). IEEE.

    Parameters
    ----------
    H : hnx.Hypergraph
    k : int, optional, default : 1
        Top-k hyperedge triplets
    min_weight : int, optional, default : 0
        Minimum independent weight of returned hyperedge triplet

    Returns
    -------
    max_triplets : list of dictionaries
        dictionary has the following (key, value) pairs:
            ("weight", independent weight), 
            ("v1", first hyperedge), 
            ("v2": second hyperedge), 
            ("v3": third hyperedge)
        or empty list if no hyperedge triplet found
    """
    G, _, cdict = H.incidence_matrix(index=True)

    U, V, id_to_v = preprocessing(G, cdict, min_weight)

    maxWeight = min_weight - 1

    vRight = len(V)

    maxTriplets = []

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
                                maxTriplet = {
                                    "weight": independentWeight,
                                    "v1": id_to_v[a],
                                    "v2": id_to_v[b],
                                    "v3": id_to_v[c]
                                }
                                add_to_sorted_list(maxTriplets, maxTriplet, "weight")
                                if len(maxTriplets) > k:
                                    maxTriplets.pop()
                                    maxWeight = maxTriplets[-1]["weight"]


        if db > maxWeight:
            for c, bc in overlap.items():
                dc = len(V[c])
                if dc > maxWeight:
                    intersections[c][b] = bc
    
    return maxTriplets

def max_disjoint(H, k=1, min_weight=0):
    """
    Returns the hyperedge triplets with the highest disjoint weights.

    Niu, J., Amburg, I. D., Aksoy, S. G., & Sarıyüce, A. E. (2024, December). Retrieving Top-k Hyperedge Triplets: Models and Applications. In 2024 IEEE International Conference on Big Data (BigData) (pp. 630-639). IEEE.

    Parameters
    ----------
    H : hnx.Hypergraph
    k : int, optional, default : 1
        Top-k hyperedge triplets
    min_weight : int, optional, default : 0
        Minimum disjoint weight of returned hyperedge triplet

    Returns
    -------
    max_triplets : list of dictionaries
        dictionary has the following (key, value) pairs:
            ("weight", disjoint weight), 
            ("v1", first hyperedge), 
            ("v2": second hyperedge), 
            ("v3": third hyperedge)
        or empty list if no hyperedge triplet found
    """
    G, _, cdict = H.incidence_matrix(index=True)

    U, V, id_to_v = preprocessing(G, cdict, min_weight)

    maxWeight = min_weight - 1

    vRight = len(V)

    maxTriplets = []

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
                                maxTriplet = {
                                    "weight": disjointWeight,
                                    "v1": id_to_v[a],
                                    "v2": id_to_v[b],
                                    "v3": id_to_v[c]
                                }
                                add_to_sorted_list(maxTriplets, maxTriplet, "weight")
                                if len(maxTriplets) > k:
                                    maxTriplets.pop()
                                    maxWeight = maxTriplets[-1]["weight"]

        if db > maxWeight:
            for c, bc in overlap.items():
                dc = len(V[c])
                if dc > maxWeight:
                    intersections[c][b] = bc
    
    return maxTriplets

def max_common(H, k=1, min_weight=0):
    """
    Returns the hyperedge triplets with the highest common weights.

    Niu, J., Amburg, I. D., Aksoy, S. G., & Sarıyüce, A. E. (2024, December). Retrieving Top-k Hyperedge Triplets: Models and Applications. In 2024 IEEE International Conference on Big Data (BigData) (pp. 630-639). IEEE.

    Parameters
    ----------
    H : hnx.Hypergraph
    k : int, optional, default : 1
        Top-k hyperedge triplets
    min_weight : int, optional, default : 0
        Minimum common weight of returned hyperedge triplet

    Returns
    -------
    max_triplets : list of dictionaries
        dictionary has the following (key, value) pairs:
            ("weight", common weight), 
            ("v1", first hyperedge), 
            ("v2": second hyperedge), 
            ("v3": third hyperedge)
        or empty list if no hyperedge triplet found
    """
    G, _, cdict = H.incidence_matrix(index=True)

    U, V, id_to_v = preprocessing(G, cdict, min_weight)

    maxWeight = min_weight - 1

    vRight = len(V)

    maxTriplets = []

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
                                maxTriplet = {
                                    "weight": commonWeight,
                                    "v1": id_to_v[a],
                                    "v2": id_to_v[b],
                                    "v3": id_to_v[c]
                                }
                                add_to_sorted_list(maxTriplets, maxTriplet, "weight")
                                if len(maxTriplets) > k:
                                    maxTriplets.pop()
                                    maxWeight = maxTriplets[-1]["weight"]

        if db > maxWeight:
            for c, bc in overlap.items():
                dc = len(V[c])
                if dc > maxWeight:
                    intersections[c][b] = bc

    return maxTriplets

def local_independent(H, target_hyperedge, k=1, min_weight=0):
    """
    Returns the hyperedge triplets containing target_hyperedge with the highest independent weights.

    Niu, J., Amburg, I. D., Aksoy, S. G., & Sarıyüce, A. E. (2024, December). Retrieving Top-k Hyperedge Triplets: Models and Applications. In 2024 IEEE International Conference on Big Data (BigData) (pp. 630-639). IEEE.

    Parameters
    ----------
    H : hnx.Hypergraph
    target_hyperedge : string
        Target hyperedge for local traversal
    k : int, optional, default : 1
        Top-k hyperedge triplets
    min_weight : int, optional, default : 0
        Minimum independent weight of returned hyperedge triplet

    Returns
    -------
    max_triplets : list of dictionaries
        dictionary has the following (key, value) pairs:
            ("weight", independent weight), 
            ("v1", first hyperedge), 
            ("v2": second hyperedge), 
            ("v3": third hyperedge)
        or empty list if no hyperedge triplet found
    """
    G, _, cdict = H.incidence_matrix(index=True)

    U, V, id_to_v, target_hyperedge = preprocessing(G, cdict, min_weight, target_hyperedge=target_hyperedge)

    maxWeight = min_weight - 1
    a = target_hyperedge

    if target_hyperedge is None or len(V[a]) <= maxWeight:
        print("Target hyperedge does not have sufficient size")
        return []

    da = len(V[a])

    vRight = len(V)

    maxTriplets = []

    overlap_a = {}
    for u in V[a]:
        for b in U[u]:
            if b == a:
                continue
            db = len(V[b])
            if db <= maxWeight:
                break
            if b not in overlap_a:
                overlap_a[b] = [u]
            else:
                overlap_a[b].append(u)

    overlap_a = {}
    for u in V[a]:
        for b in U[u]:
            if b == a:
                continue
            db = len(V[b])
            if db <= maxWeight:
                break
            if b not in overlap_a:
                overlap_a[b] = [u]
            else:
                overlap_a[b].append(u)

    for b in range(vRight):

        if b == a:
            continue
        
        db = len(V[b])
        if db <= maxWeight:
            break

        ab = overlap_a[b] if b in overlap_a else []
        abl = len(ab)
        if (float) (min(da, db) - abl) / (float) (abl + 1) <= maxWeight:
            continue

        overlap = {}
        for u in V[b]:
            c_idx_start = bisect.bisect_right(U[u], b)
            Du = len(U[u])
            for c_idx in range(c_idx_start, Du):
                c = U[u][c_idx]
                if c == a:
                    continue
                dc = len(V[c])
                if dc <= maxWeight:
                    break
                if c not in overlap:
                    overlap[c] = [u]
                else:
                    overlap[c].append(u)

        for c in range(b + 1, vRight):

            if c == a:
                continue

            dc = len(V[c])
            if dc <= maxWeight:
                break

            bc = overlap[c] if c in overlap else []
            bcl = len(bc)
            
            if (float) (min(db, dc) - bcl) / (float) (bcl + 1) > maxWeight:

                ac = overlap_a[c] if c in overlap_a else []
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
                        maxTriplet = {
                            "weight": independentWeight,
                            "v1": id_to_v[a],
                            "v2": id_to_v[b],
                            "v3": id_to_v[c]
                        }
                        add_to_sorted_list(maxTriplets, maxTriplet, "weight")
                        if len(maxTriplets) > k:
                            maxTriplets.pop()
                            maxWeight = maxTriplets[-1]["weight"]

    return maxTriplets

def local_disjoint(H, target_hyperedge, k=1, min_weight=0):
    """
    Returns the hyperedge triplets containing target_hyperedge with the highest disjoint weights.

    Niu, J., Amburg, I. D., Aksoy, S. G., & Sarıyüce, A. E. (2024, December). Retrieving Top-k Hyperedge Triplets: Models and Applications. In 2024 IEEE International Conference on Big Data (BigData) (pp. 630-639). IEEE.

    Parameters
    ----------
    H : hnx.Hypergraph
    target_hyperedge : string
        Target hyperedge for local traversal
    k : int, optional, default : 1
        Top-k hyperedge triplets
    min_weight : int, optional, default : 0
        Minimum disjoint weight of returned hyperedge triplet

    Returns
    -------
    max_triplets : list of dictionaries
        dictionary has the following (key, value) pairs:
            ("weight", disjoint weight), 
            ("v1", first hyperedge), 
            ("v2": second hyperedge), 
            ("v3": third hyperedge)
        or empty list if no hyperedge triplet found
    """
    G, _, cdict = H.incidence_matrix(index=True)

    U, V, id_to_v, target_hyperedge = preprocessing(G, cdict, min_weight, target_hyperedge=target_hyperedge)

    maxWeight = min_weight - 1
    a = target_hyperedge

    if target_hyperedge is None or math.floor(len(V[a]) / 2) <= maxWeight:
        print("Target hyperedge does not have sufficient size")
        return []

    vRight = len(V)

    maxTriplets = []

    overlap_a = {}
    for u in V[a]:
        for b in U[u]:
            if b == a:
                continue
            db = len(V[b])
            if math.floor(db / 2) <= maxWeight:
                break
            if b not in overlap_a:
                overlap_a[b] = [u]
            else:
                overlap_a[b].append(u)

    for b in range(vRight):

        if b == a:
            continue

        db = len(V[b])
        if math.floor(db / 2) <= maxWeight:
            break
        
        ab = overlap_a[b] if b in overlap_a else []
        abl = len(ab)
        if abl <= maxWeight:
            continue

        overlap = {}
        for u in V[b]:
            c_idx_start = bisect.bisect_right(U[u], b)
            Du = len(U[u])
            for c_idx in range(c_idx_start, Du):
                c = U[u][c_idx]
                if c == a:
                    continue
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
                ac = overlap_a[c] if c in overlap_a else []
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
                        maxTriplet = {
                            "weight": disjointWeight,
                            "v1": id_to_v[a],
                            "v2": id_to_v[b],
                            "v3": id_to_v[c]
                        }
                        add_to_sorted_list(maxTriplets, maxTriplet, "weight")
                        if len(maxTriplets) > k:
                            maxTriplets.pop()
                            maxWeight = maxTriplets[-1]["weight"]

    return maxTriplets

def local_common(H, target_hyperedge, k=1, min_weight=0):
    """
    Returns the hyperedge triplets containing target_hyperedge with the highest common weights.

    Niu, J., Amburg, I. D., Aksoy, S. G., & Sarıyüce, A. E. (2024, December). Retrieving Top-k Hyperedge Triplets: Models and Applications. In 2024 IEEE International Conference on Big Data (BigData) (pp. 630-639). IEEE.

    Parameters
    ----------
    H : hnx.Hypergraph
    target_hyperedge : string
        Target hyperedge for local traversal
    k : int, optional, default : 1
        Top-k hyperedge triplets
    min_weight : int, optional, default : 0
        Minimum common weight of returned hyperedge triplet

    Returns
    -------
    max_triplets : list of dictionaries
        dictionary has the following (key, value) pairs:
            ("weight", common weight), 
            ("v1", first hyperedge), 
            ("v2": second hyperedge), 
            ("v3": third hyperedge)
        or empty list if no hyperedge triplet found
    """
    G, _, cdict = H.incidence_matrix(index=True)

    U, V, id_to_v, target_hyperedge = preprocessing(G, cdict, min_weight, target_hyperedge=target_hyperedge)

    maxWeight = min_weight - 1
    a = target_hyperedge

    if target_hyperedge is None or len(V[a]) <= maxWeight:
        print("Target hyperedge does not have sufficient size")
        return []

    vRight = len(V)

    maxTriplets = []

    overlap_a = {}
    for u in V[a]:
        for b in U[u]:
            if b == a:
                continue
            db = len(V[b])
            if db <= maxWeight:
                break
            if b not in overlap_a:
                overlap_a[b] = [u]
            else:
                overlap_a[b].append(u)

    for b in range(vRight):

        if b == a:
            continue

        db = len(V[b])
        if db <= maxWeight:
            break
        
        ab = overlap_a[b] if b in overlap_a else []
        abl = len(ab)
        if abl <= maxWeight:
            continue

        overlap = {}
        for u in V[b]:
            c_idx_start = bisect.bisect_right(U[u], b)
            Du = len(U[u])
            for c_idx in range(c_idx_start, Du):
                c = U[u][c_idx]
                if c == a:
                    continue
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
                ac = overlap_a[c] if c in overlap_a else []
                acl = len(ac)
                commonWeight = min(abl, acl, bcl)
                if commonWeight > maxWeight:
                    if abl < bcl:
                        commonWeight = intersection_size(ab, ac)
                    else:
                        commonWeight = intersection_size(bc, ac)
                    
                    if commonWeight > maxWeight:
                        maxTriplet = {
                            "weight": commonWeight,
                            "v1": id_to_v[a],
                            "v2": id_to_v[b],
                            "v3": id_to_v[c]
                        }
                        add_to_sorted_list(maxTriplets, maxTriplet, "weight")
                        if len(maxTriplets) > k:
                            maxTriplets.pop()
                            maxWeight = maxTriplets[-1]["weight"]
    
    return maxTriplets

def max_triplets(H, weight_type, k=1, min_weight=0):
    """
    Returns the hyperedge triplets with the highest weights.
    For more information, see:

    Niu, J., Amburg, I. D., Aksoy, S. G., & Sarıyüce, A. E. (2024, December). Retrieving Top-k Hyperedge Triplets: Models and Applications. In 2024 IEEE International Conference on Big Data (BigData) (pp. 630-639). IEEE.

    Parameters
    ----------
    H : hnx.Hypergraph
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
    
    if weight_type == "independent":
        return max_independent(H, k, min_weight)
    elif weight_type == "disjoint":
        return max_disjoint(H, k, min_weight)
    elif weight_type == "common":
        return max_common(H, k, min_weight)
    else:
        print("Invalid weight type, must be in [independent, disjoint, common]")
        return []

def local_triplets(H, weight_type, target_hyperedge, k=1, min_weight=0):
    """
    Returns the hyperedge triplets containing target_hyperedge with the highest weights.

    Niu, J., Amburg, I. D., Aksoy, S. G., & Sarıyüce, A. E. (2024, December). Retrieving Top-k Hyperedge Triplets: Models and Applications. In 2024 IEEE International Conference on Big Data (BigData) (pp. 630-639). IEEE.

    Parameters
    ----------
    H : hnx.Hypergraph
    weight_type : string
        Type in ["independent", "disjoint", common"]
    target_hyperedge : string
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
    
    if weight_type == "independent":
        return local_independent(H, target_hyperedge, k, min_weight)
    elif weight_type == "disjoint":
        return local_disjoint(H, target_hyperedge, k, min_weight)
    elif weight_type == "common":
        return local_common(H, target_hyperedge, k, min_weight)
    else:
        print("Invalid weight type, must be in [independent, disjoint, common]")
        return []
    