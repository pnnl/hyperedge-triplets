#include <vector>
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_utils.h"
#include "portable-memory-mapping-master/MemoryMapped.h"

typedef std::vector<std::vector<int>> graph;
typedef phmap::flat_hash_map<int, std::vector<int>> mapVec;
typedef std::vector<mapVec> vecMapVec;

// reads graph file and converts to adjacency list
void readGraph(const char* file_path, graph& U, graph& V) {

    MemoryMapped f(file_path);

    uint64_t size = f.size();

    int headerEnd = 0;
    int nEdge = 0, vLeft = 0, vRight = 0;
    bool nE = true, vL = false;
    char c = f[headerEnd];
    while (c != '\n') {
        if (isdigit(c)) {
            if (nE) {
                nEdge = nEdge * 10 + c - '0';
            }
            else if (vL) {
                vLeft = vLeft * 10 + c - '0';
            }
            else {
                vRight = vRight * 10 + c - '0';
            }
        }
        else {
            if (nE) {
                nE = false;
                vL = true;
            }
            else {
                vL = false;
            }
        }
        headerEnd += 1;
        c = f[headerEnd];
    }

    U.resize(vLeft);
    V.resize(vRight);

    int u = 0, v = 0;
    bool left = true;
    for (int i = headerEnd + 1; i < size; ++i) {
        c = f[i];

        if (isdigit(c)) {
            if (left) {
                u = u * 10 + c - '0';
            }
            else {
                v = v * 10 + c - '0';
            }
        }
        else {
            if (c == ' ') {
                left = false;
            }
            // if edge has been processed
            else if (!left) {
                left = true;
                if (std::find(U[u].begin(), U[u].end(), v) == U[u].end()) {
                    U[u].emplace_back(v);
                    V[v].emplace_back(u);
                }
                u = 0; v = 0;
            }
        }
    }
}

// (𝛼, 𝛽)-core decomposition and decreasing degree ordering
void preProcessing(graph& U, graph& V, const int Kv, std::vector<int>& id_to_v) {

    const int vLeft = U.size();
    const int vRight = V.size();

    int vLeft2 = vLeft;
    int vRight2 = vRight;
    
    // (𝛼 = 1, 𝛽 = Kv)-core decomposition of (U, V)
    std::vector<int> Du(vLeft);
    std::vector<int> Dv(vRight);
    std::vector<bool> Ru(vLeft);
    std::vector<bool> Rv(vRight);

    std::vector<int> idx(vRight);

    for (int u = 0; u < vLeft; ++u) {
        Du[u] = U[u].size();
    }

    for (int v = 0; v < vRight; ++v) {
        Dv[v] = V[v].size();
        idx[v] = v;
    }

    for (int v = 0; v < vRight; ++v) {
        if (!Rv[v] && Dv[v] < Kv) {
            Rv[v] = true;
            --vRight2;
            for (int u : V[v]) {
                if (!Ru[u]) {
                    --Du[u];
                    if (Du[u] == 0) {
                        Ru[u] = true;
                        --vLeft2;
                    }
                }
            }
        }
    }

    // calculate new ids for nodes in U by removing gaps from filtered nodes
    std::vector<int> u_to_id(vLeft);
    int max_id = 0;
    for (int u = 0; u < vLeft; ++u) {
        if (!Ru[u]) {
            u_to_id[u] = max_id;
            ++max_id;
        }
    }


    // calculate new ids for hyperedges in V by decreasing degree
    std::sort(idx.begin(), idx.end(), [&Dv](int v1, int v2) {return Dv[v1] > Dv[v2];});
    std::vector<int> v_to_id(vRight);
    id_to_v.resize(vRight2);
    max_id = 0;
    for (int v : idx) {
        if (Rv[v]) {
            break;
        }
        v_to_id[v] = max_id;
        id_to_v[max_id] = v;
        ++max_id;
    }

    // update U and V with new ids and sort neighbor lists in ascending order
    graph U2(vLeft2);
    graph V2(vRight2);
    for (int v = 0; v < vRight; ++v) {
        if (!Rv[v]) {
            const int v_id = v_to_id[v];
            for (int u : V[v]) {
                if (!Ru[u]) {
                    const int u_id = u_to_id[u];
                    U2[u_id].emplace_back(v_id);
                    V2[v_id].emplace_back(u_id);
                }
            }
            std::sort(V2[v_id].begin(), V2[v_id].end());
        }
    }
    for (int u = 0; u < vLeft2; ++u) {
        std::sort(U2[u].begin(), U2[u].end());
    }
    U = std::move(U2);
    V = std::move(V2);
}

// size of the intersection between two sorted vectors
int intersectionSize (const std::vector<int>& s1, const std::vector<int>& s2) {
  int size = 0;

  std::vector<int>::const_iterator first1 = s1.begin();
  std::vector<int>::const_iterator last1 = s1.end();
  std::vector<int>::const_iterator first2 = s2.begin();
  std::vector<int>::const_iterator last2 = s2.end();

  while (first1 != last1 && first2 != last2)
  {
    if (*first1 < *first2) ++first1;
    else if (*first2 < *first1) ++first2;
    else {
        ++size;
        ++first1;
        ++first2;
    }
  }
  return size;
}

// size of the intersection between two sorted vectors with early stopping based on an upper bound
int intersectionSizeBounded (const std::vector<int>& s1, const std::vector<int>& s2, const int upper_bound_excl) {
  int size = 0;

  std::vector<int>::const_iterator first1 = s1.begin();
  std::vector<int>::const_iterator last1 = s1.end();
  std::vector<int>::const_iterator first2 = s2.begin();
  std::vector<int>::const_iterator last2 = s2.end();

  while (first1 != last1 && first2 != last2)
  {
    if (*first1 < *first2) ++first1;
    else if (*first2 < *first1) ++first2;
    else {
        ++size;
        if (size >= upper_bound_excl)
            return -1;
        ++first1;
        ++first2;
    }
  }
  return size;
}

// hyperedge triplets
struct triplet {
    float weight;
    int v1;
    int v2;
    int v3;
};

// hyperedge triplet with maximum independent weight (limited early stopping)
extern "C" triplet* basicIndependent(const char* dataset_path, float maxWeight) {

    graph U, V;

    readGraph(dataset_path, U, V);

    std::vector<int> id_to_v;

    preProcessing(U, V, maxWeight + 1, id_to_v);

    const int vRight = V.size();

    triplet* maxTriplet = new triplet();
    maxTriplet->weight = maxWeight;

    vecMapVec intersections(vRight);

    for (int a = 0; a < vRight; ++a) {
        for (int u : V[a]) {
            const int b_idx_start = std::distance(U[u].begin(), std::upper_bound(U[u].begin(), U[u].end(), a));
            const int Du = U[u].size();
            for (int b_idx = b_idx_start; b_idx < Du; ++b_idx) {
                const int b = U[u][b_idx];
                intersections[a][b].emplace_back(u);
            }
        }
    }

    for (int a = 0; a < vRight; ++a) {
        const int da = V[a].size();
        if (da <= maxWeight) {
            break;
        }
        for (int b = a + 1; b < vRight; ++b) {
            const int db = V[b].size();
            const std::vector<int>& ab = (intersections[a].contains(b)) ? intersections[a][b] : std::vector<int>();
            const int abl = ab.size();
            for (int c = b + 1; c < vRight; ++c) {
                const int dc = V[c].size();
                const std::vector<int>& ac = (intersections[a].contains(c)) ? intersections[a][c] : std::vector<int>();
                const int acl = ac.size();
                const std::vector<int>& bc = (intersections[b].contains(c)) ? intersections[b][c] : std::vector<int>();
                const int bcl = bc.size();

                int x = std::min({da - abl - acl, db - abl - bcl, dc - acl - bcl});
                int abc = std::min({abl, bcl, acl});
                float independentWeight = (float) (x + abc) / (float) (abl + acl + bcl - 2 * abc + 1);

                if (std::min({abl, bcl, acl}) > 0) {
                    if (abl < bcl) {
                        abc = intersectionSize(ab, ac);
                        independentWeight = (float) (x + abc) / (float) (abl + acl + bcl - 2 * abc + 1);
                    }
                    else {
                        abc = intersectionSize(bc, ac);
                        independentWeight = (float) (x + abc) / (float) (abl + acl + bcl - 2 * abc + 1);
                    }
                }

                if (independentWeight > maxWeight) {
                    maxTriplet->weight = independentWeight;
                    maxTriplet->v1 = id_to_v[a];
                    maxTriplet->v2 = id_to_v[b];
                    maxTriplet->v3 = id_to_v[c];
                    maxWeight = independentWeight;
                }
            }
        }
    }
    return maxTriplet;
}

// hyperedge triplet with maximum disjoint weight (limited early stopping)
extern "C" triplet* basicDisjoint(const char* dataset_path, float maxWeight) {

    graph U, V;

    readGraph(dataset_path, U, V);

    std::vector<int> id_to_v;

    preProcessing(U, V, maxWeight + 1, id_to_v);

    const int vRight = V.size();

    triplet* maxTriplet = new triplet();
    maxTriplet->weight = maxWeight;

    vecMapVec intersections(vRight);
    vecMapVec intersectionsRev(vRight);

    for (int a = 0; a < vRight; ++a) {
        for (int u : V[a]) {
            const int b_idx_start = std::distance(U[u].begin(), std::upper_bound(U[u].begin(), U[u].end(), a));
            const int Du = U[u].size();
            for (int b_idx = b_idx_start; b_idx < Du; ++b_idx) {
                const int b = U[u][b_idx];
                intersections[a][b].emplace_back(u);
                intersectionsRev[b][a].emplace_back(u);
            }
        }
    }

    for (int b = 0; b < vRight; ++b) {
        const int db = V[b].size();
        if (std::floor(db / 2) <= maxWeight) {
            break;
        }
        for (mapVec::iterator b_inter_it = intersections[b].begin(); b_inter_it != intersections[b].end(); ++b_inter_it) {
            const int c = b_inter_it->first;
            const int dc = V[c].size();
            const std::vector<int>& bc = b_inter_it->second;
            const int bcl = bc.size();
            
            for (mapVec::iterator b_interRev_it = intersectionsRev[b].begin(); b_interRev_it != intersectionsRev[b].end(); ++b_interRev_it) {
                const int a = b_interRev_it->first;
                const std::vector<int>& ab = b_interRev_it->second;
                const int abl = ab.size();
                const int da = V[a].size();
                if (intersections[a].contains(c)) {
                    const std::vector<int>& ac = intersections[a][c];
                    const int acl = ac.size();
                    int abcl = 0;
                    if (abl < bcl) {
                        abcl += intersectionSize(ab, ac);
                    }
                    else {
                        abcl += intersectionSize(bc, ac);
                    }

                    float disjointWeight = (float) (std::min({abl, bcl, acl}) - abcl) / (float) (abcl + 1);

                    if (disjointWeight > maxWeight) {
                        maxTriplet->weight = disjointWeight;
                        maxTriplet->v1 = id_to_v[a];
                        maxTriplet->v2 = id_to_v[b];
                        maxTriplet->v3 = id_to_v[c];
                        maxWeight = disjointWeight;
                    }
                }
            }
        }
    }
    return maxTriplet;
}

// hyperedge triplet with maximum disjoint weight (limited early stopping)
extern "C" triplet* basicCommon(const char* dataset_path, int maxWeight) {

    graph U, V;

    readGraph(dataset_path, U, V);

    std::vector<int> id_to_v;

    preProcessing(U, V, maxWeight + 1, id_to_v);

    const int vRight = V.size();

    triplet* maxTriplet = new triplet();
    maxTriplet->weight = maxWeight;

    vecMapVec intersections(vRight);
    vecMapVec intersectionsRev(vRight);

    for (int a = 0; a < vRight; ++a) {
        for (int u : V[a]) {
            const int b_idx_start = std::distance(U[u].begin(), std::upper_bound(U[u].begin(), U[u].end(), a));
            const int Du = U[u].size();
            for (int b_idx = b_idx_start; b_idx < Du; ++b_idx) {
                const int b = U[u][b_idx];
                intersections[a][b].emplace_back(u);
                intersectionsRev[b][a].emplace_back(u);
            }
        }
    }

    for (int b = 0; b < vRight; ++b) {
        const int db = V[b].size();
        if (db <= maxWeight) {
            break;
        }
        for (mapVec::iterator b_inter_it = intersections[b].begin(); b_inter_it != intersections[b].end(); ++b_inter_it) {
            const int c = b_inter_it->first;
            const int dc = V[c].size();
            const std::vector<int>& bc = b_inter_it->second;
            const int bcl = bc.size();
            
            for (mapVec::iterator b_interRev_it = intersectionsRev[b].begin(); b_interRev_it != intersectionsRev[b].end(); ++b_interRev_it) {
                const int a = b_interRev_it->first;
                const std::vector<int>& ab = b_interRev_it->second;
                const int abl = ab.size();
                const int da = V[a].size();
                if (intersections[a].contains(c)) {
                    const std::vector<int>& ac = intersections[a][c];
                    const int acl = ac.size();

                    int commonWeight = 0;
                    if (abl < bcl) {
                        commonWeight = intersectionSize(ab, ac);
                    }
                    else {
                        commonWeight = intersectionSize(bc, ac);
                    }
                    if (commonWeight > maxWeight) {
                        maxTriplet->weight = commonWeight;
                        maxTriplet->v1 = id_to_v[a];
                        maxTriplet->v2 = id_to_v[b];
                        maxTriplet->v3 = id_to_v[c];
                        maxWeight = commonWeight;
                    }
                }
            }
        }
    }
    return maxTriplet;
}

// hyperedge triplet with maximum independent weight (with early stopping)
extern "C" triplet* maxIndependent(const char* dataset_path, float maxWeight) {

    graph U, V;

    readGraph(dataset_path, U, V);

    std::vector<int> id_to_v;

    preProcessing(U, V, maxWeight + 1, id_to_v);

    const int vRight = V.size();

    triplet* maxTriplet = new triplet();
    maxTriplet->weight = maxWeight;

    vecMapVec intersections(vRight);

    for (int b = 0; b < vRight; ++b) {
        const int db = V[b].size();
        if (db <= maxWeight) {
            break;
        }
        mapVec overlap;
        for (int u : V[b]) {
            const int c_idx_start = std::distance(U[u].begin(), std::upper_bound(U[u].begin(), U[u].end(), b));
            const int Du = U[u].size();
            for (int c_idx = c_idx_start; c_idx < Du; ++c_idx) {
                const int c = U[u][c_idx];
                const int dc = V[c].size();
                if (dc <= maxWeight) {
                    break;
                }
                overlap[c].emplace_back(u);
            }
        }
        for (int c = b + 1; c < vRight; ++c) {
            const int dc = V[c].size();
            if (dc <= maxWeight) {
                break;
            }

            const std::vector<int>& bc = (overlap.contains(c)) ? overlap[c] : std::vector<int>();
            const int bcl = bc.size();
            
            if ((float) (std::min(db, dc) - bcl) / (float) (bcl + 1) > maxWeight) {
                for (int a = 0; a < b; ++a) {
                    const int da = V[a].size();

                    const std::vector<int>& ab = (intersections[b].contains(a)) ? intersections[b][a] : std::vector<int>();
                    const int abl = ab.size();

                    if ((float) (std::min(da, db) - abl) / (float) (abl + 1) > maxWeight) {

                        const std::vector<int>& ac = (intersections[c].contains(a)) ? intersections[c][a] : std::vector<int>();
                        const int acl = ac.size();

                        const int x = std::min({da - abl - acl, db - abl - bcl, dc - acl - bcl});
                        int abcl = std::min({abl, bcl, acl});

                        float independentWeight = (float) (x + abcl) / (float)(abl + acl + bcl - 2 * abcl + 1);

                        if (independentWeight > maxWeight) {

                            if (abcl > 0) {
                                if (abl < bcl) {
                                    abcl = intersectionSize(ab, ac);
                                    independentWeight = (float) (x + abcl) / (float)(abl + acl + bcl - 2 * abcl + 1);
                                }
                                else {
                                    abcl = intersectionSize(bc, ac);
                                    independentWeight = (float) (x + abcl) / (float)(abl + acl + bcl - 2 * abcl + 1);
                                }
                            }

                            if (independentWeight > maxWeight) {
                                maxTriplet->weight = independentWeight;
                                maxTriplet->v1 = id_to_v[a];
                                maxTriplet->v2 = id_to_v[b];
                                maxTriplet->v3 = id_to_v[c];
                                maxWeight = independentWeight;
                            }
                        }

                    }
                }
            }

        }

        if (db > maxWeight) {
            for (mapVec::iterator overlap_it = overlap.begin(); overlap_it != overlap.end(); ++overlap_it) {
                const int c = overlap_it->first;
                const std::vector<int>& bc = overlap_it->second;
                const int dc = V[c].size();
                if (dc > maxWeight) {
                    intersections[c][b] = bc;
                }
            }
        }
    }

    return maxTriplet;
}

// hyperedge triplet with maximum common weight (with early stopping)
extern "C" triplet* maxDisjoint(const char* dataset_path, float maxWeight) {

    graph U, V;

    readGraph(dataset_path, U, V);

    std::vector<int> id_to_v;

    preProcessing(U, V, maxWeight + 1, id_to_v);

    const int vRight = V.size();

    triplet* maxTriplet = new triplet();
    maxTriplet->weight = maxWeight;

    vecMapVec intersections(vRight);

    for (int b = 0; b < vRight; ++b) {
        const int db = V[b].size();
        if (std::floor(db / 2) <= maxWeight) {
            break;
        }
        mapVec overlap;
        for (int u : V[b]) {
            const int c_idx_start = std::distance(U[u].begin(), std::upper_bound(U[u].begin(), U[u].end(), b));
            const int Du = U[u].size();
            for (int c_idx = c_idx_start; c_idx < Du; ++c_idx) {
                const int c = U[u][c_idx];
                const int dc = V[c].size();
                if (std::floor(dc / 2) <= maxWeight) {
                    break;
                }
                overlap[c].emplace_back(u);
            }
        }
        for (mapVec::iterator overlap_it = overlap.begin(); overlap_it != overlap.end(); ++overlap_it) {
            const std::vector<int>& bc = overlap_it->second;
            const int bcl = bc.size();
            if (bcl > maxWeight) {
                const int c = overlap_it->first;
                const int dc = V[c].size();
                for (mapVec::iterator b_inter_it = intersections[b].begin(); b_inter_it != intersections[b].end(); ++b_inter_it) {
                    const int a = b_inter_it->first;
                    const std::vector<int>& ab = b_inter_it->second;
                    const int abl = ab.size();
                    const int da = V[a].size();
                    if (std::floor(da / 2) > maxWeight && abl > maxWeight && intersections[c].contains(a)) {
                        const std::vector<int>& ac = intersections[c][a];
                        const int acl = ac.size();
                        
                        float disjointWeight = std::min({abl, bcl, acl});

                        if (disjointWeight > maxWeight) {
                            if (abl < bcl) {
                                const int abcl = intersectionSize(ab, ac);
                                disjointWeight = (float) (disjointWeight - abcl) / (float) (abcl + 1);
                            }
                            else {
                                const int abcl = intersectionSize(bc, ac);
                                disjointWeight = (float) (disjointWeight - abcl) / (float) (abcl + 1);
                            }
                            if (disjointWeight > maxWeight) {
                                maxTriplet->weight = disjointWeight;
                                maxTriplet->v1 = id_to_v[a];
                                maxTriplet->v2 = id_to_v[b];
                                maxTriplet->v3 = id_to_v[c];
                                maxWeight = disjointWeight;
                            }
                        }
                    }
                }
            }
        }

        if (std::floor(db / 2) > maxWeight) {
            for (mapVec::iterator overlap_it = overlap.begin(); overlap_it != overlap.end(); ++overlap_it) {
                const int c = overlap_it->first;
                const std::vector<int>& bc = overlap_it->second;
                const int dc = V[c].size();
                if (std::floor(dc / 2) > maxWeight && bc.size() > maxWeight) {
                    intersections[c][b] = bc;
                }
            }
        }

    }
    
    return maxTriplet;
}


// hyperedge triplet with maximum common weight (with early stopping)
extern "C" triplet* maxCommon(const char* dataset_path, int maxWeight) {

    graph U, V;

    readGraph(dataset_path, U, V);

    std::vector<int> id_to_v;

    preProcessing(U, V, maxWeight + 1, id_to_v);

    const int vRight = V.size();

    triplet* maxTriplet = new triplet();
    maxTriplet->weight = maxWeight;

    vecMapVec intersections(vRight);

    for (int b = 0; b < vRight; ++b) {
        const int db = V[b].size();
        if (db <= maxWeight) {
            break;
        }
        mapVec overlap;
        for (int u : V[b]) {
            const int c_idx_start = std::distance(U[u].begin(), std::upper_bound(U[u].begin(), U[u].end(), b));
            const int Du = U[u].size();
            for (int c_idx = c_idx_start; c_idx < Du; ++c_idx) {
                const int c = U[u][c_idx];
                const int dc = V[c].size();
                if (dc <= maxWeight) {
                    break;
                }
                overlap[c].emplace_back(u);
            }
        }
        for (mapVec::iterator overlap_it = overlap.begin(); overlap_it != overlap.end(); ++overlap_it) {
            const std::vector<int>& bc = overlap_it->second;
            const int bcl = bc.size();
            if (bcl > maxWeight) {
                const int c = overlap_it->first;
                const int dc = V[c].size();
                for (mapVec::iterator b_inter_it = intersections[b].begin(); b_inter_it != intersections[b].end(); ++b_inter_it) {
                    const int a = b_inter_it->first;
                    const std::vector<int>& ab = b_inter_it->second;
                    const int abl = ab.size();
                    const int da = V[a].size();
                    if (abl > maxWeight && intersections[c].contains(a)) {
                        const std::vector<int>& ac = intersections[c][a];
                        const int acl = ac.size();
                        int commonWeight = std::min({abl, acl, bcl});
                        if (commonWeight > maxWeight) {
                            if (abl < bcl) {
                                commonWeight = intersectionSize(ab, ac);
                            }
                            else {
                                commonWeight = intersectionSize(bc, ac);
                            }
                            if (commonWeight > maxWeight) {
                                maxTriplet->weight = commonWeight;
                                maxTriplet->v1 = id_to_v[a];
                                maxTriplet->v2 = id_to_v[b];
                                maxTriplet->v3 = id_to_v[c];
                                maxWeight = commonWeight;
                            }
                        }
                    }
                }
            }
        }

        if (db > maxWeight) {
            for (mapVec::iterator overlap_it = overlap.begin(); overlap_it != overlap.end(); ++overlap_it) {
                const int c = overlap_it->first;
                const std::vector<int>& bc = overlap_it->second;
                if (bc.size() > maxWeight) {
                    intersections[c][b] = bc;
                }
            }
        }

    }
    
    return maxTriplet;
}
