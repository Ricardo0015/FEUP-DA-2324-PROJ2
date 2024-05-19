#include "UFDS.h"

/**
 * @brief Constructor for the UFDS (Union-Find Disjoint Sets) class.
 * Initializes the data structures for path compression and union by rank.
 *
 * @param N The number of elements in the disjoint set.
 *
 * @complexity O(N)
 */
UFDS::UFDS(unsigned int N) {
    path.resize(N);
    rank.resize(N);
    for (unsigned long i = 0; i < N; i++) {
        path[i] = i;
        rank[i] = 0;
    }
}

/**
 * @brief Finds the representative (root) of the set containing element i.
 * Uses path compression to flatten the structure for future queries.
 *
 * @param i The element to find the set for.
 * @return The representative of the set containing i.
 *
 * @complexity O(log N) amortized, due to path compression.
 */
unsigned long UFDS::findSet(unsigned int i) {
    if (path[i] != i) path[i] = findSet(path[i]);
    return path[i];
}

/**
 * @brief Checks if elements i and j belong to the same set.
 *
 * @param i The first element.
 * @param j The second element.
 * @return true if both elements are in the same set, false otherwise.
 *
 * @complexity O(log N) amortized, due to the findSet operations.
 */
bool UFDS::isSameSet(unsigned int i, unsigned int j) {
    return findSet(i) == findSet(j);
}

/**
 * @brief Unites the sets containing elements i and j.
 * Uses union by rank to keep the tree flat.
 *
 * @param i The first element.
 * @param j The second element.
 *
 * @complexity O(log N) amortized, due to path compression and rank updates.
 */
void UFDS::linkSets(unsigned int i, unsigned int j) {
    if (!isSameSet(i, j)) {
        unsigned long x = findSet(i), y = findSet(j);
        if (rank[x] > rank[y]) path[y] = x; // x becomes the root due to having a larger rank
        else {
            path[x] = y; // y becomes the root due to having a larger rank, or ...
            if (rank[x] == rank[y]) rank[y]++; // ... due to both nodes having the same rank (in order to break the tie)
        }
    }
}
