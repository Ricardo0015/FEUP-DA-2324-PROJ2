#include "Heap.h"

// Index calculations
#define parent(i) ((i) / 2)
#define leftChild(i) ((i) * 2)

/**
 * @brief Default constructor for the Heap class.
 * Initializes the heap with a sentinel value to facilitate parent/child calculations.
 */
Heap::Heap() {
    elems.push_back(0);
    // indices will be used starting in 1
    // to facilitate parent/child calculations
}

/**
 * @brief Constructor that builds a heap from a given vector of integers.
 * Uses the heapify-down process to enforce the heap property.
 * 
 * @param v The vector of integers to build the heap from.
 * 
 * @complexity O(n), where n is the number of elements in the vector.
 */
Heap::Heap(std::vector<int> v): Heap() {
    // Add the elements to the vector without respecting the heap property
    elems.insert(elems.end(), v.begin(), v.end());

    // Index of the last non-leaf node of vector v
    unsigned int startIdx = v.size() / 2;

    // Move the elements around so that the heap property is respected
    // There is no need to heapify the elements with higher indices since they are leaf nodes
    for(unsigned int i = startIdx; i > 0; i--) {
        heapifyDown(i);
    }
}

/**
 * @brief Checks if the heap is empty.
 * 
 * @return true if the heap is empty, false otherwise.
 * 
 * @complexity O(1)
 */
bool Heap::empty() {
    return elems.size() == 1;
}

/**
 * @brief Extracts and returns the minimum element from the heap.
 * Maintains the heap property after extraction by performing heapify-down.
 * 
 * @return The minimum element in the heap.
 * 
 * @complexity O(log n), where n is the number of elements in the heap.
 */
int Heap::extractMin() {
    auto x = elems[1];
    elems[1] = elems.back();
    elems.pop_back();
    if(elems.size() > 1) heapifyDown(1);
    return x;
}

/**
 * @brief Inserts a new element into the heap.
 * Maintains the heap property by performing heapify-up.
 * 
 * @param x The element to be inserted.
 * 
 * @complexity O(log n), where n is the number of elements in the heap.
 */
void Heap::insert(int x) {
    elems.push_back(x);
    heapifyUp(elems.size()-1);
}

/**
 * @brief Performs the heapify-up operation to maintain the heap property.
 * 
 * @param i The index of the element to be heapified up.
 * 
 * @complexity O(log n), where n is the number of elements in the heap.
 */
void Heap::heapifyUp(unsigned int i) {
    auto x = elems[i];
    while (i > 1 && x < elems[parent(i)]) {
        elems[i] = elems[parent(i)];
        i = parent(i);
    }
    elems[i] = x;
}

/**
 * @brief Performs the heapify-down operation to maintain the heap property.
 * 
 * @param i The index of the element to be heapified down.
 * 
 * @complexity O(log n), where n is the number of elements in the heap.
 */
void Heap::heapifyDown(unsigned int i) {
    auto x = elems[i];
    while (true) {
        unsigned k = leftChild(i);
        if (k >= elems.size())
            break;
        if (k+1 < elems.size() && elems[k+1] < elems[k])
            ++k; // right child of i
        if ( !(elems[k] < x) )
            break;
        elems[i] = elems[k];
        i = k;
    }
    elems[i] = x;
}