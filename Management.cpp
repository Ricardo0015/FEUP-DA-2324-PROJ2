//
// Created by ricardo on 25-04-2024.
//

#include "Management.h"

using namespace std;

/**
 * @brief Default constructor for the Management class.
 */
template<class T>
Management<T>::Management() = default;

/**
 * @brief Gets the set of vertices.
 *
 * @return A constant pointer to the set of vertices.
 *
 * @complexity O(1)
 */
template<class T>
const VTable<T> *Management<T>::getVertices() const {
    return vertices;
}

/**
 * @brief Sets the set of vertices.
 *
 * @param vertices A pointer to the set of vertices to be assigned.
 *
 * @complexity O(1)
 */
template<class T>
void Management<T>::setVertices(const VTable<T> *vertices) {
    this->vertices = vertices;
}

/**
 * @brief Adds a vertex with the given id, latitude, and longitude to the management system.
 * Updates the sum of latitudes and longitudes.
 *
 * @param id The id of the vertex.
 * @param latitude The latitude of the vertex.
 * @param longitude The longitude of the vertex.
 * @return The newly added vertex.
 *
 * @complexity O(log n), where n is the number of vertices, due to the insertion into the set.
 */
template<class T>
Vertex<T> Management<T>::vertexEntry(unsigned int id, double latitude, double longitude) {
    Vertex<T> v = Vertex<T>(id, points(latitude, longitude));
    vertices.insert(v);
    sumLatitude += latitude;
    sumLongitude += longitude;
    return v;
}

/**
 * @brief Finds a vertex by its id.
 *
 * @param id The id of the vertex to be found.
 * @return The found vertex, or nullptr if not found.
 *
 * @complexity O(log n), where n is the number of vertices, due to the search in the set.
 */
template<class T>
Vertex<T> Management<T>::findVertex(const unsigned int &id) {
    auto it = vertices.find(Vertex<T>(id));
    return it != vertices.end() ? *it : nullptr;
}

/**
 * @brief Gets the average latitude of all vertices.
 *
 * @return The average latitude.
 *
 * @complexity O(1)
 */
template<class T>
double Management<T>::getAverageLatitude() {
    return vertices.empty() ? 0 : sumLatitude / (double) vertices.size();
}

/**
 * @brief Gets the average longitude of all vertices.
 *
 * @return The average longitude.
 *
 * @complexity O(1)
 */
template<class T>
double Management<T>::getAverageLongitude() {
    return vertices.empty() ? 0 : sumLongitude / (double) vertices.size();
}

/**
 * @brief Finds the vertex that is furthest from the average location.
 *
 * @return A pointer to the furthest vertex.
 *
 * @complexity O(n log n), where n is the number of vertices, due to the sort operation.
 */
template<class T>
Vertex<T> *Management<T>::furVertex() {
    std::vector<Vertex<T>> o(vertices.begin(), vertices.end());
    points averageLocation(getAverageLatitude(), getAverageLongitude());
    std::sort(o.begin(), o.end(),
              [averageLocation](Vertex<T> *n1, Vertex<T> *n2) {
                  points c1(n1->getCoordinates().getLatitude(), n1->getCoordinates().getLongitude());
                  points c2(n2->getCoordinates().getLatitude(), n2->getCoordinates().getLongitude());
                  return averageLocation.haversine(c1) > averageLocation.haversine(c2);
              });
    return *o[0];
}

/**
 * @brief Clears all the data in the management system.
 * Resets the sum of latitudes and longitudes.
 *
 * @complexity O(1)
 */
template<class T>
void Management<T>::dataClear() {
    vertices = {};
    sumLongitude = 0;
    sumLatitude = 0;
}