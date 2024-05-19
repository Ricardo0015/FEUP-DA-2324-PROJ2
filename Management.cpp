//
// Created by ricardo on 25-04-2024.
//

#include "Management.h"

using namespace std;

template<class T>
Management<T>::Management() = default;

template<class T>
const VTable<T> *Management<T>::getVertices() const {
    return vertices;
}

template<class T>
void Management<T>::setVertices(const VTable<T> *vertices) {
    this->vertices = vertices;
}

template<class T>
Vertex<T> Management<T>::vertexEntry(unsigned int id, double latitude, double longitude) {
    Vertex<T> v = Vertex<T>(id, points(latitude, longitude));
    vertices.insert(v);
    sumLatitude += latitude;
    sumLongitude += longitude;
    return v;
}

template<class T>
Vertex<T> Management<T>::findVertex(const unsigned int &id) {
    auto it = vertices.find(Vertex<T>(id));
    return it != vertices.end() ? *it : nullptr;
}

template<class T>
double Management<T>::getAverageLatitude() {
    return vertices.empty() ? 0 : sumLatitude / (double) vertices.size();
}

template<class T>
double Management<T>::getAverageLongitude() {
    return vertices.empty() ? 0 : sumLongitude / (double) vertices.size();
}

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

template<class T>
void Management<T>::dataClear() {
    vertices = {};
    sumLongitude = 0;
    sumLatitude = 0;
}