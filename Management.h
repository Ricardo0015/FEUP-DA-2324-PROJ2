//
// Created by ricardo on 25-04-2024.
//

#ifndef FEUP_DA_2324_PROJ2_MANAGEMENT_H
#define FEUP_DA_2324_PROJ2_MANAGEMENT_H
#include "Graph.h"
#include <unordered_set>
#include <cmath>
#include "points.h"

template <class T>
struct VHash {
    std::size_t operator()(Vertex<T> *vertex) const {
        return std::hash<unsigned int>()(vertex->getId());
    }
};

template <class T>
struct VEquals {
    bool operator()(Vertex<T> *vertex1, Vertex<T> *vertex2) const {
        return vertex1->getId() == vertex2->getId();
    }
};

template <class T>
using VTable =  std::unordered_set<Vertex<T>*, VHash<T>, VEquals<T>>;

/** class support **/

template <class T>
class Management {
private:
    VTable<T> vertices;
    double sumLatitude = 0;
    double sumLongitude = 0;
public:
    Management();

    const VTable<T> *getVertices() const;

    void setVertices(const VTable<T> *vertices);

    Vertex<T> findVertex(const unsigned int &id);

    double getAverageLatitude();

    double getAverageLongitude();

    Vertex<T> *furVertex();
    Vertex<T> vertexEntry(unsigned int id, double latitude, double longitude);

    void dataClear();
};

#endif //FEUP_DA_2324_PROJ2_MANAGEMENT_H
