// Original code by Gonçalo Leão
// Updated by DA 2023/2024 Team

#ifndef DA_TP_CLASSES_GRAPH
#define DA_TP_CLASSES_GRAPH


#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <memory>
#include <list>
#include "MutablePriorityQueue.h"
#include "points.h"
#include "UFDS.h"
template <class T>
class Edge;

#define INF std::numeric_limits<double>::max()

/************************* Vertex  **************************/

template <class T>
class Vertex {
public:
    Vertex(T in, points c);
    bool operator<(Vertex<T> & vertex) const; // // required by MutablePriorityQueue

    T getInfo() const;
    std::vector<Edge<T> *> getAdj() const;
    bool isVisited() const;
    bool isProcessing() const;
    unsigned int getIndegree() const;
    double getDist() const;
    Edge<T> *getPath() const;
    double haversineDistance(Vertex<T> *v) const;
    const points &getCoordinates() const;
    void setCoordinates(const points &newCoordinates);
    std::vector<Edge<T> *> getIncoming() const;

    void setInfo(T info);
    void setVisited(bool visited);
    void setProcesssing(bool processing);
    void setIndegree(unsigned int indegree);
    void setDist(double dist);
    void setPath(Edge<T> *path);
    Edge<T> * addEdge(Vertex<T> *dest, double w);
    bool removeEdge(T in);
    void removeOutgoingEdges();

    friend class MutablePriorityQueue<Vertex>;
protected:
    T info;                // info node
    std::vector<Edge<T> *> adj;  // outgoing edges
    points pointex;
    // auxiliary fields
    bool visited = false; // used by DFS, BFS, Prim ...
    bool processing = false; // used by isDAG (in addition to the visited attribute)
    unsigned int indegree; // used by topsort
    double dist = 0;
    Edge<T> *path = nullptr;

    std::vector<Edge<T> *> incoming; // incoming edges

    int queueIndex = 0; 		// required by MutablePriorityQueue and UFDS

    void deleteEdge(Edge<T> *edge);
};

/********************** Edge  ****************************/

template <class T>
class Edge {
public:
    Edge(Vertex<T> *orig, Vertex<T> *dest, double w);

    Vertex<T> * getDest() const;
    double getWeight() const;
    bool isSelected() const;
    Vertex<T> * getOrig() const;
    Edge<T> *getReverse() const;
    double getFlow() const;

    void setSelected(bool selected);
    void setReverse(Edge<T> *reverse);
    void setFlow(double flow);
protected:
    Vertex<T> * dest; // destination vertex
    double weight; // edge weight, can also be used for capacity

    // auxiliary fields
    bool selected = false;

    // used for bidirectional edges
    Vertex<T> *orig;
    Edge<T> *reverse = nullptr;

    double flow; // for flow-related problems
};

/********************** Graph  ****************************/

template <class T>
class Graph {
public:
    ~Graph();
    /*
    * Auxiliary function to find a vertex with a given the content.
    */
    Vertex<T> *findVertex(const T &in) const;
    /*
     *  Adds a vertex with a given content or info (in) to a graph (this).
     *  Returns true if successful, and false if a vertex with that content already exists.
     */
    Vertex<T> addVertex(const T &in, points c);
    bool removeVertex(const T &in);

    /*
     * Adds an edge to a graph (this), given the contents of the source and
     * destination vertices and the edge weight (w).
     * Returns true if successful, and false if the source or destination vertex does not exist.
     */
    bool addEdge(const T &sourc, const T &dest, double w);
    bool removeEdge(const T &source, const T &dest);
    bool addBidirectionalEdge(const T &sourc, const T &dest, double w);

    int getNumVertex() const;
    std::vector<Vertex<T> *> getVertexSet() const;

    std:: vector<T> dfs() const;
    std:: vector<T> dfs(const T & source) const;
    void dfsVisit(Vertex<T> *v,  std::vector<T> & res) const;
    std::vector<T> bfs(const T & source) const;

    bool isDAG() const;
    bool dfsIsDAG(Vertex<T> *v) const;
    std::vector<T> topsort() const;
    void prim();
    int addTour(Vertex<T>* s);
    int preorderMSTTraversal(Vertex<T> *source);
    void triangularTSPTour();
    void printTour();
    void printTour(const std::vector<unsigned int> &tour);
    void printTour(unsigned int *tour);
    bool solution(unsigned int s, const std::vector<unsigned int> &sol, unsigned int n);
    std::pair<double, std::vector<unsigned int>> tspBT();
    void tspRecursion(std::vector<unsigned int> &cSol, double cSolDist, unsigned int cNodeIndex, double &bSolDist, std::vector<unsigned int> &bSol, unsigned int n);
    std::pair<double, std::vector<unsigned int>> nearestHeuristic(unsigned int &start);
    std::pair<unsigned int, unsigned int> getNextHeuristicEdge(std::vector<unsigned int> tour, UFDS tSets);
    std::pair<std::vector<unsigned int>, double> getInsertionEdges(std::vector<unsigned int> tour, const unsigned int newVertexId) const;
    void graphClear();
    double getTourDistance() const;
protected:
    std::vector<Vertex<T> *> vertexSet;    // vertex set
    struct tour_t {
        double distance;
        std::list<Vertex<T>> *course;
    };

    tour_t tour = {0, {}};
    unsigned int totalEdges = 0;
    double ** distMatrix = nullptr;   // dist matrix for Floyd-Warshall
    int **pathMatrix = nullptr;   // path matrix for Floyd-Warshall
    std::vector<std::vector<bool>> selectedEdges;
    std::vector<std::vector<double>> distanceMatrix;
    /*
     * Finds the index of the vertex with a given content.
     */
    int findVertexIdx(const T &in) const;
};

void deleteMatrix(int **m, int n);
void deleteMatrix(double **m, int n);


/************************* Vertex  **************************/

template <class T>
Vertex<T>::Vertex(T in, points c): info(in), pointex(c) {}
/*
 * Auxiliary function to add an outgoing edge to a vertex (this),
 * with a given destination vertex (d) and edge weight (w).
 */
template <class T>
Edge<T> * Vertex<T>::addEdge(Vertex<T> *d, double w) {
    auto newEdge = new Edge<T>(this, d, w);
    adj.push_back(newEdge);
    d->incoming.push_back(newEdge);
    return newEdge;
}

/*
 * Auxiliary function to remove an outgoing edge (with a given destination (d))
 * from a vertex (this).
 * Returns true if successful, and false if such edge does not exist.
 */
template <class T>
bool Vertex<T>::removeEdge(T in) {
    bool removedEdge = false;
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge<T> *edge = *it;
        Vertex<T> *dest = edge->getDest();
        if (dest->getInfo() == in) {
            it = adj.erase(it);
            deleteEdge(edge);
            removedEdge = true; // allows for multiple edges to connect the same pair of vertices (multigraph)
        }
        else {
            it++;
        }
    }
    return removedEdge;
}

/*
 * Auxiliary function to remove an outgoing edge of a vertex.
 */
template <class T>
void Vertex<T>::removeOutgoingEdges() {
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge<T> *edge = *it;
        it = adj.erase(it);
        deleteEdge(edge);
    }
}

template <class T>
bool Vertex<T>::operator<(Vertex<T> & vertex) const {
    return this->dist < vertex.dist;
}

template <class T>
T Vertex<T>::getInfo() const {
    return this->info;
}

template <class T>
std::vector<Edge<T>*> Vertex<T>::getAdj() const {
    return this->adj;
}

template <class T>
bool Vertex<T>::isVisited() const {
    return this->visited;
}

template <class T>
bool Vertex<T>::isProcessing() const {
    return this->processing;
}

template <class T>
unsigned int Vertex<T>::getIndegree() const {
    return this->indegree;
}

template <class T>
double Vertex<T>::getDist() const {
    return this->dist;
}

template <class T>
Edge<T> *Vertex<T>::getPath() const {
    return this->path;
}

template <class T>
std::vector<Edge<T> *> Vertex<T>::getIncoming() const {
    return this->incoming;
}

template <class T>
void Vertex<T>::setInfo(T in) {
    this->info = in;
}

template <class T>
void Vertex<T>::setVisited(bool visited) {
    this->visited = visited;
}

template <class T>
void Vertex<T>::setProcesssing(bool processing) {
    this->processing = processing;
}

template <class T>
void Vertex<T>::setIndegree(unsigned int indegree) {
    this->indegree = indegree;
}

template <class T>
void Vertex<T>::setDist(double dist) {
    this->dist = dist;
}

template <class T>
void Vertex<T>::setPath(Edge<T> *path) {
    this->path = path;
}

template<class T>
double Vertex<T>::haversineDistance(Vertex<T> *v) const{
    return pointex.haversine(v->getCoordinates());
}

template<class T>
const points &Vertex<T>::getCoordinates() const {
    return pointex;
}

template <class T>
void Vertex<T>::setCoordinates(const points &newCoordinates) {
    Vertex::pointex = newCoordinates;
}

template <class T>
void Vertex<T>::deleteEdge(Edge<T> *edge) {
    Vertex<T> *dest = edge->getDest();
    // Remove the corresponding edge from the incoming list
    auto it = dest->incoming.begin();
    while (it != dest->incoming.end()) {
        if ((*it)->getOrig()->getInfo() == info) {
            it = dest->incoming.erase(it);
        }
        else {
            it++;
        }
    }
    delete edge;
}

/********************** Edge  ****************************/

template <class T>
Edge<T>::Edge(Vertex<T> *orig, Vertex<T> *dest, double w): orig(orig), dest(dest), weight(w) {}

template <class T>
Vertex<T> * Edge<T>::getDest() const {
    return this->dest;
}

template <class T>
double Edge<T>::getWeight() const {
    return this->weight;
}

template <class T>
Vertex<T> * Edge<T>::getOrig() const {
    return this->orig;
}

template <class T>
Edge<T> *Edge<T>::getReverse() const {
    return this->reverse;
}

template <class T>
bool Edge<T>::isSelected() const {
    return this->selected;
}

template <class T>
double Edge<T>::getFlow() const {
    return flow;
}

template <class T>
void Edge<T>::setSelected(bool selected) {
    this->selected = selected;
}

template <class T>
void Edge<T>::setReverse(Edge<T> *reverse) {
    this->reverse = reverse;
}

template <class T>
void Edge<T>::setFlow(double flow) {
    this->flow = flow;
}

/********************** Graph  ****************************/

template <class T>
int Graph<T>::getNumVertex() const {
    return vertexSet.size();
}

template <class T>
std::vector<Vertex<T> *> Graph<T>::getVertexSet() const {
    return vertexSet;
}

/*
 * Auxiliary function to find a vertex with a given content.
 */
template <class T>
Vertex<T> * Graph<T>::findVertex(const T &in) const {
    for (auto v : vertexSet)
        if (v->getInfo() == in)
            return v;
    return nullptr;
}

/*
 * Finds the index of the vertex with a given content.
 */
template <class T>
int Graph<T>::findVertexIdx(const T &in) const {
    for (unsigned i = 0; i < vertexSet.size(); i++)
        if (vertexSet[i]->getInfo() == in)
            return i;
    return -1;
}
/*
 *  Adds a vertex with a given content or info (in) to a graph (this).
 *  Returns true if successful, and false if a vertex with that content already exists.
 */
template <class T>
Vertex<T> Graph<T>::addVertex(const T &in, points c) {
    Vertex<T> *newVertex = nullptr;
    if (vertexSet.size() <= in) { vertexSet.resize(in + 1); }
    newVertex = std::make_shared<Vertex>(in, c);
    vertexSet[in] = newVertex;

    return newVertex;
}

/*
 *  Removes a vertex with a given content (in) from a graph (this), and
 *  all outgoing and incoming edges.
 *  Returns true if successful, and false if such vertex does not exist.
 */
template <class T>
bool Graph<T>::removeVertex(const T &in) {
    for (auto it = vertexSet.begin(); it != vertexSet.end(); it++) {
        if ((*it)->getInfo() == in) {
            auto v = *it;
            v->removeOutgoingEdges();
            for (auto u : vertexSet) {
                u->removeEdge(v->getInfo());
            }
            vertexSet.erase(it);
            delete v;
            return true;
        }
    }
    return false;
}

/*
 * Adds an edge to a graph (this), given the contents of the source and
 * destination vertices and the edge weight (w).
 * Returns true if successful, and false if the source or destination vertex does not exist.
 */
template <class T>
bool Graph<T>::addEdge(const T &sourc, const T &dest, double w) {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr)
        return false;
    v1->addEdge(v2, w);
    return true;
}

/*
 * Removes an edge from a graph (this).
 * The edge is identified by the source (sourc) and destination (dest) contents.
 * Returns true if successful, and false if such edge does not exist.
 */
template <class T>
bool Graph<T>::removeEdge(const T &sourc, const T &dest) {
    Vertex<T> * srcVertex = findVertex(sourc);
    if (srcVertex == nullptr) {
        return false;
    }
    return srcVertex->removeEdge(dest);
}

template <class T>
bool Graph<T>::addBidirectionalEdge(const T &sourc, const T &dest, double w) {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr)
        return false;
    auto e1 = v1->addEdge(v2, w);
    auto e2 = v2->addEdge(v1, w);
    e1->setReverse(e2);
    e2->setReverse(e1);
    return true;
}

/****************** DFS ********************/

/*
 * Performs a depth-first search (dfs) traversal in a graph (this).
 * Returns a vector with the contents of the vertices by dfs order.
 */
template <class T>
std::vector<T> Graph<T>::dfs() const {
    std::vector<T> res;
    for (auto v : vertexSet)
        v->setVisited(false);
    for (auto v : vertexSet)
        if (!v->isVisited())
            dfsVisit(v, res);
    return res;
}

/*
 * Performs a depth-first search (dfs) in a graph (this) from the source node.
 * Returns a vector with the contents of the vertices by dfs order.
 */
template <class T>
std::vector<T> Graph<T>::dfs(const T & source) const {
    std::vector<int> res;
    // Get the source vertex
    auto s = findVertex(source);
    if (s == nullptr) {
        return res;
    }
    // Set that no vertex has been visited yet
    for (auto v : vertexSet) {
        v->setVisited(false);
    }
    // Perform the actual DFS using recursion
    dfsVisit(s, res);

    return res;
}

/*
 * Auxiliary function that visits a vertex (v) and its adjacent, recursively.
 * Updates a parameter with the list of visited node contents.
 */
template <class T>
void Graph<T>::dfsVisit(Vertex<T> *v, std::vector<T> & res) const {
    v->setVisited(true);
    res.push_back(v->getInfo());
    for (auto & e : v->getAdj()) {
        auto w = e->getDest();
        if (!w->isVisited()) {
            dfsVisit(w, res);
        }
    }
}

/****************** BFS ********************/
/*
 * Performs a breadth-first search (bfs) in a graph (this), starting
 * from the vertex with the given source contents (source).
 * Returns a vector with the contents of the vertices by bfs order.
 */
template <class T>
std::vector<T> Graph<T>::bfs(const T & source) const {
    std::vector<int> res;
    // Get the source vertex
    auto s = findVertex(source);
    if (s == nullptr) {
        return res;
    }

    // Set that no vertex has been visited yet
    for (auto v : vertexSet) {
        v->setVisited(false);
    }

    // Perform the actual BFS using a queue
    std::queue<Vertex<T> *> q;
    q.push(s);
    s->setVisited(true);
    while (!q.empty()) {
        auto v = q.front();
        q.pop();
        res.push_back(v->getInfo());
        for (auto & e : v->getAdj()) {
            auto w = e->getDest();
            if ( ! w->isVisited()) {
                q.push(w);
                w->setVisited(true);
            }
        }
    }
    return res;
}

/****************** isDAG  ********************/
/*
 * Performs a depth-first search in a graph (this), to determine if the graph
 * is acyclic (acyclic directed graph or DAG).
 * During the search, a cycle is found if an edge connects to a vertex
 * that is being processed in the stack of recursive calls (see theoretical classes).
 * Returns true if the graph is acyclic, and false otherwise.
 */

template <class T>
bool Graph<T>::isDAG() const {
    for (auto v : vertexSet) {
        v->setVisited(false);
        v->setProcesssing(false);
    }
    for (auto v : vertexSet) {
        if (! v->isVisited()) {
            if ( ! dfsIsDAG(v) ) return false;
        }
    }
    return true;
}

/**
 * Auxiliary function that visits a vertex (v) and its adjacent, recursively.
 * Returns false (not acyclic) if an edge to a vertex in the stack is found.
 */
template <class T>
bool Graph<T>::dfsIsDAG(Vertex<T> *v) const {
    v->setVisited(true);
    v->setProcesssing(true);
    for (auto e : v->getAdj()) {
        auto w = e->getDest();
        if (w->isProcessing()) return false;
        if (! w->isVisited()) {
            if (! dfsIsDAG(w)) return false;
        }
    }
    v->setProcesssing(false);
    return true;
}

/****************** toposort ********************/
//=============================================================================
// Exercise 1: Topological Sorting
//=============================================================================
// TODO
/*
 * Performs a topological sorting of the vertices of a graph (this).
 * Returns a vector with the contents of the vertices by topological order.
 * If the graph has cycles, returns an empty vector.
 * Follows the algorithm described in theoretical classes.
 */

template<class T>
std::vector<T> Graph<T>::topsort() const {
    std::vector<int> res;

    for (auto v : vertexSet) {
        v->setIndegree(0);
    }
    for (auto v : vertexSet) {
        for (auto e : v->getAdj()) {
            unsigned int indegree = e->getDest()->getIndegree();
            e->getDest()->setIndegree(indegree + 1);
        }
    }

    std::queue<Vertex<T> *> q;
    for (auto v : vertexSet) {
        if (v->getIndegree() == 0) {
            q.push(v);
        }
    }

    while( !q.empty() ) {
        Vertex<T> * v = q.front();
        q.pop();
        res.push_back(v->getInfo());
        for(auto e : v->getAdj()) {
            auto w = e->getDest();
            w->setIndegree(w->getIndegree() - 1);
            if(w->getIndegree() == 0) {
                q.push(w);
            }
        }
    }

    if ( res.size() != vertexSet.size() ) {
        //std::cout << "Impossible topological ordering!" << std::endl;
        res.clear();
        return res;
    }

    return res;
}

template <class T>
void Graph<T>::prim() {
    MutablePriorityQueue<Vertex<T>*> q;
    std::vector<std::vector<bool>> newMatrix(vertexSet.size(), std::vector<bool>(vertexSet.size(), false));
    this->selectedEdges = newMatrix;

    Vertex<T>* start = findVertex(0);

    for (Vertex<T>* v: vertexSet) {
        v->setDist(INF);
        v->setVisited(false);
        v->setPath(nullptr);
    }
    start->setDist(0);
    q.insert(start);

    while (!q.empty()) {
        Vertex<T>* v = q.extractMin();
        if (v->getPath() != nullptr){
            selectedEdges[findVertexIdx(v->getInfo())][findVertexIdx(v->getPath()->getDest()->getInfo())] = true;
            selectedEdges[findVertexIdx(v->getPath()->getDest()->getInfo())][findVertexIdx(v->getInfo())] = true;
        }
        v->setVisited(true);
        for (size_t i = 0; i < vertexSet.size(); i++){
            if(distMatrix[v->getInfo()][i] == INF || i == v->getInfo()) continue;
            Vertex<T>* w = findVertex(i);
            if(!w->isVisited()){
                double distAntiga = w->getDist();
                if (distMatrix[v->getInfo()][i] < distAntiga){
                    w->setPath(v);
                    w->setDist(distMatrix[v->getInfo()][i]);
                    distAntiga = INF ? q.insert(w) : q.decreaseKey(w);
                }
            }
        }
    }
}

template <class T>
int Graph<T>::addTour(Vertex<T>* s) {
    if (!tour.course.empty()) {
        if ((*(tour.course)->rbegin())->getInfo() == s->getInfo()) {
            return -2;
        }
        double edgeWeight = findEdge((*(tour.course).rbegin()), s);
        if (edgeWeight == -1) {
            return (int) edgeWeight;
        }
        tour.distance += edgeWeight;
    }
    tour.push_back(s);
    return 0;
}

template <class T>
int Graph<T>::preorderMSTTraversal(Vertex<T> *source) {
    source->setVisited(true);
    int val = addTour(source);
    if (val != 0) return val;

    for (size_t i = 0; i < vertexSet.size(); i++) {
        Vertex<T> *dest = findVertex(i);

        if (!(dest->isVisited()) && selectedEdges[source->getId()][i]) {
            preorderMSTTraversal(dest);
        }
    }

    if (source->getId() == 0) return addToTour(source);
    return 0;
}

template <class T>
void Graph<T>::triangularTSPTour() {
    tour = {0, {}};
    prim();

    for (Vertex<T> *v: vertexSet) v->setVisited(false);

    int exec_val = preorderMSTTraversal(findVertex(0));
    switch (exec_val) {
        case -1:
            printf("Couldn't calculate approximation of TSP for this graph!\n");
            break;
        case -2:
            printf("Course would contain self-loop!\n");
            break;
        default:
            break;
    }
}

 template<class T>
void Graph<T>::printTour() {
    printf("Path choosed: ");
    for (Vertex<T> *v: tour.course) {
        printf(" %d", v->getId());
    }
    printf("\n");
}


template<class T>
void Graph<T>::printTour(const std::vector<unsigned int> &tour) {
    printf("Path choosed: ");
    for (const unsigned int &v: tour) {
        printf(" %d", v);
    }
    printf("\n");
}

template<class T>
void Graph<T>::printTour(unsigned int *tour) {
    printf("Path choosed: ");
    for (int i = 0; i < vertexSet.size(); i++) {
        printf(" %d", tour[i]);
    }
    printf("\n");
}

template <class T>
bool Graph<T>::solution(unsigned int s, const std::vector<unsigned int> &sol, unsigned int n) {
    for (int i = 0; i < n; i++) {
        if (sol[i] == s) {
            return true;
        }
    }
    return false;
}

template <class T>
std::pair<double, std::vector<unsigned int>> Graph<T>::tspBT() {
    unsigned int n = this->vertexSet.size();
    std::vector<unsigned int> currentSolution(n + 1);
    std::vector<unsigned int> path(n + 1);
    currentSolution[0] = 0;
    double bestSolutionDist = INF;
    tspRecursion(currentSolution, 0, 1, bestSolutionDist, path, n);
    path[n] = 0;
    return {bestSolutionDist, path};
}

template <class T>
void Graph<T>::tspRecursion(std::vector<unsigned int> &cSol, double cSolDist,
                    unsigned int cNodeIndex,
                    double &bSolDist, std::vector<unsigned int> &bSol, unsigned int n) {
    if (cNodeIndex == n) {
        //Could need to verify here if last node connects to first
        if (this->distanceMatrix[cSol[cNodeIndex - 1]][0] != INF) {
            //Add dist from last node back to zero and check if it's an improvement
            if (cSolDist + this->distanceMatrix[cSol[cNodeIndex - 1]][0] < bSolDist) {
                bSolDist = bSolDist + this->distanceMatrix[cSol[cNodeIndex - 1]][0];
                for (int i = 0; i < n; i++) {
                    bSol[i] = cSol[i];
                }
            }
        }
        return;

    }
    for (int i = 1; i < n; i++) {
        if (i == 12) {
            int a = 1;
        }
        if (this->distanceMatrix[cSol[cNodeIndex - 1]].size() > i &&
            this->distanceMatrix[cSol[cNodeIndex - 1]][i] + cSolDist < bSolDist) {
            if (!solution(i, cSol, cNodeIndex)) {
                cSol[cNodeIndex] = i;
                tspRecursion(cSol,
                             this->distanceMatrix[cSol[cNodeIndex - 1]][i] + cSolDist,
                             cNodeIndex + 1, bSolDist, bSol, n);
            }
        }
    }
}

template <class T>
std::pair<double, std::vector<unsigned int>> Graph<T>::nearestHeuristic(unsigned int &start) {
    double distance = 0;
    std::vector<unsigned int> tour;
    UFDS tSets(vertexSet.size());

    //Get shortest adjacent edge
    std::vector<double> adjacent = distanceMatrix[start];
    auto it = std::min_element(adjacent.begin(), adjacent.end());
    unsigned int minEdgeIndex = std::distance(adjacent.begin(), it);

    //Initialize the partial tour with the chosen vertex and its closest neighbour
    tour.push_back(start);
    tour.push_back(minEdgeIndex);
    tSets.linkSets(start, minEdgeIndex);
    distance += distanceMatrix[start][minEdgeIndex];

    //Two cities are already in the tour, repeat for the leftover cities
    for (int i = 2; i < vertexSet.size(); i++) {
        std::pair<unsigned int, unsigned int> nextEdge = getNextHeuristicEdge(tour, tSets);
        unsigned int newVertexId = nextEdge.second;

        auto insertionEdges = getInsertionEdges(tour, newVertexId);

        unsigned int closingVertex = insertionEdges.first.back();
        auto closingVertexIt = std::find(tour.begin(), tour.end(), closingVertex);
        tour.insert(closingVertexIt, newVertexId); //Insert new vertex in between the two old vertices

        tSets.linkSets(start, newVertexId);

        //Remove the length of the edge that was replaced, and add the length of the two new edges
        distance =
                distance - distanceMatrix[insertionEdges.first[0]][insertionEdges.first.back()] + insertionEdges.second;
    }
    //The two untied edges will always be the starting two vertices
    tour.push_back(start);
    return {distance + distanceMatrix[minEdgeIndex][start], tour};
}

template <class T>
std::pair<unsigned int, unsigned int> Graph<T>::getNextHeuristicEdge(std::vector<unsigned int> tour, UFDS tSets) {
    double smallestLength = INF;
    std::pair<unsigned int, unsigned int> edgeExtremities;

    for (auto id: tour) {
        for (int i = 0; i < distanceMatrix[id].size(); i++) {
            if (!tSets.isSameSet(tour[0], i)) {
                if (distanceMatrix[id][i] < smallestLength) {
                    smallestLength = distanceMatrix[id][i];
                    edgeExtremities = {id, i};
                }
            }
        }
    }
    return edgeExtremities;
}



template <class T>
std::pair<std::vector<unsigned int>, double> Graph<T>::getInsertionEdges(std::vector<unsigned int> tour, const unsigned int newVertexId) const {
    std::pair<std::vector<unsigned int>, double> result = {{}, INF};

    for (int i = 0; i < tour.size() - 1; i++) {
        //If there are two edges that could replace the current one, connecting its ends to the new vertex
        double currentDistance = distanceMatrix[tour[i]][newVertexId] + distanceMatrix[newVertexId][tour[i + 1]];
        if (currentDistance < result.second) {
            result.second = currentDistance;
            result.first = {tour[i], newVertexId, tour[i + 1]};
        }
    }
    return result;
}


template <class T>
void Graph<T>::graphClear() {
    distanceMatrix = {};
    vertexSet = {};
    totalEdges = 0;
}

template <class T>
double Graph<T>::getTourDistance() const {
    return tour.distance;
}

inline void deleteMatrix(int **m, int n) {
    if (m != nullptr) {
        for (int i = 0; i < n; i++)
            if (m[i] != nullptr)
                delete [] m[i];
        delete [] m;
    }
}

inline void deleteMatrix(double **m, int n) {
    if (m != nullptr) {
        for (int i = 0; i < n; i++)
            if (m[i] != nullptr)
                delete [] m[i];
        delete [] m;
    }
}

template <class T>
Graph<T>::~Graph() {
    deleteMatrix(distMatrix, vertexSet.size());
    deleteMatrix(pathMatrix, vertexSet.size());
}



#endif /* DA_TP_CLASSES_GRAPH */