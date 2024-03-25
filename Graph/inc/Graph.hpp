#ifndef GRAPH_1_HPP
#define GRAPH_1_HPP

#include <vector>
#include "Edge.hpp"
#include "Vertex.hpp"

class Graph
{  
public:
    Graph(int V);

    void createCompleteGraph(std::vector<Vertex>& vertices);

    void addEdge(int src, int dest, int weight); 

    int kruskalMST();

    void DFS(int v, std::vector<bool>& visited, std::vector<int>& cycle);
    
    int getNoOfGraphEdges() const;

    std::vector<Edge> getGraphEdges() const;

    std::vector<Edge> getMSTEdges() const;

    ~Graph() = default;

private:
    int V;
    std::vector<Edge> graphEdges;
    std::vector<Edge> MSTEdges;
};

#endif