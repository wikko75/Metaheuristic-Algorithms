#include <algorithm>
#include <memory>
#include "Graph.hpp"
#include "DisjoinSets.hpp"
#include "utilsFunc.hpp"

Graph::Graph(int V) : V {V}, graphEdges {}, MSTEdges {}
{
    graphEdges.reserve(V / 2 * (V - 1));
}


void Graph::createCompleteGraph(std::vector<Vertex> &vertices)
{    
    //construct complete graph with weight as Euclidian norm
    for (int i{0}; i < vertices.size() - 1; ++i)
    {
        Vertex vertex1 = vertices[i];
        for (int j {i + 1}; j < vertices.size(); ++j)
        {
            Vertex vertex2 = vertices[j];
            int weight = euclideanNorm(vertex1, vertex2);

            //add edge between vertex1 and vertex2 with proper weight 
            this->addEdge(vertex1.getNumber(), vertex2.getNumber(), weight);
        }
    }
}


void Graph::addEdge(int src, int dest, int weight) 
{
    graphEdges.emplace_back(src, dest, weight);
}


int Graph::kruskalMST() 
{
    int mst_wt {0};

    std::sort(graphEdges.begin(), graphEdges.end(), [](const Edge& a, const Edge& b) {
        return a.getWeight() < b.getWeight();
    });

    std::unique_ptr<DisjoinSets> sets = std::make_unique<DisjoinSets>(V);

    for (const Edge& edge : graphEdges) 
    {
        int u = edge.getSrc();
        int v = edge.getDest();

        //find sets, where vertex u/v belongs to
        int set_u = sets->find(u);
        int set_v = sets->find(v);

        if (set_u != set_v) 
        {
            // std::cout << u << " - " << v << "\t\tweight: " << edge.getWeight() <<  std::endl;
            mst_wt += edge.getWeight();
            sets->merge(set_u, set_v);
            MSTEdges.push_back(edge);
        }
    }
    return mst_wt;
}


void Graph::DFS(int v, std::vector<bool> &visited, std::vector<int> &cycle)
{
    visited[v - 1] = true;
    cycle.emplace_back(v);

    // look for edges in MSTEdges which are conneted to v
    for (const Edge& edge : MSTEdges) {
        if (edge.hasVertex(v))
        {
            int nextVertex = edge.getOppositeVertex(v);
            if (!visited[nextVertex - 1])
            {
                DFS(nextVertex, visited, cycle);
            }
        }
    } 
}


int Graph::getNoOfGraphEdges() const 
{
    return this->graphEdges.size();
}


std::vector<Edge> Graph::getGraphEdges() const 
{
    return graphEdges;
}


std::vector<Edge> Graph::getMSTEdges() const
{
    return MSTEdges;
}

