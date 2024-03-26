#ifndef LOCAL_SEARCH_1_HPP
#define LOCAL_SEARCH_1_HPP

#include <random>
#include <vector>
#include <string>
#include "Graph.hpp"

struct LocalSearchResult
{
    std::vector<int> cycle;
    int cost;
    int noOfInverts;
};

//can do better - work all the time on initial vector and dont return copy of it, just cycle cost
//another improvement -- instead of passing vertices and cycle, pass cycle in form of vector<Vertex>
//and use calculateCycle that takes vector<Vertex> as an argument
LocalSearchResult localSearch(const std::vector<Vertex> &vertices, std::vector<int>& initialCycle);

struct ExperimentData
{
    int MST{};
    int bestLSCycleWeight{};
    int avgLSCycleWeight{};
    int avgNoOfInverts{};
};

/// @brief Local Search done floor(sqrt(n)) times for each set of data with DFS as initial cycle
ExperimentData experiment1(Graph& graph, const std::vector<Vertex>& vertices, std::mt19937& mt);

/// @brief Experiment 2 - Local Search for n random permutations of each data set 
ExperimentData experiment2(Graph& graph, std::vector<Vertex>& vertices, std::mt19937& mt);

void saveExperimentData(const std::string& path, const ExperimentData& data);



#endif