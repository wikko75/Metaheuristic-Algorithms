#ifndef ALGORITHMS_1_HPP
#define ALGORITHMS_1_HPP

#include <vector>
#include "Vertex.hpp"
#include <optional>

void localAnnealing(const std::vector<Vertex>& verticies, std::vector<int>& cycle, float T, float deltaT, int epoch, int samples);


/// @brief Computes invert indexes (i and j) of neighbors of cycle. Neighbor in form of std::pair<int, int> having indexes i and j which are used for computing inverse
/// @param cycle initial cycle which neighbors are being computed of 
/// @param verticies vector of Vertex objects containing data regarding verticies 
/// @param neighborsListCapacity number of neightbors that are computed
/// @return indices i and j that form inverse of cycle. Sorted in worse to best order.
std::vector<std::pair<int, int>> generateNeighbors(const std::vector<int>& cycle, const std::vector<Vertex>& verticies, int neighborsListCapacity);


std::optional<std::vector<int>> findBestNeighbor(std::vector<int> initialCycle, const std::vector<Vertex>& verticies,
                                                 std::vector<std::pair<int, int>>& neighbors, std::vector<std::vector<int>>& tabuList);

struct TabuSearchResult
{
    std::vector<int> cycle;
    int cost;
};


//can do better - work all the time on initial vector and dont return copy of it, just cycle cost
//another improvement -- instead of passing vertices and cycle, pass cycle in form of vector<Vertex>
//and use calculateCycle that takes vector<Vertex> as an argument
TabuSearchResult tabuSearch(const std::vector<Vertex> &verticies, std::vector<int>& initialCycle,
                             int tabuListCapacity, int maxRepeats, int neighborsListCapacity);


#endif