#ifndef UTILS_FUNC_1_HPP
#define UTILS_FUNC_1_HPP

#include "Vertex.hpp"
#include <filesystem>
#include <random> // for std::mt19937
#include <chrono> // for std::chrono
#include "Graph.hpp"


int euclideanNorm(const Vertex& v1, const Vertex& v2);

void getGraphDataFromFile(const std::filesystem::path& filePath, std::vector<Vertex>& vertices, int startLine = 9);

void saveMST(const std::string& path, const Graph& graph);

void saveDFSCycle(const std::string& path, const std::vector<int>& cycle);

//computes weigth of cycle formed by DFS
int calculateCycle(const std::vector<Vertex>& vertices, const std::vector<int>& cycle);

int calculatePermutationCycle(const std::vector<Vertex>& permutation);

std::vector<Vertex> verticesPermutation(std::vector<Vertex>& vertices, std::mt19937& mt, std::uniform_int_distribution<int>& randomVertex);


#endif