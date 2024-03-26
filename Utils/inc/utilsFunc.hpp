#ifndef UTILS_FUNC_1_HPP
#define UTILS_FUNC_1_HPP

#include "Vertex.hpp"
#include <iostream>
#include <filesystem>
#include <random> // for std::mt19937
#include <chrono> // for std::chrono
#include "Graph.hpp"



/// @brief Calculates norm between v1 and v2
/// @param v1 
/// @param v2 
/// @return norm value
int euclideanNorm(const Vertex& v1, const Vertex& v2);


/// @brief Loads graph data (vertices' properties) from filePath and stores it in vertices vector
/// @param filePath path to file where to look for data
/// @param vertices vector to save data to (out parameter)
/// @param startLine start reading data from this line (inclusive), default = 9
void getGraphDataFromFile(const std::filesystem::path& filePath, std::vector<Vertex>& vertices, int startLine = 9);


/// @brief 
/// @param path 
/// @param graph 
void saveMST(const std::string& path, const Graph& graph);


/// @brief Saves cycle edges to provided file
/// @param path path to file, if file does not exists it will be created
/// @param cycle vector storing vertices' numbers that form cycle
void saveCycle(const std::string& path, const std::vector<int>& cycle);


/// @brief Computes weigth of cycle formed by DFS
/// @param vertices vector with data describing each vertex properties
/// @param cycle cycle which weight is computed (out parameter)
/// @return  weight of  cycle
int calculateCycle(const std::vector<Vertex>& vertices, const std::vector<int>& cycle);


/// @brief Calculates weight of cycle formed by permutation of vertices
/// @param permutation vector with data required for weight calculation
/// @return weight of cycle
int calculatePermutationCycle(const std::vector<Vertex>& permutation);


/// @brief Permutates provided vertices using random generators
/// @param vertices vector of Vertex objects based on which permutation will be done
/// @param mt random generator
/// @param randomVertex uniform number generator that generates random vertices from graph  
/// @return vector of Vertex objects after permutation
std::vector<Vertex> permutateVertices(const std::vector<Vertex>& vertices, std::mt19937& mt, std::uniform_int_distribution<int>& randomVertex);


/// @brief Performs invert operations on provided cycle between i and j (inclusive)
/// @param cycle 
/// @param i 
/// @param j 
void invert(std::vector<int>& cycle, int i, int j);


struct InvertWeightDiff
{
    int cost{};
    int i{};
    int j{};
};

//possible change -- instead of providing std::vector<Vertex> &vertices, std::vector<int>& cycle, just provide std::vector<Vertex>
//and use calculatePermutationCycle

/// @brief Smart way of computation difference in weights after performing invert
/// @param vertices 
/// @param cycle 
/// @param startIdx 
/// @param stopIdx 
/// @return 
InvertWeightDiff weightDiff(const std::vector<Vertex> &vertices, const std::vector<int>& cycle, int startIdx, int stopIdx);


struct Timer
{
    Timer(const std::string& name): name{name}
    {   
        start = std::chrono::high_resolution_clock::now();
        std::cout << "["<< this->name  << "] " << " Timer started!\n";
    }

    ~Timer()
    {
        end = std::chrono::high_resolution_clock::now();
        std::cout << "[" << this->name << "]" <<  " Timer ended!\n";
        std::cout << "Elapsed time: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;
    }

private:
    std::chrono::_V2::system_clock::time_point start;
    std::chrono::_V2::system_clock::time_point end;
    std::string name;
};


#endif