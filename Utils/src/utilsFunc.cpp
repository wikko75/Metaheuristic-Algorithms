#include "utilsFunc.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include "Graph.hpp"


int euclideanNorm(const Vertex& v1, const Vertex& v2) 
{
    double dx = v1.getX() - v2.getX();
    double dy = v1.getY() - v2.getY();
    double distance = std::sqrt(dx * dx + dy * dy);
    
    return static_cast<int>(std::round(distance));
}


void getGraphDataFromFile(const std::filesystem::path& filePath, std::vector<Vertex>& vertices, int startLine)
{
    std::fstream fileStream{};
    fileStream.open(filePath, std::ios::in);

    if (!fileStream.is_open())
    {
        std::cerr << "Can't open a file\n";
        return;
    }

    int currentLineCount{0};
    std::string line{};
    
    while (!fileStream.eof())
    {
        ++currentLineCount;
        std::getline(fileStream, line);

        if (line == "EOF")
        {
            fileStream.close();
            return;
        }
        
        if (currentLineCount >= startLine)
        {
            std::stringstream strStream(line);
            int vertexNumber{};
            int x{};
            int y{};

            if (strStream >> vertexNumber >> x >> y)
            {
                vertices.emplace_back(vertexNumber, x, y);
            }
            else
            {
                std::cerr << "Can't convert line to numbers :(\n";
            }
        }
    }

    fileStream.close();   
}


void saveMST(const std::string& path, const Graph& graph)
{
    std::ofstream ostream{};
    std::cout << path << "\n\n";
    ostream.open(path, std::ios::out);

    if (!ostream.is_open())
    {
        std::cerr << "Can't write MST edges to file!\nAborting..." << std::endl;
        ostream.close();
        return;
    }

    std::cout << "Saving MST edges...\n";
    for (const Edge &edge : graph.getMSTEdges())
    {
        ostream <<  edge.getSrc() << " " << edge.getDest() << "\n";
    }

    ostream.close();
}


void saveDFSCycle(const std::string& path, const std::vector<int>& cycle)
{
    std::ofstream fstream {};

    fstream.open(path, std::ios::out);

    std::cout << "Path to DFS edges: " <<  path << "\n";

    if(!fstream.is_open())
    {
        std::cerr << "Can't write DFS cycle edges\nAborting...\n";
        fstream.close();
        return;
    }

    std::cout << "Saving DFS edges...\n";
    for(int i {0}; i < cycle.size() - 1; ++i)
    {
        
        fstream << cycle[i] << " " << cycle[i + 1] << "\n";   
        
    }

    fstream << cycle[cycle.size() - 1] << " " << cycle[0] << "\n";

    fstream.close();
}


int calculateCycle(const std::vector<Vertex>& vertices, const std::vector<int>& cycle)
{
    int cycleWeight {};

    //add all weights between first and last element in cycle
    for(int index{0}; index < cycle.size() - 1;  ++index)
    {
        cycleWeight += euclideanNorm(vertices[cycle[index] - 1], vertices[cycle[index + 1] - 1]);
    }

    //add weight between last and first element in cycle (cause we have to come back, since it's a cycle)
    cycleWeight += euclideanNorm(vertices[cycle[0] - 1], vertices[cycle[cycle.size() - 1] - 1]);

    return cycleWeight;
}


int calculatePermutationCycle(const std::vector<Vertex>& permutation)
{
    int cycleWeight {};

    //add all weights between first and last element in cycle
    for(int index{0}; index < permutation.size() - 1;  ++index)
    {
        cycleWeight += euclideanNorm(permutation[index], permutation[index + 1]);
    }

    //add weight between last and first element in cycle (cause we have to come back, since it's a cycle)
    cycleWeight += euclideanNorm(permutation[permutation.size() - 1], permutation[0]);

    return cycleWeight;
}


std::vector<Vertex> verticesPermutation(std::vector<Vertex>& vertices, std::mt19937& mt, std::uniform_int_distribution<int>& randomVertex) {
    std::vector<Vertex> permutation = vertices;  // Creating a copy to generate permutations

    // Use a different seed for better randomness
    std::mt19937 gen(mt());
    
    // Rearrange the elements in the vector to create permutations
    for (int i = static_cast<int>(permutation.size()) - 1; i >= 1; --i) {
        int j = randomVertex(gen, std::uniform_int_distribution<int>::param_type(0, i));  // Generating a random index directly from the distribution

        // Swap elements to create a permutation
        std::swap(permutation[i], permutation[j]);
    }

    return permutation;
}