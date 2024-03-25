#include <iostream>
#include <algorithm>
#include "Graph.hpp"
#include "utilsFunc.hpp"


struct Timer
{

    Timer(const std::string& name)
    {   
        this->name = name;
        start = std::chrono::high_resolution_clock::now();
        std::cout << "["<< this->name  << "] " << "Timer started!\n";
    }

    ~Timer()
    {
        end = std::chrono::high_resolution_clock::now();
        std::cout << "[" << this->name << "]" <<  "ended!\n";
        std::cout << "Elapsed time: " << std::chrono::duration<float, std::milli>(end - start).count() << std::endl;
    }

private:
    std::chrono::_V2::system_clock::time_point start;
    std::chrono::_V2::system_clock::time_point end;
    std::string name;
};



int test(std::vector<Vertex>& vertices, std::mt19937& mt, std::uniform_int_distribution<int>& randomVertex, int noOfGroups, int noOfRepeats, int& minValue)
{
    double avg {};
    std::vector<int> minValueFromGroup {};
    minValueFromGroup.reserve(noOfGroups);

    for (int i {0}; i < noOfRepeats; ++i)
    {
        std::vector<int> group{};
        group.reserve(noOfGroups);
        for (int j {0}; j < noOfGroups; ++j)
        {
            std::vector<Vertex> permutation = verticesPermutation(vertices, mt, randomVertex);
            int permutationCycleWeight = calculatePermutationCycle(permutation);
            group.emplace_back(permutationCycleWeight);
        }

        std::sort(group.begin(), group.end(), [](int a, int b){
             return a < b;
             });

        avg += group[0];
        minValueFromGroup.emplace_back(group[0]);
    }

    std::sort(minValueFromGroup.begin(), minValueFromGroup.end(), [](int a, int b){
             return a < b;
             });
    
    minValue = minValueFromGroup[0];
    return static_cast<int>(std::round(avg/noOfRepeats));
}



int main()
{
    //path to files containing vertices data in format (no. of vertex, x coor, y coor)
    const std::string pathToVerticesData {std::filesystem::current_path() / "test" / "verticesData"};

    for (auto &verticesDataFile : std::filesystem::directory_iterator(pathToVerticesData))
    {       
        std::cout << "==========================\n\n";
        std::cout << verticesDataFile.path() << "\n";
        std::vector<Vertex> vertices {};

        // populate vertices with data and create graph
        getGraphDataFromFile(verticesDataFile.path(), vertices);
        
        std::unique_ptr<Graph> graph { std::make_unique<Graph>(vertices.size()) };
        graph->createCompleteGraph(vertices);

        // Create MST 
        // -------------------------------
        const int MST { graph->kruskalMST() };
        std::cout << "MST = " << MST  << "\n";

        // Save MST edges
        // -------------------------------
        std::string pathToMSTEdges {std::filesystem::current_path() / "res" / "MSTEdges/"};
        pathToMSTEdges.append(verticesDataFile.path().filename().string());
        saveMST(pathToMSTEdges, *graph);

        // Perform DFS
        // -------------------------------
        std::vector<bool> visited(vertices.size(), false);
        std::vector<int> cycle{};
        cycle.reserve(vertices.size());

        // Seed our Mersenne Twister using steady_clock
        std::mt19937 mt{ static_cast<std::mt19937::result_type>(
            std::chrono::steady_clock::now().time_since_epoch().count())};

        // Create a reusable random number generator that generates random vertices from graph 
        std::uniform_int_distribution randomVertex{ 1, static_cast<int>(vertices.size()) - 1 }; 
            
        graph->DFS(randomVertex(mt), visited, cycle);
        int cycleWeight = calculateCycle(vertices, cycle);
        std::cout << "Cycle weight for DFS: " << cycleWeight << "\n";

        // Save Cycle
        // -------------------------------
        std::string pathToDFSCycleEdges {std::filesystem::current_path() / "res" / "DFSCycleEdges/"};
        pathToDFSCycleEdges.append(verticesDataFile.path().filename().string());
        saveDFSCycle(pathToDFSCycleEdges, cycle);


        // TESTS
        // -------------------------------
        int MIN_VALUE_10_GROUPS {};
        int MIN_VALUE_20_GROUPS {};
        std::cout << "TEST 1 [100 repeats] [10 groups] Average min " << test(vertices, mt, randomVertex, 10, 100, MIN_VALUE_10_GROUPS) << "\n"; 
        std::cout << "TEST 2 [50 repeats] [20 groups] Average min " << test(vertices, mt, randomVertex, 20, 50, MIN_VALUE_20_GROUPS) << "\n"; 

        std::cout << "MIN_VALUE_10_GROUPS: " << MIN_VALUE_10_GROUPS << "\n";
        std::cout << "MIN_VALUE_20_GROUPS: " << MIN_VALUE_20_GROUPS << "\n\n\n";

        std::cout << "==========================\n\n";

    }  

    return 0;
}


// TODO - refactor of DFS (put it in Graph.cpp) DONE
// TODO - put all utility functions (saving etc.) into designated folder  DONE
// TODO - euclidian norm (just taking int, not Vertex) not really reasonable to do
// TODO - createCompleteGraph also part of Graph.cpp  DONE