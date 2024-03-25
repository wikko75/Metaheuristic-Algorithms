#include <iostream>
#include <algorithm>
#include <filesystem>
#include "utilsFunc.hpp"
#include "localSearch.hpp"



int main() 
{
    //path to files containing vertices data in format (no. of vertex, x coor, y coor)
    const std::string pathToVerticesData { std::filesystem::current_path() / "test" / "verticesData" };
    
    {
        Timer outer {"Verticies >= 1000"};


        for (auto &verticesDataFile : std::filesystem::directory_iterator(pathToVerticesData))
        {       
            std::cout << "==========================\n\n";
            std::cout << verticesDataFile.path() << "\n";
            std::vector<Vertex> vertices {};

            // populate vertices with data
            getGraphDataFromFile(verticesDataFile.path(), vertices);
            std::unique_ptr<Graph> graph { std::make_unique<Graph>(vertices.size()) };
            graph->createCompleteGraph(vertices);

            // Create MST 
            // -------------------------------
            const int MST { graph->kruskalMST() };

            std::cout << "MST = " << MST  << "\n";

            // Save MST edges data
            // -------------------------------
            std::string pathToMSTEdges { std::filesystem::current_path() / "test" / "MSTEdges/" };
            pathToMSTEdges.append(verticesDataFile.path().filename().string());
            saveMST(pathToMSTEdges, *graph);


            // Seed our Mersenne Twister using steady_clock
            std::mt19937 mt{ static_cast<std::mt19937::result_type>(
                std::chrono::steady_clock::now().time_since_epoch().count()
                )};

            // Create a reusable random number generator that generates random vertices from graph 
            std::uniform_int_distribution randomVertex{ 1, static_cast<int>(vertices.size()) - 1 };


            // Initial (DFS and Local Search done once for each set of data, results stored in respective /res/ folders)
            // --------------------------------------------------------------------------------------
            

            // Perform DFS
            // -------------------------------
            std::vector<bool> visited(vertices.size(), false);
            std::vector<int> cycle{};
            cycle.reserve(vertices.size());

            graph->DFS(randomVertex(mt), visited, cycle);                
            int cycleWeight = calculateCycle(vertices, cycle);
            std::cout << "Cycle weight for DFS: " << cycleWeight << "\n";


            // save DFS Cycle
            // ------------------------------
            std::string pathToDFSCycleEdges { std::filesystem::current_path() / "test" / "DFSCycleEdges/" };
            pathToDFSCycleEdges.append(verticesDataFile.path().filename().string());
            saveCycle(pathToDFSCycleEdges, cycle);


            // Perform Local Search
            // -------------------------------
            LocalSearchResult optimalCycle = localSearch(vertices, cycle);
            std::cout << "Optimal cycle weight = " << optimalCycle.cost << "\n";

            // Save Local Search cycle 
            std::string pathToLocalSearchCycleEdges { std::filesystem::current_path() / "test" / "LocalSearchCycleEdges/" };
            pathToLocalSearchCycleEdges.append(verticesDataFile.path().filename().string());
            saveCycle(pathToLocalSearchCycleEdges, optimalCycle.cycle);


            // --------------------------------------------------------------------------------------


            // Experiment 1 - Local Search done floor(sqrt(n)) times for each set of data with DFS as initial cycle
            // --------------------------------------------------------------------------------------
            
            ExperimentData data { experiment1(*graph, vertices, mt) };
            data.MST = MST;

            std::string experiment1DataPath { std::filesystem::current_path() / "test" / "experiment1Data/" };
            experiment1DataPath.append("EX1_DATA_").append(verticesDataFile.path().filename().string());

            saveExperimentData(experiment1DataPath, data);

            std::cout << "MST: " << MST << "\n";
            std::cout << "Best LS cycle weight: " << data.bestLSCycleWeight << "\n";
            std::cout << "Average LS cycle weight: " << data.avgLSCycleWeight << "\n";
            std::cout << "Average number of inverts: " << data.avgNoOfInverts << "\n";

            // --------------------------------------------------------------------------------------


            // Experiment 2 - Local Search for n random permutations of each data set
            // --------------------------------------------------------------------------------------

            ExperimentData data2 { experiment2(*graph, vertices, mt) };
            std::string path { std::filesystem::current_path() / "res" / "experiment2Data/" };
            path.append("EX2_DATA_").append(verticesDataFile.path().filename().string());
            saveExperimentData(path, data);
            
            std::cout << "Best LS cycle weight: " << data2.bestLSCycleWeight << "\n";
            std::cout << "Average LS cycle weight: " << data2.bestLSCycleWeight << "\n";
            std::cout << "Average number of inverts: " << data2.avgNoOfInverts << "\n";
        }
    }
    return 0;
}