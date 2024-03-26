#include <iostream>
#include <vector>
#include <random> // for std::mt19937
#include <chrono> // for std::chrono
#include "utilsFunc.hpp"
#include "experiments.hpp"
#include "algorithms.hpp"




int main() {
    //path to files containing vertices data in format (no. of vertex, x coor, y coor)
    const std::filesystem::path pathToVerticesData {std::filesystem::current_path() / "res" / "verticesData"};

    // Seed our Mersenne Twister using steady_clock
    std::mt19937 mt { static_cast<std::mt19937::result_type>(
        std::chrono::steady_clock::now().time_since_epoch().count()
        )};

    {
        Timer outer {"Verticies < 1000"};


        for (auto &verticesDataFile : std::filesystem::directory_iterator(pathToVerticesData))
        {       
            std::cout << verticesDataFile.path() << "\n";
            std::vector<Vertex> vertices {};

            //populate vertices with SAData
            getGraphDataFromFile(verticesDataFile.path(), vertices);
            
            std::vector<int> cycle{};
            cycle.reserve(vertices.size());

            for (int i {1}; i < vertices.size() + 1; ++i)
            {
                cycle.push_back(i);
            }
            
            //------ SIMULATED ANNEALING -----//

            const ExperimentData SAData  { simulatedAnnealingExperiment(vertices, cycle, 1, mt) };
           

            std::filesystem::path path { std::filesystem::current_path() / "res" / "SimulatedAnnealingData/" };
            path.append("SA_DATA_" + verticesDataFile.path().filename().string());

            std::cout << path << "\n";
            saveExperimentData(path, SAData);
            
            //---------------------------//



            //--------TABU SEARCH--------//

            const ExperimentData TSData { tabuSearchExperiment(vertices, cycle, 1, mt) };
           
            path = std::filesystem::current_path() / "res" / "TabuSearchData/";
            path.append("TS_DATA_" + verticesDataFile.path().filename().string());

            std::cout << "Tabu Search cost (best): " << TSData.bestCycleWeight << "\n";
            std::cout << "Tabu Search cost (avg): " << TSData.avgCycleWeight << "\n";

            saveExperimentData(path, TSData);
        }
    }

    return 0;
}