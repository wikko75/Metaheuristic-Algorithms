#include <iostream>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <random> // for std::mt19937
#include <limits>
#include <fmt/core.h>
#include <fmt/ranges.h>

#include "utilsFunc.hpp"
#include "experiments.hpp"
#include "geneticTSPSolver.hpp"




int main() {
    //path to files containing verticies data in format (vertex no., x coor, y coor)
    const std::filesystem::path pathToverticiesData {std::filesystem::current_path() / "test" /  "verticiesData" };

    //random gen. used throughout program
    std::mt19937 mt {std::random_device{}()};

    constexpr int iter {20000};
    constexpr int maxIterWithoutImprov {1000};
    constexpr int populationSize {50};
    constexpr int noOfPairs {20};
    constexpr float mutationProb {.15};
    constexpr int populationSample {10};

    constexpr int  repeats {1};

    {
        for (auto &verticiesDataFile : std::filesystem::directory_iterator(pathToverticiesData))
        {       
            fmt::print("Path: {}\n", verticiesDataFile.path().string());
            
            std::vector<Vertex> verticies {};

            //populate verticies with data
            getGraphDataFromFile(verticiesDataFile.path(), verticies);
            
            std::vector<int> cycle{};
            cycle.reserve(verticies.size());

            for (int i {1}; i < static_cast<int>(verticies.size() + 1); ++i)
            {
                cycle.push_back(i);
            }

            const ExperimentData data { geneticTest2(cycle, iter, maxIterWithoutImprov, repeats, populationSize, noOfPairs, mutationProb, verticies, mt) };
           
            fmt::print("\nResults after {} repeats:\nMinimal cost: {}\nAvg. cost {}\n", repeats, data.bestCycleWeight, data.avgCycleWeight);

            const std::filesystem::path path {std::filesystem::current_path().append("res/resultsOPC")
                                            .append("GenTSP_DATA" + verticiesDataFile.path().filename().string())};

            fmt::print("Save path: {}\n", path.string());
            saveExperimentData(path, data);

            const ExperimentData data2 { geneticTest1(cycle, iter, maxIterWithoutImprov, repeats, populationSize, noOfPairs, populationSample, mutationProb, verticies, mt) };
  
            fmt::print("\nResults after {} repeats:\nMinimal cost: {}\nAvg. cost {}\n", repeats, data2.bestCycleWeight, data2.avgCycleWeight);

            const std::filesystem::path path2 {std::filesystem::current_path().append("res/resultsTPC")
                                            .append("GenTSP_DATA" + verticiesDataFile.path().filename().string())};

            fmt::print("Save path: {}\n", path2.string());
            saveExperimentData(path2, data2);
        }
    }

    return 0;
}