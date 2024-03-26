#include "experiments.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <fmt/core.h>
#include "geneticTSPSolver.hpp"

ExperimentData geneticTest1(std::vector<int> &initialCycle, const int iterations, const int maxIterWithoutImprov, const int repeats, const int populationSize, const int noOfPairs, const int sampleSize, const float mutationProb, const std::vector<Vertex> &verticies, std::mt19937 &mt)
{
    fmt::print("\n Starting experiment...\n");
    int minCost {std::numeric_limits<int>::max()};
    int avgCost {};

    for (int i {0}; i < repeats; ++i)
    {

        int currMinCost { geneticTSPSolver(initialCycle, iterations, maxIterWithoutImprov, populationSize, noOfPairs, sampleSize, mutationProb, verticies, mt) };
        if (currMinCost < minCost)
        {
            minCost = currMinCost;
        }

        avgCost += currMinCost;
    }

    fmt::print("\nResults after {} repeats:\nMinimal cost: {}\nAvg. cost {}\n", repeats, minCost, avgCost/repeats);

    const ExperimentData data {"GenTSP", minCost, avgCost/repeats};

    return data;

}

ExperimentData geneticTest2(std::vector<int> &initialCycle, const int iterations, const int maxIterWithoutImprov, const int repeats,
                            const int populationSize, const int noOfPairs, const float mutationProb, const std::vector<Vertex> &verticies, std::mt19937 &mt)
{

    int minCost {std::numeric_limits<int>::max()};
    int avgCost {};

    for (int i {0}; i < repeats; ++i)
    {

        int currMinCost { geneticTSPSolver(initialCycle, iterations, maxIterWithoutImprov, populationSize, noOfPairs, mutationProb, verticies, mt) };
        if (currMinCost < minCost)
        {
            minCost = currMinCost;
        }

        avgCost += currMinCost;
    }

    fmt::print("\nResults after {} repeats:\nMinimal cost: {}\nAvg. cost {}\n", repeats, minCost, avgCost/repeats);

    const ExperimentData data {"GenTSP", minCost, avgCost/repeats};

    return data;
}

void saveExperimentData(const std::filesystem::path& path, const ExperimentData& data)
{   
    std::ofstream outStream{};
    outStream.open(path, std::ios::out);

    if(outStream.is_open())
    {
        outStream << "Best " + std::string(data.label) + " cycle weight: " << data.bestCycleWeight << "\n";
        outStream << "Average " + std::string(data.label) + " cycle weight: " << data.avgCycleWeight << "\n";
    }
    else
    {
        std::cerr << "Can't save experiment data!\nAborting...";
        outStream.close();
    }

    outStream.close();
}