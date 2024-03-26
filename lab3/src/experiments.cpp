#include "experiments.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include "algorithms.hpp"
#include "utilsFunc.hpp"


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


ExperimentData simulatedAnnealingExperiment(const std::vector<Vertex> &vertices, std::vector<int> &cycle, const int iterations, std::mt19937& mt)
{
    int avgWeight {};
    int min { std::numeric_limits<int>::max() };

    for (int i {0}; i < iterations; ++i)
    {
    
        std::shuffle(cycle.begin(), cycle.end(), mt);

        //Local annealing 
        localAnnealing(vertices, cycle, 100, 0.85f, 4000, 2500);

        int SAWeight { calculateCycle(vertices, cycle) };
        std::cout << "Cycle weight for local annealing: " << SAWeight << "\n";

        if (SAWeight < min)
        {
            min = SAWeight;
        }
        
        avgWeight += SAWeight;
    }

    ExperimentData SAData { "SA", min, avgWeight / iterations };

    return SAData;
}



ExperimentData tabuSearchExperiment(const std::vector<Vertex>& vertices, std::vector<int>& cycle, const int iterations, std::mt19937& mt)
{

    int avgWeight {};
    int min { std::numeric_limits<int>::max() }; 
    for (int  i {0}; i < iterations; ++i)
    {
        std::shuffle(cycle.begin(), cycle.end(), mt);
        TabuSearchResult TS {tabuSearch(vertices, cycle, pow(vertices.size(), 2), vertices.size(), 100)};

        if (TS.cost < min)
        {
            min = TS.cost;
        }

        avgWeight += TS.cost;
    }

    ExperimentData TSData {"TS", min, avgWeight / iterations};

    return TSData;
}