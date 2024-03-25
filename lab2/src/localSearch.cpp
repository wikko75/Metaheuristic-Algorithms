#include <iostream>
#include "localSearch.hpp"
#include "utilsFunc.hpp"
#include <fstream>


LocalSearchResult localSearch(const std::vector<Vertex> &vertices, std::vector<int>& initialCycle)
{
    std::vector<int> currentOptimalCycle { initialCycle };
    int currentOptimalCycleCost { calculateCycle(vertices, currentOptimalCycle) };
    bool canImprove { true };
    int noOfInverts {};

    while (canImprove)
    {      
        InvertWeightDiff maxDiff { 0, -1, -1};

        for (int i{0}; i < initialCycle.size() - 1; ++i)
        {
            for (int j{i + 1}; j < initialCycle.size(); ++j)
            {                
                InvertWeightDiff diff { weightDiff(vertices, currentOptimalCycle, i, j) };
                if (diff.cost > maxDiff.cost)
                {
                    maxDiff = diff;
                }
            }
        }

        if (maxDiff.cost <= 0)
        {
            return { currentOptimalCycle, currentOptimalCycleCost, noOfInverts };
        }

        invert(currentOptimalCycle, maxDiff.i, maxDiff.j);
        ++noOfInverts;
        currentOptimalCycleCost -= maxDiff.cost;
    }

    return { currentOptimalCycle, currentOptimalCycleCost, noOfInverts };
}


void saveExperimentData(const std::string &path, const ExperimentData &data)
{   
    std::ofstream outStream{};
    outStream.open(path, std::ios::out);

    if(outStream.is_open())
    {
        outStream << "MST: " << data.MST << "\n";
        outStream << "Best LS cycle weight: " << data.bestLSCycleWeight << "\n";
        outStream << "Average LS cycle weight: " << data.avgLSCycleWeight << "\n";
        outStream << "Average number of inverts: " << data.avgNoOfInverts<< "\n";
    }
    else
    {
        std::cerr << "Can't save experiment data!\nAborting...";
        outStream.close();
    }

    outStream.close();
}


ExperimentData experiment1(Graph& graph, const std::vector<Vertex>& vertices, std::mt19937& mt)
{
    std::uniform_int_distribution randomVertex{ 1, static_cast<int>(vertices.size()) - 1 };

    int noOfInverts {0};
    int sumOfCyclesWeights {0};
    int bestCycleWeight {std::numeric_limits<int>::max()};
    const int SQRT_N { static_cast<int>(ceil(sqrt(vertices.size()))) };

    std::vector<bool> visited(vertices.size(), false);
    std::vector<int> cycle{};
    cycle.reserve(vertices.size());

    {
        Timer timer {"floor(sqrt(n)) timer"};
        for(int i{0}; i < SQRT_N; ++i)
        {
            //Perform DFS
            //-------------------------------
            graph.DFS(randomVertex(mt), visited, cycle);
    
            //Perform Local Search
            //-------------------------------
            LocalSearchResult optimalCycle { localSearch(vertices, cycle) };
            
            sumOfCyclesWeights += optimalCycle.cost;
            noOfInverts += optimalCycle.noOfInverts;

            if (optimalCycle.cost < bestCycleWeight)
            {
                bestCycleWeight = optimalCycle.cost;
            }

            //if I dont assign i get an error cause i dont have reserve before loop
            visited.assign(vertices.size(), false);
            cycle.clear();
            cycle.reserve(vertices.size());
        }
    }

    ExperimentData data { 
        -1,
        bestCycleWeight,
        sumOfCyclesWeights / SQRT_N,
        noOfInverts / SQRT_N 
    };

    return data;
}


ExperimentData experiment2(Graph &graph, std::vector<Vertex> &vertices, std::mt19937 &mt)
{

    std::uniform_int_distribution randomVertex{ 1, static_cast<int>(vertices.size()) - 1 };

    int noOfInverts {0};
    int sumOfCyclesWeights {0};
    int bestCycleWeight { std::numeric_limits<int>::max() };
    std::vector<int> cycle{};
    cycle.reserve(vertices.size());

    {
        Timer timer {"floor(sqrt(n)) timer"};

        for (int i{0}; i < vertices.size(); ++i)
        {
            std::vector<Vertex> randomCycle = permutateVertices(vertices, mt, randomVertex);
            /*
                TODO REFACTOR

                ! bad - because i dont return std::vector<int> from  permutateVertices() but std::vector<Vertex> instead
                ! That's why i have to convert :(
            */

            for(const Vertex& vertex : randomCycle)
            {
                cycle.push_back(vertex.getNumber());
            }

            LocalSearchResult optimalCycle = localSearch(vertices, cycle);

            if (optimalCycle.cost < bestCycleWeight)
            {
                bestCycleWeight = optimalCycle.cost;
            }

            noOfInverts += optimalCycle.noOfInverts;
            sumOfCyclesWeights += optimalCycle.cost;

            cycle.clear();
        }
    }

        ExperimentData data {
            -1,
            bestCycleWeight,
            sumOfCyclesWeights / static_cast<int>(vertices.size()),
            noOfInverts / static_cast<int>(vertices.size())
        };

    return data;
}