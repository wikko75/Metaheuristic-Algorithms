#include "algorithms.hpp"
#include <random>
#include <chrono>
#include <algorithm>
#include "utilsFunc.hpp"


void localAnnealing(const std::vector<Vertex>& verticies, std::vector<int>& cycle, float T, float deltaT, int epoch, int samples)
{

    std::mt19937 mt{ static_cast<std::mt19937::result_type>(
    std::chrono::steady_clock::now().time_since_epoch().count()
    )};

    std::uniform_int_distribution randomVertex { 1, static_cast<int>(verticies.size())};
    std::uniform_real_distribution<float> uniformRandomFloat {0.0, 1.0};

    int i {};
    int j {};
    int temp {};
    float randomFloat {};
    int deltaWeight {};

    while (T > 0.0 && epoch > 0)
    {
        for(int sampleCounter {0};  sampleCounter < samples; ++sampleCounter)
        {
            i = randomVertex(mt);
            do
            {
                j = randomVertex(mt);
            }
            while (i == j);

            InvertWeightDiff diff = weightDiff(verticies, cycle, i - 1, j - 1);
            deltaWeight = diff.cost;

            if (deltaWeight > 0)
            {
               invert(cycle, i - 1, j - 1);
            }

            else
            {
                randomFloat = uniformRandomFloat(mt);

                if(randomFloat < exp(deltaWeight / T))
                {
                    invert(cycle, i - 1, j - 1);
                }
            }
        }

        T *= deltaT;
        --epoch;
    }    
}


// @brief Computes invert indexes (i and j) of neighbors of cycle. Neighbor in form of std::pair<int, int> having indexes i and j which are used for computing inverse
/// @param cycle initial cycle which neighbors are being computed of 
/// @param verticies vector of Vertex objects containing data regarding verticies 
/// @param neighborsListCapacity number of neightbors that are computed
/// @return indices i and j that form inverse of cycle. Sorted in worse to best order.
std::vector<std::pair<int, int>> generateNeighbors(
    const std::vector<int>& cycle,
    const std::vector<Vertex>& verticies,
    int neighborsListCapacity
    )
{
    std::mt19937 mt{ static_cast<std::mt19937::result_type>(
    std::chrono::steady_clock::now().time_since_epoch().count()
    )};

    std::uniform_int_distribution randomVertexIndex { 0, static_cast<int>(verticies.size()) - 1};

    int i {};
    int j {};
    InvertWeightDiff maxDiff { 0, -1, -1};
    std::vector<std::pair<int, int>> neighborsInvertIndices{};

    for (int k {0}; k < neighborsListCapacity; ++k)
    {
        i = randomVertexIndex(mt);
        do
        {
            j = randomVertexIndex(mt);
        } while (i == j);


        InvertWeightDiff diff { weightDiff(verticies, cycle, i, j) };

        if (diff.cost > maxDiff.cost)
        {
            maxDiff = diff;
            neighborsInvertIndices.emplace_back(maxDiff.i, maxDiff.j);
        }
    }

    return neighborsInvertIndices;
}


std::optional<std::vector<int>> findBestNeighbor(
    std::vector<int> initialCycle,
    const std::vector<Vertex>& verticies,
    std::vector<std::pair<int, int>>& neighbors,
    std::vector<std::vector<int>>& tabuList
    )
{
    while (neighbors.size() != 0)
    {
        std::vector<int> cycle {initialCycle};
        auto [i, j] {neighbors.back()};
        invert(cycle, i, j);

        auto it = std::find_if(tabuList.begin(), tabuList.end(), [&cycle](const std::vector<int>& tabu)
        {
            return std::equal(tabu.begin(), tabu.end(), cycle.begin(), cycle.end());
        });


        if ( it == tabuList.end() || tabuList.size() == 0)
        {
            return cycle;
        }

        neighbors.pop_back();
    }
    //all in tabu - return nothing
    return {};
}


TabuSearchResult tabuSearch(const std::vector<Vertex> &verticies, std::vector<int>& initialCycle,
                             int tabuListCapacity, int maxRepeats, int neighborsListCapacity)
{
    std::vector<int> currentOptimalCycle { initialCycle };
    int repeats {};

    std::vector<std::vector<int>> tabuList{};
    tabuList.reserve(tabuListCapacity);

    while (tabuList.size() < tabuListCapacity && repeats < maxRepeats)
    {      
        std::vector<std::pair<int, int>> neighbors { generateNeighbors(currentOptimalCycle, verticies, neighborsListCapacity) };
        auto bestNeighbor { findBestNeighbor(currentOptimalCycle, verticies, neighbors, tabuList) };

        //no bestNeightbor found (all in tabu list)
        if (!bestNeighbor.has_value())
        {
            ++repeats;
        }
        else
        {
            currentOptimalCycle = bestNeighbor.value();
            tabuList.emplace_back(bestNeighbor.value());
            repeats = 0;
        }
    }

    return { currentOptimalCycle, calculateCycle(verticies, currentOptimalCycle) };
}