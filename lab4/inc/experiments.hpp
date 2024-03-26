#ifndef EXPERIMENTS_2_HPP
#define EXPERIMENTS_2_HPP

#include <filesystem>
#include <vector>
#include <Vertex.hpp>
#include <random>


struct ExperimentData
{
    const char* label;
    int bestCycleWeight{};
    int avgCycleWeight{};
};

void saveExperimentData(const std::filesystem::path& path, const ExperimentData& data);

ExperimentData geneticTest1(std::vector<int>& initialCycle, const int iterations, const int maxIterWithoutImprov,
    const int repeats, const int populationSize, const int noOfPairs, const int sampleSize, const float mutationProb,
    const std::vector<Vertex>& verticies, std::mt19937& mt);

ExperimentData geneticTest2( std::vector<int>& initialCycle, const int iterations, const int maxIterWithoutImprov,
    const int repeats, const int populationSize, const int noOfPairs,
    const float mutationProb, const std::vector<Vertex>& verticies, std::mt19937& mt);

#endif