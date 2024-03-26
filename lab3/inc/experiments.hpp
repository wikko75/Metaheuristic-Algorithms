#ifndef EXPERIMENTS_1_HPP
#define EXPERIMENTS_1_HPP

#include <filesystem>
#include <random>
#include <vector>
#include "Vertex.hpp"

struct ExperimentData
{
    const char* label;
    int bestCycleWeight{};
    int avgCycleWeight{};
};

void saveExperimentData(const std::filesystem::path& path, const ExperimentData& data);

ExperimentData simulatedAnnealingExperiment(const std::vector<Vertex>& vertices, std::vector<int>& cycle, const int iterations, std::mt19937& mt);

ExperimentData tabuSearchExperiment(const std::vector<Vertex>& vertices, std::vector<int>& cycle, const int iterations, std::mt19937& mt);

#endif