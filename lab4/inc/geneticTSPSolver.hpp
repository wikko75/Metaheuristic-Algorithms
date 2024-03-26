#ifndef GENETIC_TSP_SOLVER_1_HPP
#define GENETIC_TSP_SOLVER_1_HPP


#include <vector>
#include <random>
#include "Vertex.hpp"


std::vector<std::vector<int>> createPopulation(std::vector<int>& cycle, int populationSize, std::mt19937& mt);

int evaluatePopulation(const std::vector<std::vector<int>>& population, const std::vector<Vertex> verticies);

std::vector<int> getBestSpecimenFromPopulation(const std::vector<Vertex>& verticies,  std::vector<std::vector<int>>& population);

std::vector<std::vector<int>> getPopulationSample(std::vector<std::vector<int>>& population, const int sampleSize, std::mt19937& mt);

/// @brief Selects specimens from population. Each specimen is chosen as best specimen from random sample of population. 
/// @param population initial population parents are chosen from
/// @param verticies data of verticies each specimen of population is built upon
/// @param noOfPairs no. of parents selected = 2 * noOfPairs
/// @param sampleSize size of sample each specimen is chosen from
/// @param mt random numbers generator
/// @return vector of specimens
std::vector<std::vector<int>> chooseParents(std::vector<std::vector<int>>& population, const std::vector<Vertex>& verticies, const int noOfPairs, const int sampleSize, std::mt19937& mt);

/// @brief Randomly selects specimens from population
/// @param population initial population parents are choosen from
/// @param noOfPairs no. of parents selected = 2 * noOfPairs
/// @param mt random numbers generator
/// @return vector of specimens
std::vector<std::vector<int>> chooseParents(std::vector<std::vector<int>>& population, const int noOfPairs, std::mt19937& mt);

/// @brief Creates child from parents. parent2 contributes range: [startIdx + 1, stopIdx - 1] of its genes
/// @param parent1 first element from pair
/// @param parent2 second element form pair
/// @param startIdx 
/// @param stopIdx 
/// @return std::vector<int> representing child
std::vector<int> breedChild(const std::vector<int>& parent1, const std::vector<int>& parent2, const int startIdx, const int stopIdx);

/// @brief Creates child from parents.
/// @param parent1 first element from pair
/// @param parent2 second element from pair 
/// @param crossingPoint parent2 contributes range: [0, crossingPoint -1] of its genes
/// @return std::vector<int> representing child
std::vector<int> breedChild(const std::vector<int>& parent1, const std::vector<int>& parent2, const int crossingPoint);

/// @brief Performs crossover of specimens in population using parents and two crossing points.
/// @param population initial population
/// @param setOfParents vector of parents
/// @param verticies data of verticies each specimen of population is built upon
/// @param mt random number generator 
/// @param childrenPerPair no. of children each pair of parents breeds; 2 by default
void twoPointCrossover(std::vector<std::vector<int>>& population, const std::vector<std::vector<int>>& setOfParents, const std::vector<Vertex>& verticies, std::mt19937& mt, const int childrenPerPair = 2);

/// @brief Performs crossover of specimens in population using parents and one crossing point
/// @param population initial population
/// @param setOfParents vector of parents
/// @param verticies data of verticies each specimen of population is built upon
/// @param mt random number generator 
/// @param childrenPerPair no. of children each pair of parents breeds; 2 by default
void onePointCrossover(std::vector<std::vector<int>>& population, const std::vector<std::vector<int>>& setOfParents, const std::vector<Vertex>& verticies, std::mt19937& mt, const int childrenPerPair = 2);

/// @brief Mutates population by performing 2-OPT invert
/// @param population target of mutation
/// @param probOfMutation probability of mutation of specimen in population; range [0.0, 1.0].
/// @param mt random numbers generator
void mutatePopulation(std::vector<std::vector<int>>& population, const float probOfMutation, std::mt19937& mt);



/// @brief Genetic algorithm solving TSP.
/// 
/// Algorithm details:
/// - Crossover procedure: Two-point crossover with use of partially-mapped crossover (PMX).
/// - Parents selection: Tournament - best specimen from random sample of population.
/// - After each crossover the worst part of population is rejected, so that population size stayes the same
/// 
/// @param initialCycle cycle algorithm starts optimization from.
/// @param iterations number of iterations (generations).
/// @param maxIterWithoutImprov safety guard preventing from finishing optimization too early.
/// @param populationSize size of population algorithm operates on (if high, significantly increases computation time).
/// @param noOfPairs number of pairs of parents (each pair breeds children).
/// @param sampleSize size of population parents are selected from during 'parents selection process'.
/// @param mutationProb probability of mutation for each specimen in population. Should be in range [0.0, 1.0].
/// @param verticies data of vertices each specimen of population is built upon.
/// @param mt random numbers generator.
/// @return Minimal cost of cycle.
int geneticTSPSolver(
    std::vector<int>& initialCycle, const int iterations, const int maxIterWithoutImprov,
    const int populationSize, const int noOfPairs, const int sampleSize,
    const float mutationProb, const std::vector<Vertex>& verticies, std::mt19937& mt);



/// @brief Genetic algorithm solving TSP.
/// 
/// Algorithm details:
/// - Crossover procedure: One-point crossover.
/// - Parents selection: Random specimen from population.
/// - After each crossover the worst part of population is rejected, so that population size stayes the same
///
/// @param initialCycle cycle algorithm starts optimization from.
/// @param iterations number of iterations (generations).
/// @param maxIterWithoutImprov safety guard preventing from finishing optimization too early.
/// @param populationSize size of population algorithm operates on (if high, significantly increases computation time).
/// @param noOfPairs number of pairs of parents (each pair breeds children).
/// @param mutationProb probability of mutation for each specimen in  population. Should be in range [0.0, 1.0].
/// @param verticies data of vertices each specimen of population is built upon.
/// @param mt random numbers generator.
/// @return Minimal cost of cycle.
int geneticTSPSolver(
    std::vector<int>& initialCycle, const int iterations, const int maxIterWithoutImprov,
    const int populationSize, const int noOfPairs, const float mutationProb,
    const std::vector<Vertex>& verticies, std::mt19937& mt);


#endif