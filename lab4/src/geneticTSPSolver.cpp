#include "geneticTSPSolver.hpp"
#include <algorithm>
#include <fmt/core.h>
#include "utilsFunc.hpp"


std::vector<std::vector<int>> createPopulation(std::vector<int>& cycle, int populationSize, std::mt19937& mt)
{
    std::vector<std::vector<int>> population {};
    population.reserve(populationSize);

    for (int i {0}; i < populationSize; ++i)
    {
        std::shuffle(cycle.begin(), cycle.end(), mt);
        population.emplace_back(cycle);
    }
    return population;
}


int evaluatePopulation(const std::vector<std::vector<int>>& population, const std::vector<Vertex> verticies)
{
    int minCost {std::numeric_limits<int>::max()};
    int specimenCost {};

    for (const auto& specimen : population)
    {
        specimenCost = calculateCycle(verticies, specimen);

        if (specimenCost < minCost)
        {
            minCost = specimenCost;
        }
    }

    return minCost;
}


std::vector<int> getBestSpecimenFromPopulation(const std::vector<Vertex>& verticies,  std::vector<std::vector<int>>& population)
{
    std::sort(population.begin(), population.end(), [&verticies](const auto& specimen1, const auto& specimen2)
    {
        return calculateCycle(verticies, specimen1) < calculateCycle(verticies, specimen2);
    });

    return population[0];
}


std::vector<std::vector<int>> getPopulationSample(std::vector<std::vector<int>>& population, const int sampleSize, std::mt19937& mt)
{    
    std::shuffle(population.begin(), population.end(), mt);

    return {population.begin(), population.begin() + sampleSize};
}


std::vector<std::vector<int>> chooseParents(
    std::vector<std::vector<int>>& population,
    const std::vector<Vertex>& verticies,
    const int noOfPairs,
    const int sampleSize,
    std::mt19937& mt)
{
    if (sampleSize <= 0 || sampleSize > population.size())
    {
        fmt::print("Can't choose parents. Wrong sampleSize [{}] specified! Population size: {}\n", sampleSize, population.size());
        return {};
    }

    if (noOfPairs <= 0 || noOfPairs > population.size() / 2)
    {
        fmt::print("Can't choose parents. Wrong noOfPairs [{}] specified! Population size: {}\n", noOfPairs, population.size());
        return {};
    }
    
    std::vector<std::vector<int>> setOfParents {};
    setOfParents.reserve(noOfPairs * 2);

    for (int i{0}; i < noOfPairs; ++i)
    {
        //get random sample of specimens of population 
        std::vector<std::vector<int>> populationSample {getPopulationSample(population, sampleSize, mt)};
        
        //get best specimen from population sample - it's our parent
        std::vector<int> firstParent {getBestSpecimenFromPopulation(verticies, populationSample)};
        
        
        std::vector<int> secondParent {};
        secondParent.reserve(verticies.size());

        //repeat for second parent
        //ensure both parenst are unique
        do
        {
            std::vector<std::vector<int>> populationSample {getPopulationSample(population, sampleSize, mt)};
            secondParent = getBestSpecimenFromPopulation(verticies, populationSample);
        } while (std::equal(firstParent.begin(), firstParent.end(), secondParent.begin(), secondParent.end()));

        setOfParents.emplace_back(firstParent);
        setOfParents.emplace_back(secondParent);
    }

    return setOfParents;
}


std::vector<std::vector<int>> chooseParents(
    std::vector<std::vector<int>>& population,
    const int noOfPairs,
    std::mt19937& mt)
{
    if (noOfPairs <= 0 || noOfPairs > population.size() / 2)
    {
        fmt::print("Can't choose parents. Wrong noOfPairs [{}] specified! Population size: {}\n", noOfPairs, population.size());
        return {};
    }
    
    std::vector<std::vector<int>> setOfParents {};
    setOfParents.reserve(noOfPairs * 2);

    std::uniform_int_distribution randomIdx {0, static_cast<int>(population.size()) - 1};
    int j {};
    int k {};
    for (int i{0}; i < noOfPairs; ++i)
    {

        j = randomIdx(mt);

        do
        {
            k = randomIdx(mt);
        } while (j == k);

        setOfParents.emplace_back(population[j]);
        setOfParents.emplace_back(population[k]);

    }

    return setOfParents;
}


std::vector<int> breedChild(const std::vector<int>& parent1, const std::vector<int>& parent2, const int startIdx, const int stopIdx)
{
    const int cycleLen {static_cast<int>(parent1.size())};

    std::vector<int> child {};
    child.reserve(cycleLen);

    std::copy(parent1.begin(), parent1.end(), std::back_inserter(child));

    //copy verticies from parent2 into child (startIdx, stopIdx)
    for(int i {startIdx + 1}; i < stopIdx; ++i)
    {
        auto it {std::find(child.begin(), child.end(), parent2[i])};

        child[it - child.begin()] = child[i];
        child[i] = parent2[i];
    }
    
    return child;
}


std::vector<int> breedChild(const std::vector<int>& parent1, const std::vector<int>& parent2, const int crossingPoint)
{
    int cycleLen {static_cast<int>(parent1.size())};
    std::vector<int> child {};
    child.reserve(cycleLen);

    //copy parent1 cycle to child
    std::copy(parent1.begin(), parent1.end(), std::back_inserter(child));

    for (int i {0}; i < crossingPoint; ++i)
    {
        auto it {std::find(child.begin(), child.end(), parent2[i])};

        child[it - child.begin()] = child[i];
        child[i] = parent2[i];
    }
    
    return child;
}


void twoPointCrossover(std::vector<std::vector<int>>& population, 
               const std::vector<std::vector<int>>& setOfParents,
               const std::vector<Vertex>& verticies, 
               std::mt19937& mt,
               const int childrenPerPair)
{
    int cycleLen {static_cast<int>(population[0].size())};
    const int initialPopulationSize {population.size()};
    
    //verticies idx  count: 1...cycleLen
    std::uniform_int_distribution randomVertexIdx {0, cycleLen - 1};
    std::uniform_real_distribution randomDouble {0.0, 1.0};

    //(startIdx to stopIdx) is taken from second parent
    //[0 to startIdx] and [stopIdx to end] is taken from first parent
    int startIdx {};
    int stopIdx {};
    int j {};

    for (int i {0}; i < setOfParents.size() - 1; i += 2)
    { 
        j = i + 1; //index of specimen right next to specimen with idx i

        for (int k {0}; k < childrenPerPair; ++k)
        {
            startIdx = randomVertexIdx(mt);
            do
            {
                stopIdx = randomVertexIdx(mt);
            } while (
                abs(startIdx - stopIdx) < 2     //abs(startIdx - stopIdx) < 2 because we want at least one vertex taken from second parent
                );

            if (startIdx > stopIdx)
            {
                int temp {startIdx};
                startIdx = stopIdx;
                stopIdx = temp;
            }
            
            double transmissionProbability {randomDouble(mt)};

            if (transmissionProbability <= 0.5)
            {
                std::vector<int> child {breedChild(setOfParents[i], setOfParents[j], startIdx, stopIdx)};
                population.emplace_back(child);
            }
            else
            {
                std::vector<int> child {breedChild(setOfParents[j], setOfParents[i], startIdx, stopIdx)};
                population.emplace_back(child);
            }
        }
    }

    if (population.size() == initialPopulationSize)
    {
        return;
    }
    
    //delete worse specimens from nextGeneration so that population size stayes the same
    const int offset {population.size() - initialPopulationSize};

    std::sort(population.begin(), population.end(), [&verticies](const auto& specimen1, const auto specimen2)
    {
        return calculateCycle(verticies, specimen1) < calculateCycle(verticies, specimen2);
    });

    for (int i {0}; i < offset; ++i)
    {
        population.pop_back();
    }

    population.shrink_to_fit();

}


void onePointCrossover(std::vector<std::vector<int>>& population,
                const std::vector<std::vector<int>>& setOfParents,
                const std::vector<Vertex>& verticies,
                std::mt19937& mt,
                const int childrenPerPair)
{
    int cycleLen {static_cast<int>(population[0].size())};
    const int initialPopulationSize {population.size()};
    
    std::uniform_int_distribution randomVertexIdx {1, cycleLen - 2};
    std::uniform_real_distribution randomDouble {0.0, 1.0};

    //[0 to crossingIdx] is taken from second parent
    int crossingIdx {};
    int j {};

    for (int i {0}; i < setOfParents.size() - 1; i += 2)
    { 
        j = i + 1;
        for (int k {0}; k < childrenPerPair; ++k)
        {
            crossingIdx = randomVertexIdx(mt);
            
            double transmissionProbability {randomDouble(mt)};
            if (transmissionProbability <= 0.5)
            {
                std::vector<int> child {breedChild(setOfParents[i], setOfParents[j], crossingIdx)};
                population.emplace_back(child);
            }
            else
            {
                std::vector<int> child {breedChild(setOfParents[j], setOfParents[i], crossingIdx)};
                population.emplace_back(child);
            }
        }
    }


    if (population.size() == initialPopulationSize)
    {
        return;
    }

    //delete worse specimens from population so that its size stayes the same
    const int offset {population.size() - initialPopulationSize};

    std::sort(population.begin(), population.end(), [&verticies](const auto& specimen1, const auto specimen2)
    {
        return calculateCycle(verticies, specimen1) < calculateCycle(verticies, specimen2);
    });

    for (int i {0}; i < offset; ++i)
    {
        population.pop_back();
    }

    population.shrink_to_fit();
}


void mutatePopulation(std::vector<std::vector<int>>& population, const float probOfMutation, std::mt19937& mt)
{
    const int cycleLen {static_cast<int>(population[0].size())};
    std::uniform_real_distribution randomDouble {0.0, 1.0};
    std::uniform_int_distribution randomIdx {0, cycleLen - 1};
    float num {};
    int i {};
    int j {};

    for (auto& specimen : population)
    {
        num = static_cast<float>(randomDouble(mt));
        if (num > probOfMutation)
        {
            continue;
        }

        i = randomIdx(mt);

        do
        {
            j = randomIdx(mt);
        } while (i == j);

        invert(specimen, i, j);
    }
}


int geneticTSPSolver(
    std::vector<int>& initialCycle, const int iterations, const int maxIterWithoutImprov,
    const int populationSize, const int noOfPairs, const int sampleSize,
    const float mutationProb, const std::vector<Vertex>& verticies, std::mt19937& mt)
{

    if (noOfPairs <= 0 || noOfPairs > populationSize / 2)
    {
        fmt::print("Wrong noOfPairs [{}] specified! Population size: {}\n", noOfPairs, populationSize);
        return -1;
    }

    if (sampleSize <= 0 || sampleSize > populationSize)
    {
        fmt::print("Wrong sampleSize [{}] specified! Population size: {}\n", sampleSize, populationSize);
        return -1;
    }

    std::vector<std::vector<int>> population {createPopulation(initialCycle, populationSize, mt)};
    std::vector<std::vector<int>> setOfParents {};
    setOfParents.reserve(noOfPairs * 2);

    int minCost {std::numeric_limits<int>::max()};
    int currIter {0};
    int currIterWithoutImprov {0};
    int currMinCost {evaluatePopulation(population, verticies)};
    fmt::print("Evaluating initial population...\nCost: {}\n\n", currMinCost);


    while (currIter < iterations && currIterWithoutImprov < maxIterWithoutImprov)
    {
        minCost = currMinCost;

        setOfParents = chooseParents(population, verticies, noOfPairs, sampleSize, mt);

        twoPointCrossover(population, setOfParents, verticies, mt);

        mutatePopulation(population, mutationProb, mt);

        currMinCost = evaluatePopulation(population, verticies);
        fmt::print("Evaluating population after mutation...\nCost: {}\n\n", currMinCost);
        
        if (currMinCost == minCost)
        {
            ++currIterWithoutImprov;
        }
        else
        {
            currIterWithoutImprov = 0;
        }
        
        ++currIter;
    }

    fmt::print("Performed: {} iterations\n\n", currIter);

    return minCost;
}


int geneticTSPSolver(
    std::vector<int>& initialCycle, const int iterations, const int maxIterWithoutImprov,
    const int populationSize, const int noOfPairs, const float mutationProb,
    const std::vector<Vertex>& verticies, std::mt19937& mt)
{

    if (noOfPairs <= 0 || noOfPairs > populationSize / 2)
    {
        fmt::print("Wrong noOfPairs [{}] specified! Population size: {}\n", noOfPairs, populationSize);
        return -1;
    }
    

    std::vector<std::vector<int>> population {createPopulation(initialCycle, populationSize, mt)};
    std::vector<std::vector<int>> setOfParents {};
    setOfParents.reserve(noOfPairs * 2);

    int minCost {std::numeric_limits<int>::max()};
    int currIter {0};
    int currIterWithoutImprov {0};
    int currMinCost {evaluatePopulation(population, verticies)};
    fmt::print("Evaluating initial population...\nCost: {}\n\n", currMinCost);


    while (currIter < iterations && currIterWithoutImprov < maxIterWithoutImprov)
    {
        minCost = currMinCost;

        setOfParents = chooseParents(population, noOfPairs, mt);

        onePointCrossover(population, setOfParents, verticies, mt);

        mutatePopulation(population, mutationProb, mt);

        currMinCost = evaluatePopulation(population, verticies);
        fmt::print("Evaluating population after mutation...\nCost: {}\n\n", currMinCost);
        
        if (currMinCost == minCost)
        {
            ++currIterWithoutImprov;
        }
        else
        {
            currIterWithoutImprov = 0;
        }
        
        ++currIter;
    }

    fmt::print("Performed: {} iterations\n\n", currIter);

    return minCost;
}