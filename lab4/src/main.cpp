#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cmath>
#include <random> // for std::mt19937
#include <optional>
#include <limits>
#include <fmt/core.h>
#include <fmt/ranges.h>



struct Timer
{
    Timer(const std::string& name): name{name}
    {   
        start = std::chrono::high_resolution_clock::now();
        std::cout << "["<< this->name  << "] " << " Timer started!\n";
    }

    ~Timer()
    {
        end = std::chrono::high_resolution_clock::now();
        std::cout << "[" << this->name << "]" <<  " Timer ended!\n";
        //changed to double
        std::cout << "Elapsed time: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;
    }

private:
    std::chrono::_V2::system_clock::time_point start;
    std::chrono::_V2::system_clock::time_point end;
    std::string name;
};


class Vertex
{
private:
    int number;
    int x;
    int y;

public:
    Vertex(int n, int xCoord, int yCoord)
        : number{n}, x{xCoord}, y{yCoord} {}

    int getNumber() const
    {
        return number;
    }

    int getX() const
    {
        return x;

    }
    int getY() const
    {
        return y;
    }
};


class Edge 
{
private:
    int src;
    int dest;
    int weight;
    
public:
    Edge(int src, int dest, int weight)
        : src{src}, dest{dest}, weight{weight} {}

public:
    int getSrc() const { return src;}
    int getDest() const { return dest;}
    int getWeight() const { return weight;}

    bool hasVertex(int v) const
    {
        if (src == v || dest == v)
        {
            return true;
        }

        return false;
    }

    int getOppositeVertex(int v) const
    {
        if(v == src)
            return dest;

        return src;
    }
};


class DisjoinSets 
{
private:
    std::vector<int> parent;
    std::vector<int> rank;
    int n;

public:
    DisjoinSets(int V) : parent{}, rank{}, n{V} 
    {
        //verticies count starts from 1 and goes up to n ex. n = 5,  {1,2,3,4,5}
        rank.reserve(n+1);
        parent.reserve(n+1);
        rank.emplace_back(-1);
        parent.emplace_back(-1);

        for (int i = 1; i < n + 1; ++i)
        {
            // Initially, all verticies are in 
            // different sets and have rank 0.
            rank.emplace_back(0);
            //every element is parent of itself 
            parent.emplace_back(i);
        }   
    }

    // Find the parent of a node 'u' 
    // Path Compression 
    int find(int u) 
    {                                                                   
        if (u != parent[u]) 
        {
            parent[u] = find(parent[u]);
        }
        return parent[u];
    }


    void merge(int u, int v) 
    {
        int uRoot = find(u);
        int vRoot = find(v);

        if (rank[uRoot] > rank[vRoot])
        {
            parent[vRoot] = uRoot;
        }
        else
        {
            parent[uRoot] = vRoot;
        }

        if (rank[uRoot] == rank[vRoot])
        {   
            ++rank[vRoot];
        }  
    }
};


class Graph
{
private:
    int V;
    std::vector<Edge> graphEdges;
    std::vector<Edge> MSTEdges;

public:
    Graph(int V) : V{V}, graphEdges{}, MSTEdges{}
    {
        graphEdges.reserve(V / 2 * (V - 1));
    }

    void addEdge(int src, int dest, int weight) 
    {
        graphEdges.emplace_back(src, dest, weight);
    }

    int kruskalMST() 
    {
        int mst_wt {0};

        std::sort(graphEdges.begin(), graphEdges.end(), [](const Edge& a, const Edge& b) {
            return a.getWeight() < b.getWeight();
        });

        std::unique_ptr<DisjoinSets> sets = std::make_unique<DisjoinSets>(V);

        for (const Edge& edge : graphEdges) 
        {
            int u = edge.getSrc();
            int v = edge.getDest();

            //find sets, where vertex u/v belongs to
            int set_u = sets->find(u);
            int set_v = sets->find(v);

            if (set_u != set_v) 
            {
                // std::cout << u << " - " << v << "\t\tweight: " << edge.getWeight() <<  std::endl;
                mst_wt += edge.getWeight();
                sets->merge(set_u, set_v);
                MSTEdges.push_back(edge);
            }
        }
        return mst_wt;
    }
    
    int getNoOfGraphEdges() const 
    {
        return this->graphEdges.size();
    }

    std::vector<Edge> getGraphEdges() const 
    {
        return graphEdges;
    }

    std::vector<Edge> getMSTEdges() const
    {
        return MSTEdges;
    }
};


void getGraphDataFromFile(const std::filesystem::path& filePath, std::vector<Vertex>& verticies, int startLine = 9)
{
    std::fstream fileStream{};
    fileStream.open(filePath, std::ios::in);

    if (!fileStream.is_open())
    {
        std::cerr << "Can't open a file\n";
        return;
    }

    int currentLineCount{0};
    std::string line{};
    
    while (!fileStream.eof())
    {
        ++currentLineCount;
        std::getline(fileStream, line);

        if (line == "EOF")
        {
            fileStream.close();
            return;
        }
        
        if (currentLineCount >= startLine)
        {
            std::stringstream strStream(line);
            int vertexNumber{};
            int x{};
            int y{};

            if (strStream >> vertexNumber >> x >> y)
            {
                verticies.emplace_back(vertexNumber, x, y);
            }
            else
            {
                std::cerr << "Can't convert line to numbers :(\n";
            }
        }
    }

    fileStream.close();   
}


int euclideanNorm(const Vertex &v1, const Vertex &v2) {
    double dx = v1.getX() - v2.getX();
    double dy = v1.getY() - v2.getY();
    double distance = std::sqrt(dx * dx + dy * dy);
    
    return static_cast<int>(std::round(distance));
}


std::unique_ptr<Graph> createCompleteGraph(std::vector<Vertex>& verticies)
{
    std::unique_ptr<Graph> graph = std::make_unique<Graph>(verticies.size());
    
    //construct complete graph with weight as Euclidian norm
    for (int i{0}; i < verticies.size() - 1; ++i)
    {
        Vertex vertex1 = verticies[i];
        for (int j{i + 1}; j < verticies.size(); ++j)
        {
            Vertex vertex2 = verticies[j];
            int weight = euclideanNorm(vertex1, vertex2);

            //add edge between vertex1 and vertex2 with proper weight 
            graph->addEdge(vertex1.getNumber(), vertex2.getNumber(), weight);
            // std::cout << "Edge formed between (" << vertex1.getNumber() << ", " << vertex2.getNumber() << ")\n";
            // std::cout <<"Weight = " <<  weight << "\n";
        }
    }

    return std::move(graph);
}


void DFS(int v, const Graph& graph, const std::vector<Edge>& MSTEdges, std::vector<bool>& visited, std::vector<int>& cycle) {
    visited[v - 1] = true;
    cycle.emplace_back(v);

    // look for edges in MSTEdges which are conneted to v
    for (const Edge& edge : MSTEdges) {
        if (edge.hasVertex(v))
        {
            int nextVertex = edge.getOppositeVertex(v);
            if (!visited[nextVertex - 1])
            {
                DFS(nextVertex, graph, MSTEdges, visited, cycle);
            }
        }
    }
}


int calculateCycle(const std::vector<Vertex>& verticies, const std::vector<int>& cycle)
{
    int cycleWeight {};

    //add all weights between first and last element in cycle
    for(int index{0}; index < cycle.size() - 1;  ++index)
    {
        cycleWeight += euclideanNorm(verticies[cycle[index] - 1], verticies[cycle[index + 1] - 1]);
    }

    //add weight between last and first element in cycle (cause we have to come back, since it's a cycle)
    cycleWeight += euclideanNorm(verticies[cycle[0] - 1], verticies[cycle[cycle.size() - 1] - 1]);

    return cycleWeight;
}


/// @brief Saves cycle edges to provided file
/// @param path path to file, if file does not exist it will be created
/// @param cycle vector storing verticies' numbers that form cycle
void saveCycle(const std::string& path, const std::vector<int>& cycle)
{
    std::ofstream fstream {};

    fstream.open(path, std::ios::out);

    std::cout << "Path to file: " <<  path << "\n";

    if(!fstream.is_open())
    {
        std::cerr << "Can't write cycle edges\nAborting...\n";
        fstream.close();
    }
    else
    {
        std::cout << "Saving edges...\n";
        for(int i {0}; i < cycle.size() - 1; ++i)
        {     
            fstream << cycle[i] << " " << cycle[i + 1] << "\n";   
        }

        fstream << cycle[cycle.size() - 1] << " " << cycle[0] << "\n";
    }

    fstream.close();
}


struct ExperimentData
{
    const char* label;
    int bestCycleWeight{};
    int avgCycleWeight{};
};
        

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


void invert(std::vector<int>& cycle, int i, int j)
{
    i <= j ? std::reverse(cycle.begin() + i, cycle.begin() + j + 1)
           : std::reverse(cycle.begin() + j, cycle.begin() + i + 1);
}


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

        // fmt::print("Specimen cost: {}\n", specimenCost);
        if (specimenCost < minCost)
        {
            minCost = specimenCost;
        }
    }

    // fmt::print("MinCost: {}\n", minCost);
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

std::vector<std::pair<std::vector<int>, std::vector<int>>> chooseParents(
    std::vector<std::vector<int>>& population,
    const std::vector<Vertex>& verticies,
    const int noOfPairs,
    const int sampleSize,
    std::mt19937& mt)
{
    if (sampleSize <= 0 || sampleSize > population.size())
    {
        fmt::print("Can't choose parents. Wrong sampleSize [{}] specified! Population size: {}\n", sampleSize, population.size());
    }

    if (noOfPairs <= 0 || noOfPairs > population.size() / 2)
    {
        fmt::print("Can't choose parents. Wrong noOfPairs [{}] specified! Population size: {}\n", noOfPairs, population.size());
    }
    
    std::vector<std::pair<std::vector<int>, std::vector<int>>> setOfParents {};
    setOfParents.reserve(noOfPairs);

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

        setOfParents.emplace_back(std::make_pair(firstParent, secondParent));
          
        // fmt::print("Parents:\n{}\n\n{}\n------\n", firstParent, secondParent); 
    }

    return setOfParents;
}


/// @brief Randomly
/// @param population 
/// @param verticies 
/// @param noOfPairs 
/// @param mt 
/// @return 
std::vector<std::vector<int>> chooseParents(
    std::vector<std::vector<int>>& population,
    const int noOfPairs,
    std::mt19937& mt)
{
    if (noOfPairs <= 0 || noOfPairs > population.size() / 2)
    {
        fmt::print("Can't choose parents. Wrong noOfPairs [{}] specified! Population size: {}\n", noOfPairs, population.size());
    }
    
    std::vector<std::vector<int>> setOfParents {};
    setOfParents.reserve(noOfPairs);

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

        // fmt::print("Parents:\n{}\n\n{}\n------\n", population[j], population[k]); 
    }

    return setOfParents;
}


std::vector<int> breedChild(const std::vector<int>& parent1, const std::vector<int>& parent2, const int startIdx, const int stopIdx)
{
    int cycleLen {static_cast<int>(parent1.size())};
    std::vector<int> sequenceMap {};
    sequenceMap.resize(cycleLen);

    std::vector<int> child {};
    child.reserve(cycleLen);


    //copy parent1 cycle to child
    std::copy(parent1.begin(), parent1.end(), std::back_inserter(child));

    for (int i {0}; i < cycleLen; ++i)
    {
        //example: (Vertex 4) - 1 is on index 5 in child cycle    // -1 because we wanna idx of vertex 4. (again verticies are 1...n; indicies are 0,..,n-1)
        //with permutation 1..n it simply states that (vertex 1) - 1 is on index 0, (vertex 2) - 1 is on index 1 and so on...
        sequenceMap[child[i] - 1] = i;  
    }

    //copy verticies from parent2 into child (taking map info into consideration)
    for(int i {startIdx + 1}; i < stopIdx; ++i)
    {
        int vertex {parent2[i]};
        int vertexIdxInMap = sequenceMap[vertex - 1];
        
        //swap verticies 
        child[vertexIdxInMap] = child[i];
        child[i] = parent1[vertexIdxInMap];

        //swap indicies
        int temp = sequenceMap[child[vertexIdxInMap] - 1];
        sequenceMap[child[vertexIdxInMap] - 1] = vertexIdxInMap;
        sequenceMap[vertex - 1] = temp; 
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
    
    // fmt::print("Child created:\n {}", child);
    return child;
}


std::vector<std::vector<int>> crossover(std::vector<std::pair<std::vector<int>, std::vector<int>>>& setOfParents, const int populationSize, std::mt19937& mt)
{
    int noOfPairs {static_cast<int>(setOfParents.size())};
    int childrenPerPair {static_cast<int>(std::ceil(populationSize / (float)noOfPairs))};
    int cycleLen {static_cast<int>(setOfParents[0].first.size())};
    
    std::vector<std::vector<int>> nextGeneration {};
    nextGeneration.reserve(noOfPairs * childrenPerPair);

    //verticies idx  count: 1...cycleLen
    std::uniform_int_distribution randomVertexIdx {0, cycleLen - 1};
    std::uniform_real_distribution randomDouble {0.0, 1.0};

    //(startIdx to stopIdx) is taken from second parent
    //[0 to startIdx] and [stopIdx to end] is taken from first parent
    int startIdx {};
    int stopIdx {};

    for (const auto& [parent1, parent2] : setOfParents)
    { 
        std::vector<int> sequenceMap {};
        sequenceMap.reserve(cycleLen + 1);

        for (int i {0}; i < childrenPerPair; ++i)
        {
            startIdx = randomVertexIdx(mt);
            do
            {
                stopIdx = randomVertexIdx(mt);
            } while (
                abs(startIdx - stopIdx) <= 2              //abs(startIdx - stopIdx) < 2 because we want at least one vertex taken from second parent
            );

            if (startIdx > stopIdx)
            {
                int temp {startIdx};
                startIdx = stopIdx;
                stopIdx = temp;
            }
            
            double transmissionProbability {randomDouble(mt)};
            // fmt::print("startIdx: {}, stopIdx: {}\n", startIdx, stopIdx);
            if (transmissionProbability < 0.5)
            {
                std::vector<int> firstChild {breedChild(parent1, parent2, startIdx, stopIdx)};
                nextGeneration.emplace_back(firstChild);
                // fmt::print("Child:\n{}\n", firstChild);
            }
            else
            {
                std::vector<int> secondChild {breedChild(parent2, parent1, startIdx, stopIdx)};
                nextGeneration.emplace_back(secondChild);
                // fmt::print("Child:\n{}\n", secondChild);
            }
        }
    }

    if (nextGeneration.size() == populationSize)
    {
        return nextGeneration;
    }
    
    //delete random specimens from nextGeneration so that population size stayes the same
    const int offset {nextGeneration.size() - populationSize};
    std::shuffle(nextGeneration.begin(), nextGeneration.end(), mt);

    for (int i {0}; i < offset; ++i)
    {
        nextGeneration.erase(nextGeneration.begin() + i);
    }
    
    nextGeneration.shrink_to_fit();

    return nextGeneration;
}


void crossoverV2(std::vector<std::vector<int>>& population, const std::vector<std::vector<int>>& setOfParents, std::mt19937& mt, const std::vector<Vertex>& verticies)
{
    int noOfPairs {static_cast<int>(setOfParents.size() / 2)};
    int cycleLen {static_cast<int>(setOfParents[0].size())};
    const int initialPopulationSize {population.size()};
    
    //verticies idx  count: 1...cycleLen
    std::uniform_int_distribution randomVertexIdx {1, cycleLen - 2};
    std::uniform_real_distribution randomDouble {0.0, 1.0};

    //[0 to crossingIdx] is taken from first second parent
    int crossingIdx {};

    for (int i {0}; i < setOfParents.size() - 1; ++i)
    { 

        for (int j {1}; j < setOfParents.size(); ++j)
        {
            crossingIdx = randomVertexIdx(mt);
            
            // fmt::print("crossingIdx: {}\n", crossingIdx);
            
            std::vector<int> firstChild {breedChild(setOfParents[i], setOfParents[j], crossingIdx)};
            population.emplace_back(firstChild);
            // fmt::print("Child:\n{}\n", firstChild);
            
            std::vector<int> secondChild {breedChild(setOfParents[j], setOfParents[i], crossingIdx)};
            population.emplace_back(secondChild);
            // fmt::print("Child:\n{}\n", secondChild);
        }
    }


    
    //delete worse specimens from nextGeneration so that population size stayes the same
    const int offset {population.size() - initialPopulationSize};
    std::shuffle(population.begin(), population.end(), mt);
    std::sort(population.begin(), population.end(), [&verticies](const auto& specimen1, const auto specimen2)
    {
        return calculateCycle(verticies, specimen1) < calculateCycle(verticies, specimen2);
    });

    for (int i {0}; i < offset; ++i)
    {
        population.pop_back();
    }
    // fmt::print("Population size after crossover: {}\n", population.size());
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

    std::vector<std::vector<int>> population {createPopulation(initialCycle, populationSize, mt)};

    int minCost {std::numeric_limits<int>::max()};
    int currIter {0};
    int currIterWithoutImprov {0};
    int currMinCost {evaluatePopulation(population, verticies)};
    // fmt::print("Evaluating initial population...\nCost: {}\n\n", currMinCost);

    while (currIter < iterations && currIterWithoutImprov < maxIterWithoutImprov)
    {
        minCost = currMinCost;
        // auto setOfParents {chooseParents(population, verticies, noOfPairs, sampleSize, mt)};
        std::vector<std::vector<int>> setOfParents {chooseParents(population, noOfPairs, mt)};
       
        // std::vector<std::vector<int>> nextGeneration {crossover(setOfParents, population.size(), mt)};
        // std::vector<std::vector<int>> nextGeneration {crossoverV2(setOfParents, population.size(), mt)};
        crossoverV2(population, setOfParents, mt, verticies);
        //fmt::print("Performing crossover...\n");
        currMinCost = evaluatePopulation(population, verticies);
        // fmt::print("Evaluating population after crossover...\nCost: {}\n", currMinCost);
        
        //fmt::print("Performing mutation...\n");
        mutatePopulation(population, mutationProb, mt);
        currMinCost = evaluatePopulation(population, verticies);
        // fmt::print("Evaluating population after mutation...\nCost: {}\n\n", currMinCost);
        
        //fmt::print("New population size: {}\n", nextGeneration.size());

        if (currMinCost == minCost)
        {
            ++currIterWithoutImprov;
        }
        else
        {
            currIterWithoutImprov = 0;
        }
        
        // population = nextGeneration;
        ++currIter;
    }

    fmt::print("Performed: {} iterations\n\n", currIter);

    return minCost;
}

int main() {
    //path to files containing verticies data in format (vertex no., x coor, y coor)
    const std::filesystem::path pathToverticiesData {std::filesystem::current_path().append("test/verticiesData")};

    //random gen. used throughout program
    std::mt19937 mt {std::random_device{}()};

    {
        for (auto &verticiesDataFile : std::filesystem::directory_iterator(pathToverticiesData))
        {       
            fmt::print("Path: {}\n", verticiesDataFile.path().string());
            
            std::vector<Vertex> verticies {};

            //populate verticies with data
            getGraphDataFromFile(verticiesDataFile.path(), verticies);
            
            std::vector<int> cycle{};
            cycle.reserve(verticies.size());

            for (int i {1}; i < verticies.size() + 1; ++i)
            {
                cycle.push_back(i);
            }

            int minCost {std::numeric_limits<int>::max()};
            int avgCost {};
            int repeats {10};
            for (int i {0}; i < repeats; ++i)
            {
                const int iter {300};
                const int populationSize {30};
                const int noOfPairs {10};
                const int populationSample {20};
                const float mutationProb {.3};
                const int maxIterWithoutImprov {30};
                // const int iter {100};
                // const int populationSize {50};
                // const int noOfPairs {10};
                // const int populationSample {20};
                // const float mutationProb {.20};
                int currMinCost {geneticTSPSolver(cycle, iter, maxIterWithoutImprov, populationSize, noOfPairs, populationSample, mutationProb, verticies, mt)};
                // fmt::print("Cost after {} iterations: {}\n", iter, currMinCost);

                if (currMinCost < minCost)
                {
                    minCost = currMinCost;
                }

                avgCost += currMinCost;
            }

            fmt::print("\nResults after {} repeats:\nMinimal cost: {}\nAvg. cost {}\n", repeats, minCost, avgCost/repeats);

            const std::filesystem::path path {std::filesystem::current_path().append("res/resultsV2")
                                            .append("GenTSP_DATA" + verticiesDataFile.path().filename().string())};

            fmt::print("Save path: {}\n", path.string());
            const ExperimentData data {"GenTSP", minCost, avgCost/repeats};
            saveExperimentData(path, data);
        }
    }

    return 0;
}



/*//TODO problemen jest chyba parenSet !!! trzeba zmienic koncepcje. Generujemy 

*/