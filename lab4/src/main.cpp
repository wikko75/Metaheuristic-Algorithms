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
    DisjoinSets(int V): n{V}, parent{}, rank{} 
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


int weightDiff(const std::vector<Vertex> &verticies, std::vector<int>& cycle, int startIdx, int stopIdx)
{

    if (startIdx > stopIdx)
    {
        int temp {startIdx};
        startIdx = stopIdx;
        stopIdx = temp;
    }
    
    //no change at all
    if (startIdx == 0 && stopIdx == verticies.size() - 1)
    {
        return 0;
    }

    //edges: [start --- (stop+1)] and [stop --- lastVertex] are changed
    if (startIdx == 0)
    { 
        //add weight between [start --- lastVertex] and [stop --- (stop + 1)]  (initial state)
        int initialCycle = euclideanNorm(verticies[cycle[startIdx] - 1], verticies[cycle[cycle.size() - 1] - 1])
                            +   euclideanNorm(verticies[cycle[stopIdx] - 1], verticies[cycle[stopIdx + 1] - 1]);

        //add weight between [start --- (stop + 1)] and [stop --- lastVertex]   (state after inversion)
        int newCycleChange = euclideanNorm(verticies[cycle[startIdx] - 1], verticies[cycle[stopIdx + 1] - 1])
                            +   euclideanNorm(verticies[cycle[stopIdx] - 1], verticies[cycle[cycle.size() - 1] - 1]);
        
        return initialCycle - newCycleChange;
    }

    //edges: [(start - 1) --- start] and [stop --- firstVertex] are changed
    if(stopIdx == verticies.size() - 1)
    {
         //add weight between [(start - 1) --- start] and [stop --- firstVertex]  (initial state)
        int initialCycle = euclideanNorm(verticies[cycle[startIdx - 1] - 1], verticies[cycle[startIdx] - 1])
                            +   euclideanNorm(verticies[cycle[stopIdx] - 1], verticies[cycle[0] - 1]);

        //add weight between [start --- firstVertex] and [stop --- (start - 1)]   (state after inversion)
        int newCycleChange = euclideanNorm(verticies[cycle[startIdx] - 1], verticies[cycle[0] - 1])
                            +   euclideanNorm(verticies[cycle[stopIdx] - 1], verticies[cycle[startIdx - 1] - 1]);
        
        return initialCycle - newCycleChange;
    }

    //add weight between [(start - 1) --- start]  and [stop --- (stop + 1)]  (initial state)
    int initialCycle = euclideanNorm(verticies[cycle[startIdx - 1] - 1], verticies[cycle[startIdx] - 1])
                        +   euclideanNorm(verticies[cycle[stopIdx] - 1], verticies[cycle[stopIdx + 1] - 1]);

    //add weight between [(start) --- (stop + 1)] and [stop --- (start - 1)]   (state after inversion)
    int newCycleChange = euclideanNorm(verticies[cycle[startIdx] - 1], verticies[cycle[stopIdx + 1] - 1])
                        +   euclideanNorm(verticies[cycle[stopIdx] - 1], verticies[cycle[startIdx - 1] - 1]);    

    return initialCycle - newCycleChange;
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
    //prob. not neccessary, compiler's gonna perform optimalization
    return std::move(population);
}


int evaluatePopulation(const std::vector<std::vector<int>>& population, const std::vector<Vertex> verticies)
{
    int minCost {std::numeric_limits<int>::max()};
    int specimenCost {};

    for (const auto& specimen : population)
    {
        specimenCost = calculateCycle(verticies, specimen);

        fmt::print("Specimen cost: {}\n", specimenCost);
        if (specimenCost < minCost)
        {
            minCost = specimenCost;
        }
    }

    fmt::print("MinCost: {}\n", minCost);
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
        //get random smaple of specimens of population 
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
    }

    return setOfParents;
}


void crossover(std::vector<std::pair<std::vector<int>, std::vector<int>>>& setOfParents, int populationSize)
{
    int noOfPairs {setOfParents.size()};
    int childrenPerPair {static_cast<int>(std::ceil(populationSize / noOfPairs))};
    
    //children of parents set
    std::vector<std::vector<int>> newGeneration {};
    newGeneration.reserve(noOfPairs * childrenPerPair);

    for (const auto [parent1, parent2] : setOfParents)
    {
        for (int i {0}; i < childrenPerPair; ++i)
        {
            //TODO crossover
        }
    }


    //TODO delete random specimens from newGeneration so that population size stayes the same
}


void mutation()
{
    //TODO
}



int main() {
    //path to files containing verticies data in format (vertex no., x coor, y coor)
    const std::filesystem::path pathToverticiesData {std::filesystem::current_path().append("test/verticiesData")};

    // random gen. used throughout program
    std::mt19937 mt { std::random_device{}()};

    {
        for (auto &verticiesDataFile : std::filesystem::directory_iterator(pathToverticiesData))
        {       
            // std::cout << verticiesDataFile.path() << "\n";
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

            std::vector<std::vector<int>> population {createPopulation(cycle, 100, mt)};

            fmt::print("Evaluating population...\nCost: {}\n\n", evaluatePopulation(population, verticies));
            auto setOfParents {chooseParents(population, verticies, 10, 40, mt)};


        }
    }

    return 0;
}