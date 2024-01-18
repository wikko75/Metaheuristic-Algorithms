#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <fstream>
#include<sstream>
#include <filesystem>
#include <cmath>
#include <random> // for std::mt19937
#include <chrono> // for std::chrono
#include <optional>



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
    {
        this->src = src;
        this->dest = dest;
        this->weight = weight;
    }

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
            // Initially, all vertices are in 
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

private:
    int V;
    std::vector<Edge> graphEdges;
    std::vector<Edge> MSTEdges;
};


void getGraphDataFromFile(const std::filesystem::path& filePath, std::vector<Vertex>& vertices, int startLine = 9)
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
                vertices.emplace_back(vertexNumber, x, y);
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


std::unique_ptr<Graph> createCompleteGraph(std::vector<Vertex>& vertices)
{
    std::unique_ptr<Graph> graph = std::make_unique<Graph>(vertices.size());
    
    //construct complete graph with weight as Euclidian norm
    for (int i{0}; i < vertices.size() - 1; ++i)
    {
        Vertex vertex1 = vertices[i];
        for (int j{i + 1}; j < vertices.size(); ++j)
        {
            Vertex vertex2 = vertices[j];
            int weight = euclideanNorm(vertex1, vertex2);

            //add edge between vertex1 and vertex2 with proper weight 
            graph->addEdge(vertex1.getNumber(), vertex2.getNumber(), weight);
            // std::cout << "Edge formed between (" << vertex1.getNumber() << ", " << vertex2.getNumber() << ")\n";
            // std::cout <<"Weight = " <<  weight << "\n";
        }
    }
    //std::move
    return graph;
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


int calculateCycle(const std::vector<Vertex>& vertices, const std::vector<int>& cycle)
{
    int cycleWeight {};

    //add all weights between first and last element in cycle
    for(int index{0}; index < cycle.size() - 1;  ++index)
    {
        // std::cout << "Calculating norm between: [ " << cycle[index] << ", " << cycle[index + 1] << " ]\n";
        cycleWeight += euclideanNorm(vertices[cycle[index] - 1], vertices[cycle[index + 1] - 1]);
    }

    //add weight between last and first element in cycle (cause we have to come back, since it's a cycle)
    // std::cout << "Calculating norm between: [ " << cycle[cycle.size() - 1] << ", " << cycle[0] << " ]\n";

    cycleWeight += euclideanNorm(vertices[cycle[0] - 1], vertices[cycle[cycle.size() - 1] - 1]);

    return cycleWeight;
}


int calculatePermutationCycle(const std::vector<Vertex>& permutation)
{
    int cycleWeight {};

    //add all weights between first and last element in cycle
    for(int index{0}; index < permutation.size() - 1;  ++index)
    {
        // std::cout << "Calculating norm between: [ " << permutation[index].getNumber() << ", " << permutation[index + 1].getNumber() << " ]\n";
        cycleWeight += euclideanNorm(permutation[index], permutation[index + 1]);
    }

    //add weight between last and first element in cycle (cause we have to come back, since it's a cycle)
    // std::cout << "Calculating norm between: [ " << permutation[permutation.size() - 1].getNumber() << ", " << permutation[0].getNumber() << " ]\n";

    cycleWeight += euclideanNorm(permutation[permutation.size() - 1], permutation[0]);

    return cycleWeight;
}


//I don't like this function//
/// @brief Permutates provided vertices using random generators
/// @param vertices vector of Vertex objects based on which permutation will be done
/// @param mt random generator
/// @param randomVertex uniform number generator that generates random vertices from graph  
/// @return vector of Vertex objects after permutation
std::vector<Vertex> permutateVertices(std::vector<Vertex>& vertices, std::mt19937& mt, std::uniform_int_distribution<int>& randomVertex) {
    std::vector<Vertex> permutation = vertices;  // Creating a copy to generate permutations

    // Use a different seed for better randomness
    std::mt19937 gen(mt());
    // Rearrange the elements in the vector to create permutations
    for (int i = static_cast<int>(permutation.size()) - 1; i >= 1; --i) {
        int j = randomVertex(gen, std::uniform_int_distribution<int>::param_type(0, i));  // Generating a random index directly from the distribution

        // Swap elements to create a permutation
        std::swap(permutation[i], permutation[j]);
    }

    return permutation;
}



/// @brief Saves cycle edges to provided file
/// @param path path to file, if file does not exists it will be created
/// @param cycle vector storing vertices' numbers that form cycle
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



int weightDiff(const std::vector<Vertex> &vertices, std::vector<int>& cycle, int startIdx, int stopIdx)
{

    if (startIdx > stopIdx)
    {
        int temp {startIdx};
        startIdx = stopIdx;
        stopIdx = temp;
    }
    
    //no change at all
    if (startIdx == 0 && stopIdx == vertices.size() - 1)
    {
        return 0;
    }

    //edges: [start --- (stop+1)] and [stop --- lastVertex] are changed
    if (startIdx == 0)
    { 
        //add weight between [start --- lastVertex] and [stop --- (stop + 1)]  (initial state)
        int initialCycle = euclideanNorm(vertices[cycle[startIdx] - 1], vertices[cycle[cycle.size() - 1] - 1])
                            +   euclideanNorm(vertices[cycle[stopIdx] - 1], vertices[cycle[stopIdx + 1] - 1]);

        //add weight between [start --- (stop + 1)] and [stop --- lastVertex]   (state after inversion)
        int newCycleChange = euclideanNorm(vertices[cycle[startIdx] - 1], vertices[cycle[stopIdx + 1] - 1])
                            +   euclideanNorm(vertices[cycle[stopIdx] - 1], vertices[cycle[cycle.size() - 1] - 1]);
        
        return initialCycle - newCycleChange;
    }

    //edges: [(start - 1) --- start] and [stop --- firstVertex] are changed
    if(stopIdx == vertices.size() - 1)
    {
         //add weight between [(start - 1) --- start] and [stop --- firstVertex]  (initial state)
        int initialCycle = euclideanNorm(vertices[cycle[startIdx - 1] - 1], vertices[cycle[startIdx] - 1])
                            +   euclideanNorm(vertices[cycle[stopIdx] - 1], vertices[cycle[0] - 1]);

        //add weight between [start --- firstVertex] and [stop --- (start - 1)]   (state after inversion)
        int newCycleChange = euclideanNorm(vertices[cycle[startIdx] - 1], vertices[cycle[0] - 1])
                            +   euclideanNorm(vertices[cycle[stopIdx] - 1], vertices[cycle[startIdx - 1] - 1]);
        
        return initialCycle - newCycleChange;
    }

    //add weight between [(start - 1) --- start]  and [stop --- (stop + 1)]  (initial state)
    int initialCycle = euclideanNorm(vertices[cycle[startIdx - 1] - 1], vertices[cycle[startIdx] - 1])
                        +   euclideanNorm(vertices[cycle[stopIdx] - 1], vertices[cycle[stopIdx + 1] - 1]);

    //add weight between [(start) --- (stop + 1)] and [stop --- (start - 1)]   (state after inversion)
    int newCycleChange = euclideanNorm(vertices[cycle[startIdx] - 1], vertices[cycle[stopIdx + 1] - 1])
                        +   euclideanNorm(vertices[cycle[stopIdx] - 1], vertices[cycle[startIdx - 1] - 1]);    

    return initialCycle - newCycleChange;
}




int main() {
    //path to files containing vertices data in format (no. of vertex, x coor, y coor)
    const std::filesystem::path pathToVerticesData {std::filesystem::current_path().append("res/verticesData")};

    // Seed our Mersenne Twister using steady_clock
    std::mt19937 mt { static_cast<std::mt19937::result_type>(
        std::chrono::steady_clock::now().time_since_epoch().count()
        )};

    {
        for (auto &verticesDataFile : std::filesystem::directory_iterator(pathToVerticesData))
        {       
            std::cout << verticesDataFile.path() << "\n";
            std::vector<Vertex> vertices {};

            //populate vertices with data
            getGraphDataFromFile(verticesDataFile.path(), vertices);
            
            std::vector<int> cycle{};
            cycle.reserve(vertices.size());

            for (int i {1}; i < vertices.size() + 1; ++i)
            {
                cycle.push_back(i);
            }
        }
    }

    return 0;
}