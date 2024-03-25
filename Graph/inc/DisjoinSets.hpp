#ifndef DISJOIN_SETS_1_HPP
#define DISJOIN_SETS_1_HPP

#include <vector>

class DisjoinSets 
{
private:
    std::vector<int> parent;
    std::vector<int> rank;
    int n;

public:
    DisjoinSets(int V);

    // Find the parent of a node 'u' 
    // Path Compression 
    int find(int u);

    void merge(int u, int v);

    ~DisjoinSets() = default;
};

#endif