#include "DisjoinSets.hpp"

DisjoinSets::DisjoinSets(int V): n{V}, parent{}, rank{} 
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


int DisjoinSets::find(int u) 
{                                                                   
    if (u != parent[u]) 
    {
        parent[u] = find(parent[u]);
    }
    return parent[u];
}


void DisjoinSets::merge(int u, int v) 
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