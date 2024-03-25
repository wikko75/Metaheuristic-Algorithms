#ifndef EDGE_1_HPP
#define EDGE_1_HPP

class Edge 
{
private:
    int src;
    int dest;
    int weight;
    
public:
    Edge(int src, int dest, int weight);
    
    int getSrc() const noexcept;
    
    int getDest() const noexcept;
    
    int getWeight() const noexcept;

    bool hasVertex(int v) const noexcept;

    int getOppositeVertex(int v) const noexcept;

    ~Edge() = default;
};

#endif