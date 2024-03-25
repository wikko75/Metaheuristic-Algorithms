#include "Edge.hpp"

Edge::Edge(int src, int dest, int weight)
    : src {src}, dest {dest}, weight {weight} {}

int Edge::getSrc() const noexcept
{
    return src;
}

int Edge::getDest() const noexcept
{ 
    return dest;
}

int Edge::getWeight() const noexcept
{ 
    return weight;
}

bool Edge::hasVertex(int v) const noexcept
{
    if (src == v || dest == v)
    {
        return true;
    }

    return false;
}

int Edge::getOppositeVertex(int v) const noexcept
{
    if(v == src)
        return dest;

    return src;
}