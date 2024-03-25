#include "Vertex.hpp"

Vertex::Vertex(int n, int xCoord, int yCoord)
    : number{n}, x{xCoord}, y{yCoord} {}


int Vertex::getNumber() const noexcept
{
    return number;
}

int Vertex::getX() const noexcept
{
    return x;
}

int Vertex::getY() const noexcept
{
    return y;
}