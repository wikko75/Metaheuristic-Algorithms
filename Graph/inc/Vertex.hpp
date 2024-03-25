#ifndef VERTEX_1_HPP
#define VERTEX_1_HPP

class Vertex
{
private:
    int number;
    int x;
    int y;

public:
    Vertex(int n, int xCoord, int yCoord);

    int getNumber() const noexcept;

    int getX() const noexcept;

    int getY() const noexcept;

    ~Vertex() = default;
};

#endif