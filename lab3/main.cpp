#include "../include/Polinom.hpp"


using type_t = double;


static constexpr size_t NUMBER_OF_NODES = 5;


enum class Grid {
    UNIFORM = 0,
    CHEBYSHEVSKAYA
};


template<typename V>
V square(V x)
{
    // return x * x;
    // return 1 / (1 + x * x);
    return 1 / std::atan(1 + 10 * x * x);
}


template<typename V>
std::vector<V> createGrid(size_t number_of_nodes, std::pair<V, V> segment, Grid grid)
{
    switch (grid)
    {
    case Grid::UNIFORM:
    {
        std::vector<V> nodes(number_of_nodes);
        auto step = (std::abs(segment.first) + std::abs(segment.second)) / number_of_nodes;

        for (size_t i = 0; i < number_of_nodes; i++)
            nodes[i] = (segment.first + step * i);

        return nodes;
    }
    case Grid::CHEBYSHEVSKAYA:
    {
        std::vector<V> nodes(number_of_nodes);
        /*TODO*/
        return nodes;
    }
    default:
        break;
    }
    return std::vector<V>(number_of_nodes, 0);
}


int main()
{
    std::pair<type_t, type_t> segment(-3, 3);
    auto nodes = createGrid(NUMBER_OF_NODES, segment, Grid::UNIFORM);
    InterpolationPolinom<type_t> Lagrange(nodes, square<type_t>);
    auto polinom = Lagrange.getPolinomCoeffs();
    for (auto coeff : polinom)
        std::cout << coeff << " ";
    std::cout << std::endl;
}