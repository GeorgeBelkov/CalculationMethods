#include "../include/Polinom.hpp"


using type_t = double;
constexpr static size_t NODES_COUNT = 5;


int main()
{
    std::pair<type_t, type_t> uniform_segment(1, 3);
    std::pair<type_t, type_t> chebyshev_segment(1, 3);
    Grid<type_t> uniform_grid(NODES_COUNT, uniform_segment, GridType::UNIFORM);
    Grid<type_t> chebyshev_grid(NODES_COUNT, chebyshev_segment, GridType::CHEBYSHEVSKAYA);


    uniform_grid.print();
    chebyshev_grid.print();
}