#include "../include/Polinom.hpp"


using type_t = double;


template<typename V>
V testFunction(V x)
{
    return 1.0 / std::atan(1 + 10 * std::pow(x, 2));
    // return 1.0 / (1 + std::pow(x, 2));
    // return std::pow(x, 2);
}


int main()
{
    std::pair<type_t, type_t> uniform_segment(-3, 3);
    std::pair<type_t, type_t> chebyshev_segment(-3, 3);
    Grid<type_t> uniform_grid(NODES_COUNT, uniform_segment, GridType::UNIFORM);
    Grid<type_t> chebyshev_grid(NODES_COUNT, chebyshev_segment, GridType::CHEBYSHEVSKAYA);

    // Вывод получившихся сеток в файл.
    uniform_grid.print();
    chebyshev_grid.print();

    InterpolationTable<type_t> uniform_grid_table(uniform_grid, testFunction<type_t>);
    InterpolationTable<type_t> chebyshev_grid_table(chebyshev_grid, testFunction<type_t>);

    // Вывод получившихся таблиц в файл.
    uniform_grid_table.printTable();
    chebyshev_grid_table.printTable();

    Polinom<type_t> Lagrange_with_uniform_grid(NODES_COUNT - 1);
    Polinom<type_t> Lagrange_with_chebyshev_grid(NODES_COUNT - 1);

    LagrangeInterpolation<type_t>(Lagrange_with_uniform_grid, uniform_grid, uniform_grid_table);
    LagrangeInterpolation<type_t>(Lagrange_with_chebyshev_grid, chebyshev_grid, chebyshev_grid_table);

    Lagrange_with_uniform_grid.print();
    Lagrange_with_chebyshev_grid.print();
}