#include "../include/Polinom.hpp"


using type_t = double;


enum class Task {
    MAIN = 0,
    ADDITIONAL1,
    ADDITIONAL2,
    ADDITIONAL3
};



template<typename V>
V testFunction(V x)
{
    return 1.0 / std::atan(1 + 10 * std::pow(x, 2));
    // return 1.0 / (1 + std::pow(x, 2));
    // return std::pow(x, 2);
}


template<typename V>
V square(V x)
{
    return std::pow(x, 2);
}


template<typename V>
V constant(V x)
{
    return 1;
}


template<typename V>
V thirdAdditionalTestFunction(V x)
{

}


int main()
{
    Task task = Task::MAIN;
    switch (task)
    {
    case Task::MAIN:
    {
        std::pair<type_t, type_t> uniform_segment(-1, 1);
        std::pair<type_t, type_t> chebyshev_segment(-1, 1);
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

        std::vector<Polinom<type_t>> Spline_with_uniform_grid;
        std::vector<Polinom<type_t>> Spline_with_chebyshev_grid;

        CubicSplineInterpolation(Spline_with_uniform_grid, uniform_grid_table);
        CubicSplineInterpolation(Spline_with_chebyshev_grid, chebyshev_grid_table);

        for (auto& polinom : Spline_with_uniform_grid)
            polinom.printSpline();

        for (auto& polinom : Spline_with_chebyshev_grid)
            polinom.printSpline();
        
        break;
    }
    case Task::ADDITIONAL1:
    {
        std::pair<type_t, type_t> uniform_segment(-5, 5);
        std::pair<type_t, type_t> chebyshev_segment(-5, 5);
        Grid<type_t> uniform_grid(NODES_COUNT, uniform_segment, GridType::UNIFORM);
        Grid<type_t> chebyshev_grid(NODES_COUNT, chebyshev_segment, GridType::CHEBYSHEVSKAYA);

        // Вывод получившихся сеток в файл.
        uniform_grid.print();
        chebyshev_grid.print();

        InterpolationTable<type_t> uniform_grid_table(uniform_grid, square<type_t>);
        InterpolationTable<type_t> chebyshev_grid_table(chebyshev_grid, square<type_t>);

        std::vector<Polinom<type_t>> Spline_with_uniform_grid;
        std::vector<Polinom<type_t>> Spline_with_chebyshev_grid;

        CubicSplineInterpolation(Spline_with_uniform_grid, uniform_grid_table);
        CubicSplineInterpolation(Spline_with_chebyshev_grid, chebyshev_grid_table);

        for (auto& polinom : Spline_with_uniform_grid)
            polinom.printSpline();

        for (auto& polinom : Spline_with_chebyshev_grid)
            polinom.printSpline();

        break;
    }
    case Task::ADDITIONAL2:
    {
        std::pair<type_t, type_t> uniform_segment(-10, 10);
        Grid<type_t> uniform_grid(NODES_COUNT, uniform_segment, GridType::UNIFORM);

        // Вывод получившихся сеток в файл.
        uniform_grid.print();

        InterpolationTable<type_t> uniform_grid_table(uniform_grid, constant<type_t>);

        Polinom<type_t> Lagrange_with_uniform_grid(NODES_COUNT - 1);

        LagrangeInterpolation<type_t>(Lagrange_with_uniform_grid, uniform_grid, uniform_grid_table);

        Lagrange_with_uniform_grid.print();

        break;
    }
    default:
        break;
    }
}