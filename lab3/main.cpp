#include "../include/Polinom.hpp"


using type_t = double;


enum class Task {
    MAIN_LAGRANGE_UNIFORM = 0,
    MAIN_LAGRANGE_CHEBYSHEV,
    MAIN_SPLINE_UNIFORM,
    MAIN_SPLINE_CHEBYSHEV,
    ADDITIONAL1,
    ADDITIONAL2,
};


template<typename V>
V complexFunction(V x)
{
    if (x > (-M_PI - 0.1) && x < 0)
    {
        return std::sin(x);
    }
    else if (x >= 0 && x <= (1 + 0.1))
    {
        return x - std::pow(x, 3) / 3;
    }
    
    
    return 5;
}


template<typename V>
V testFunction(V x)
{
    // return std::pow(x, 2);                                                                                                               //  +
    // return 1.0 / (1 + std::pow(x, 2));                                                                                                   //  +
    // return 1.0 / std::atan(1 + 10 * std::pow(x, 2));                                                                                     //  +
    // return std::pow(4 * std::pow(x, 3) + 2 * std::pow(x, 2) - 4 * x + 2, std::sqrt(2)) + std::asin(1 / (5 + x - std::pow(x, 2))) - 5;    //  +
    // return 1.0 / (1 + 25 * std::pow(x, 2));                                                                                                 +
    return std::sin(x);
    // return complexFunction(x);
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
void interpolant(std::function<V(V)>const& func, std::pair<V, V> segment)
{
    constexpr static size_t ITERS = 100000;
    static std::ofstream fout("functionInfo.txt");
    auto dist = std::abs(segment.first - segment.second) / (ITERS - 1);
    for (size_t i = 0; i < ITERS; i++)
        fout << (segment.first + i * dist) << " " << func(segment.first + i * dist) << "\n";

    fout.close();
}



int main()
{
    Task task = Task::MAIN_LAGRANGE_UNIFORM;
    std::string bash_command;
    switch (task)
    {
    case Task::MAIN_LAGRANGE_UNIFORM:
    {
        std::pair<type_t, type_t> uniform_segment(0, 10);
        Grid<type_t> uniform_grid(NODES_COUNT, uniform_segment, GridType::UNIFORM);

        interpolant<type_t>(testFunction<type_t>, uniform_segment);

        // Вывод получившейся сетки в файл.
        uniform_grid.print();

        InterpolationTable<type_t> uniform_grid_table(uniform_grid, testFunction<type_t>);

        // Вывод получившейся таблицы в файл.
        uniform_grid_table.printTable();

        Polinom<type_t> Lagrange_with_uniform_grid(NODES_COUNT - 1);
        LagrangeInterpolation<type_t>(Lagrange_with_uniform_grid, uniform_grid, uniform_grid_table);
        Lagrange_with_uniform_grid.print();
        
        bash_command = "python3 plotter.py lagrange";

        break;
    }
    case Task::MAIN_LAGRANGE_CHEBYSHEV:
    {
        std::pair<type_t, type_t> chebyshev_segment(-M_PI, 1);
        Grid<type_t> chebyshev_grid(NODES_COUNT, chebyshev_segment, GridType::CHEBYSHEVSKAYA);

        interpolant<type_t>(testFunction<type_t>, chebyshev_segment);

        // Вывод получившейся сетки в файл.
        chebyshev_grid.print();

        InterpolationTable<type_t> chebyshev_grid_table(chebyshev_grid, testFunction<type_t>);

        // Вывод получившейся таблицы в файл.
        chebyshev_grid_table.printTable();

        Polinom<type_t> Lagrange_with_chebyshev_grid(NODES_COUNT - 1);
        LagrangeInterpolation<type_t>(Lagrange_with_chebyshev_grid, chebyshev_grid, chebyshev_grid_table);
        Lagrange_with_chebyshev_grid.print();
        
        bash_command = "python3 plotter.py lagrange";

        break;
    }
    case Task::MAIN_SPLINE_UNIFORM:
    {
        std::pair<type_t, type_t> uniform_segment(-M_PI, 1);
        Grid<type_t> uniform_grid(NODES_COUNT, uniform_segment, GridType::UNIFORM);

        interpolant<type_t>(testFunction<type_t>, uniform_segment);

        // Вывод получившейся сетки в файл.
        uniform_grid.print();

        InterpolationTable<type_t> uniform_grid_table(uniform_grid, testFunction<type_t>);

        // Вывод получившейся таблицы в файл.
        uniform_grid_table.printTable();

        std::vector<Polinom<type_t>> Spline_with_uniform_grid;

        CubicSplineInterpolation(Spline_with_uniform_grid, uniform_grid_table);

        for (auto& polinom : Spline_with_uniform_grid)
           polinom.printSpline();

        bash_command = "python3 plotter.py spline";

        break;
    }
    case Task::MAIN_SPLINE_CHEBYSHEV:
    {
        std::pair<type_t, type_t> chebyshev_segment(-1, 1);
        Grid<type_t> chebyshev_grid(NODES_COUNT, chebyshev_segment, GridType::CHEBYSHEVSKAYA);

        interpolant<type_t>(testFunction<type_t>, chebyshev_segment);

        // Вывод получившейся сетки в файл.
        chebyshev_grid.print();

        InterpolationTable<type_t> chebyshev_grid_table(chebyshev_grid, testFunction<type_t>);

        // Вывод получившейся таблицы в файл.
        chebyshev_grid_table.printTable();

        std::vector<Polinom<type_t>> Spline_with_chebyshev_grid;

        CubicSplineInterpolation(Spline_with_chebyshev_grid, chebyshev_grid_table);

        for (auto& polinom : Spline_with_chebyshev_grid)
            polinom.printSpline();
        
        bash_command = "python3 plotter.py spline";

        break;
    }
    case Task::ADDITIONAL1:
    {
        // Флажок, 1 - если интерполяция сплайном проводится на равномерной сетке.
        bool is_grid_uniform = false;

        if (is_grid_uniform)
        {
            std::pair<type_t, type_t> uniform_segment(-100, 100);
            Grid<type_t> uniform_grid(NODES_COUNT, uniform_segment, GridType::UNIFORM);
            InterpolationTable<type_t> uniform_grid_table(uniform_grid, square<type_t>);

            interpolant<type_t>(square<type_t>, uniform_segment);

            // Вывод получившейся сетки в файл.
            uniform_grid.print();
            
            // Вывод получившейся таблицы в файл.
            uniform_grid_table.printTable();

            std::vector<Polinom<type_t>> Spline_with_uniform_grid;
            CubicSplineInterpolation(Spline_with_uniform_grid, uniform_grid_table);

            for (auto& polinom : Spline_with_uniform_grid)
                polinom.printSpline();

        }
        else
        {
            std::pair<type_t, type_t> chebyshev_segment(-5, 5);
            Grid<type_t> chebyshev_grid(NODES_COUNT, chebyshev_segment, GridType::CHEBYSHEVSKAYA);
            InterpolationTable<type_t> chebyshev_grid_table(chebyshev_grid, square<type_t>);

            interpolant<type_t>(square<type_t>, chebyshev_segment);

            // Вывод получившейся сетки в файл.
            chebyshev_grid.print();

            // Вывод получившейся таблицы в файл.
            chebyshev_grid_table.printTable();

            std::vector<Polinom<type_t>> Spline_with_chebyshev_grid;
            CubicSplineInterpolation(Spline_with_chebyshev_grid, chebyshev_grid_table);

            for (auto& polinom : Spline_with_chebyshev_grid)
                polinom.printSpline();

        }

        bash_command = "python3 plotter.py spline";

        break;
    }
    case Task::ADDITIONAL2:
    {
        std::pair<type_t, type_t> uniform_segment(-10, 10);
        Grid<type_t> uniform_grid(NODES_COUNT, uniform_segment, GridType::UNIFORM);
        InterpolationTable<type_t> uniform_grid_table(uniform_grid, constant<type_t>);

        interpolant<type_t>(constant<type_t>, uniform_segment);

        // Вывод получившейся сетки в файл.
        uniform_grid.print();

        // Вывод получившейся таблицы в файл.
        uniform_grid_table.printTable();

        Polinom<type_t> Lagrange_with_uniform_grid(NODES_COUNT - 1);

        LagrangeInterpolation<type_t>(Lagrange_with_uniform_grid, uniform_grid, uniform_grid_table);

        Lagrange_with_uniform_grid.print();

        bash_command = "python3 plotter.py lagrange";

        break;
    }
    default:
        break;
    }

    std::cout << "python!\n";
    std::system(bash_command.c_str());
}