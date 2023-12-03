#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

using type_t = double;


constexpr static size_t NODES_COUNT = 15;


enum class Task {
    FIRST = 0,
    SECOND,
    THIRD,
    FOURTH,
    FIFTH
};


// Тестовые функции


template<typename T>
T firstTestFunc(T x)
{
    return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}


template<typename T>
T firstTestFuncDif(T x)
{
    return 5 * std::pow(x, 4) - 9.28 * std::pow(x, 3) + 5.9535 * std::pow(x, 2) - 1.5119 * x + 0.121495;
}


template<typename T>
T secondTestFunc(T x)
{
    return std::sqrt(x + 1) - 1;
}


template<typename T>
T secondTestFuncDif(T x)
{
    return 1 / (std::sqrt(x + 1) * 2);
}


template<typename T>
T thirdTestFunc(T x)
{
    return (35 * std::pow(x, 3) - 67 * std::pow(x, 2) - 3 * x + 3);
}


template<typename T>
T thirdTestFuncDif(T x)
{
    return (105 * std::pow(x, 2) - 134 * x - 3);
}


// Нужные


template<typename T>
void localization(std::vector<std::pair<T, T>>& local_segs, std::vector<T> const& grid, std::function<T(T)> func)
{
    size_t i = 0, j = i + 1;
    while (i < grid.size() - 1)
    {
        if (func(grid[i]) * func(grid[j]) < 0)
        {
            local_segs.push_back(std::make_pair(grid[j - 1], grid[j]));
            i = j;
            ++j;
        }
        else if (func(grid[i]) * func(grid[j]) == 0)
        {
            local_segs.push_back(std::make_pair(grid[j - 1], grid[j]));
            i = j + 1;
            j = i + 1;
        }
        else
        {
            ++i; ++j;
        }
    }
}


void makeGrid(std::vector<double>& nodes, std::pair<double, double> segment, size_t nodes_number)
{
    auto dist = std::abs(segment.first - segment.second) / (nodes_number - 1);
    int iter = 0;
    for (auto& node : nodes)
    {
        node = segment.first + iter * dist;
        iter++;
    }
}


template<typename T>
T bisectionsMethod(std::pair<T, T> segment, std::function<T(T)> const& func, size_t& iterations)
{
    static double epsilon = 1e-6;
    auto middle = (segment.first + segment.second) / 2;
    if (std::abs(segment.first - segment.second) < (epsilon * 2) || std::abs(func(middle)) < epsilon)
        return middle;

    iterations += 1;
    if (func(segment.first) * func(middle) < 0)
        return bisectionsMethod(std::make_pair(segment.first, middle), func, iterations);
    return bisectionsMethod(std::make_pair(middle, segment.second), func, iterations);
}


template<typename T>
T newton(std::function<T(T)>const& func, std::function<T(T)>const& func_dif, size_t& iterations, T first_approx)
{
    static double epsilon = 1e-6;
    T prev = first_approx, next, temp;
    do
    {
        if (std::abs(func(prev)) < epsilon)
            return prev;

        iterations++;
        temp = prev;
        next = prev - func(prev) / func_dif(prev);
        prev = next;
    } while (std::abs(next - temp) > epsilon);
    return next;
}


template<typename T>
T modifideNewton(std::function<T(T)>const& func, std::function<T(T)> const& func_dif,
                 size_t& iterations, T first_approx, std::pair<T, T> segment)
{
    static double epsilon = 1e-6;
    T prev = first_approx, next, temp;
    do
    {
        if (std::abs(func(prev)) < epsilon)
            return prev;
        
        iterations++;
        temp = prev;
        next = prev - func(prev) / func_dif(prev);
        if (next < segment.first || next > segment.second)
        {
            next = (func(segment.first) * segment.second - func(segment.second) * segment.first) / (func(segment.second) - func(segment.first));
            if (func(segment.first) * func(next) < 0)
                segment.second = next;
            else
                segment.first = next;
        }
        prev = next;
    } while (std::abs(next - temp) > epsilon);
    return next;
}


int main()
{

    std::function<type_t(type_t)> function;
    std::function<type_t(type_t)> funcDif;
    std::pair<type_t, type_t> segment;

    Task task = Task::SECOND;
    switch (task)
    {
    case Task::FIRST:
    {
        segment.first = 0;
        segment.second = 1;
        function = firstTestFunc<type_t>;
        funcDif = firstTestFuncDif<type_t>;
        break;
    }
    case Task::SECOND:
    {
        segment.first = -1;
        segment.second = 10;
        function = secondTestFunc<type_t>;
        funcDif = secondTestFuncDif<type_t>;
        break;
    }
    case Task::THIRD:
    {
        segment.first = 0;
        segment.second = 1;
        function = thirdTestFunc<type_t>;
        funcDif = thirdTestFuncDif<type_t>;
        break;
    }
    default:
        break;
    }


    std::vector<type_t> grid(NODES_COUNT);
    size_t newton_iters = 0, bisection_iters = 0, modified_newton_iters = 0;

    makeGrid(grid, segment, NODES_COUNT);

    for (auto& node : grid)
    {
        std::cout << node << " ";
    }
    std::cout << "\n\n";

    std::vector<std::pair<type_t, type_t>> local_segments;
    localization<type_t>(local_segments, grid, function);

    for (auto& pair : local_segments)
    {
        std::cout << "(" << pair.first << ";" << pair.second << ") ";
    }
    std::cout << "\n\n";


    std::cout << "Bisections:\n\n";
    for (auto& pair : local_segments)
    {
        static size_t id = 0;
        std::cout << "root " << id++ << " " << bisectionsMethod<type_t>(pair, function, bisection_iters)
                  << "; iters = " << bisection_iters << std::endl;
        bisection_iters = 0;
    }
    std::cout << "\nNewton:\n\n";
    for (auto& pair : local_segments)
    {
        static size_t id = 0;
        type_t approx;

        if (task == Task::FIRST)
            approx = (pair.second + pair.first) / 2;
        else if (task == Task::SECOND)
            approx = 8;
        else if (task == Task::THIRD)
        {
            // approx = 0;
            approx = (pair.second + pair.first) / 2;
        }
        
        std::cout << "root " << id++ << " " << newton<type_t>(function, funcDif, newton_iters, approx)
                  << "; iters = " << newton_iters << std::endl;
        newton_iters = 0;
    }
    std::cout << "\nModified Newton:\n\n";
    for (auto& pair : local_segments)
    {
        static size_t id = 0;
        auto approx = (pair.second + pair.first) / 2;
        std::cout << "root " << id++ << " "
                  << modifideNewton<type_t>(function, funcDif, modified_newton_iters, pair.second, segment)
                  << "; iters = " << modified_newton_iters << std::endl;
        modified_newton_iters = 0;
    }

}