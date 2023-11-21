#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <numeric>
#include <iostream>
#include <algorithm>


constexpr static size_t NODES_COUNT = 200;
constexpr static double epsilon = 1e-8;


template<typename T>
class Grid;

template<typename T>
struct InterpolationTable;


enum class GridType {
    UNIFORM = 0,
    CHEBYSHEVSKAYA
};


// Функция возвращает прогоночные коэффиценты alpha[i], betta[i] для нахождения c_i в кубическом сплайне.
template<typename V>
std::pair<V, V> getRunThroughCoeffs(InterpolationTable<V> const& table, std::vector<std::pair<V, V>> const& coeffs, size_t index)
{
    // можно оптимальнее (избавиться от двух аддитивных и двух мультипликативных)
    // g_i and g_(i-1)
    auto g_i = (table.table[index].second - table.table[index - 1].second) / (table.table[index].first - table.table[index - 1].first);
    auto g_i_prev = (table.table[index - 1].second - table.table[index - 2].second) / (table.table[index - 1].first - table.table[index - 2].first);
    if (index == 2)
    {
        auto alpha = (table.table[index].first - table.table[index - 1].first) / (-2 * (table.table[index].first - table.table[index - 2].first));
        auto betta = (3 * (g_i - g_i_prev)) / (2 * (table.table[index].first - table.table[index - 2].first));

        return std::make_pair(alpha, betta);
    }
    else if (index == table.table.size() - 1)
    {
        auto betta = (-3 * (g_i - g_i_prev) + (table.table[index - 1].first - table.table[index - 2].first) * coeffs[index - 3].second) /
                     (-2 * (table.table[index].first - table.table[index - 2].first) -
                      (table.table[index - 1].first - table.table[index - 2].first) * coeffs[index - 3].first);
        return std::make_pair(0, betta);
    }
    else
    {
        auto alpha = (table.table[index].first - table.table[index - 1].first) /
                     (-2 * (table.table[index].first - table.table[index - 2].first) -
                      (table.table[index - 1].first - table.table[index - 2].first) * coeffs[index - 3].first);

        auto betta = (-3 * (g_i - g_i_prev) + (table.table[index - 1].first - table.table[index - 2].first) * coeffs[index - 3].second) /
                     (-2 * (table.table[index].first - table.table[index - 2].first) +
                      -1 * (table.table[index - 1].first - table.table[index - 2].first) * coeffs[index - 3].first);
        return std::make_pair(alpha, betta);
    }
}


template<typename T>
struct InterpolationTable
{
    std::vector<std::pair<T, T>> table;

    InterpolationTable(Grid<T> const& grid, std::function<T(T)> const& function)
    {
        table.resize(grid.getGrid().size());
        for (size_t i = 0; i < table.size(); i++)
            table[i] = std::make_pair(grid.getGrid()[i], function(grid.getGrid()[i]));
    }

    void printTable() const
    {
        static std::ofstream table_fout("tableInfo.txt");
        static bool is_closed = false;

        if (is_closed)
            table_fout.open("tableInfo.txt", std::ios::out);

        for (auto& pair : table)
            table_fout << "x: " << pair.first << " y: " << pair.second << "\n";
        table_fout << "\n\n";

        table_fout.close();
        is_closed = true;
    }
};


template<typename T>
class Grid
{
private:
    std::vector<T> grid;
    GridType grid_type;
public:
    // Конструкторы
    Grid(size_t nodes_count, std::pair<T, T> segment, GridType grid_tp);
    Grid(Grid<T> const& other): grid(other.grid), grid_type(other.grid_type) {}
    Grid(Grid<T>&& other) noexcept: grid(std::move(other.grid)), grid_type(other.grid_type) {}

    // Деструкторы
    ~Grid() = default;
    
    // Оператор присваивания
    Grid<T>& operator=(Grid<T> const& other);
    // Оператор перемещения
    Grid<T>& operator=(Grid<T>&& other) noexcept;

    // Методы класса
    void print() const;

    // Get, Set методы
    std::vector<T> getGrid() const { return grid; }
};


template<typename T>
Grid<T>::Grid(size_t nodes_count, std::pair<T, T> segment, GridType grid_tp)
{
    grid.resize(nodes_count);
    grid_type = grid_tp;
    switch (grid_tp)
    {
    // Если нужна равномерная сетка
    case GridType::UNIFORM:
    {
        auto dist = std::abs(segment.first - segment.second) / (nodes_count - 1);
        uint16_t iter = 0;
        for (auto& node : grid)
        {
            node = segment.first + iter * dist;
            iter++;
        }
        break;
    }
    // Если нужна Чебышевская сетка
    case GridType::CHEBYSHEVSKAYA:
    {
        std::vector<T> chebyshev_nodes;
        for (size_t k = 0; k < nodes_count; k++)
        {
            auto chebyshev_node = std::cos(((2 * (k + 1) - 1) * M_PI) / (2 * nodes_count));
            grid[k] = (0.5 * (segment.first + segment.second) + 0.5 * (segment.second - segment.first) * chebyshev_node);
        }
        std::sort(grid.begin(), grid.end());
        break;
    }
    default:
    {
        std::cout << "\nChoose correct grid type!\n";
        exit(1);
    }
    }
}


template<typename T>
Grid<T>& Grid<T>::operator=(Grid<T> const& other)
{
    grid.clear();
    grid = other.grid;
    grid_type = other.grid_type;
    return *this;
}


template<typename T>
Grid<T>& Grid<T>::operator=(Grid<T>&& other) noexcept
{
    grid.clear();
    grid = std::move(other.grid);
    grid_type = other.grid_type;
    return *this;
}


// Методы класса


template<typename T>
void Grid<T>::print() const
{
    static std::ofstream fout("gridInfo.txt");

    if (grid_type == GridType::CHEBYSHEVSKAYA)
        fout << "Chebyshevskaya grid:\n";
    else
        fout << "Uniform grid:\n";
    
    for (auto& node : grid)
        fout << node << " ";

    fout << "\n\n";
}


template<typename T>
class Polinom
{
private:
    std::vector<T> polinom_coeffs;

public:
    // Конструкторы
    Polinom() = default;
    explicit Polinom(size_t deg) { polinom_coeffs.reserve(deg + 1); }
    Polinom(Polinom<T> const& other): polinom_coeffs(other.polinom_coeffs) {}
    Polinom(Polinom<T>&& other) noexcept: polinom_coeffs(std::move(other.polinom_coeffs)) {}

    // Деструкторы
    ~Polinom() = default;
    
    // Оператор присваивания
    Polinom<T>& operator=(Polinom<T> const& other);
    // Оператор перемещения
    Polinom<T>& operator=(Polinom<T>&& other) noexcept;


    // Методы класса
    void print() const;
    void printSpline() const;


    // Дружественные функции класса

    // Комбинирует индексы от 0 до polinom_deg с помощью ф-лы: C из polinom_deg по number
    template<typename V>
    friend V getCoeff(Grid<V> const& grid, uint16_t polinom_deg, uint16_t number, size_t excluded);

    template<typename V>
    friend void LagrangeInterpolation(Polinom<V>& Lagrange_polinom, Grid<V> const& grid, InterpolationTable<V> const& table);

    template<typename V>
    friend void LagrangeInterpolationFake(Grid<V> const& grid, InterpolationTable<V> const& table);

    template<typename V>
    friend void CubicSplineInterpolation(std::vector<Polinom<V>>& spline, InterpolationTable<V> const& table);

};


template<typename T>
Polinom<T>& Polinom<T>::operator=(Polinom<T> const& other)
{
    polinom_coeffs.clear();
    polinom_coeffs(other.polinom_coeffs);
    return *this;
}


template<typename T>
Polinom<T>& Polinom<T>::operator=(Polinom<T>&& other) noexcept
{
    polinom_coeffs.clear();
    polinom_coeffs = std::move(other.polinom_coeffs);
    return *this;
}


template<typename T>
void Polinom<T>::print() const
{
    static std::ofstream fout("LagrangeInterpolationInfo.txt");
    static bool is_closed = false;

    if (is_closed)
        fout.open("LagrangeInterpolationInfo.txt", std::ios::out);
    

    for (size_t i = 0; i < polinom_coeffs.size(); i++)
    {
        if (i == polinom_coeffs.size() - 1)
            fout << polinom_coeffs[i];
        else
            fout << polinom_coeffs[i] << " ";
    }

    fout << "\n\n";

    fout.close();
    is_closed = true;
}


template<typename T>
void Polinom<T>::printSpline() const
{
    static std::ofstream fout("SplineInterpolationInfo.txt");
    static size_t id = 1;
    static bool is_closed = false;

    if (is_closed)
        fout.open("SplineInterpolationInfo.txt", std::ios::out);

    fout << id << ": ";
    for (size_t i = 0; i < polinom_coeffs.size(); i++)
    {
        if (i == polinom_coeffs.size() - 1)
            fout << polinom_coeffs[i];
        else
            fout << polinom_coeffs[i] << " ";
    }
    
    fout << "\n";

    id++;
    if (id == NODES_COUNT)
    {
        id = 1;
        fout << "\n\n";
        fout.close();
        is_closed = true;
    }
}


template<typename V>
V getCoeff(Grid<V> const& grid, uint16_t polinom_deg, uint16_t number, size_t excluded)
{
    // Хранит всевозможные комбинации произведений длины number
    std::vector<V> combinations;

    // Битмаска
    std::string bitmask(number, 1);
    bitmask.resize(polinom_deg + 1, 0);
    auto grid_copy = grid.getGrid();
    
    do {
        V combination = 1;
        bool push_flag = true;
        for (size_t i = 0; i <= polinom_deg; ++i)
            if (bitmask[i])
            {
                if (i == excluded)
                {
                    push_flag = false;
                    break;
                }
                combination *= grid_copy[i];
            }

        if (push_flag)
            combinations.push_back(combination);
        

    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    // Возвращаем сумму (то есть коэффицент) при (polinom_deg - number) степени.
    V sum = 0;
    for (auto& elem : combinations)
        if (std::abs(elem) > epsilon)
            sum += elem;
    return sum;
}


template<typename V>
void LagrangeInterpolation(Polinom<V>& Lagrange_polinom, Grid<V> const& grid, InterpolationTable<V> const& table)
{
    size_t polinom_deg = table.table.size();
    Lagrange_polinom.polinom_coeffs.resize(polinom_deg);
    auto grid_copy = grid.getGrid();
    // Хранит базис, по которому строится полином Лагранжа
    std::vector<Polinom<V>> basic_functions(polinom_deg);

    for (size_t k = 0; k < polinom_deg; k++)
    {
        // берем k-ю базисную функцию
        auto& base_function = basic_functions[k];
        base_function.polinom_coeffs.resize(polinom_deg, 1);
        
        V multiplier = 1, sum = 0;
        for (size_t i = 0; i < polinom_deg; i++)
            if (i != k)
            {
                multiplier *= (grid_copy[k] - grid_copy[i]);
                sum += grid_copy[i];
            }
        

        // заполняем ее коэффиценты
        base_function.polinom_coeffs[0] *= (1 / multiplier) * table.table[k].second;
        base_function.polinom_coeffs[1] *= (-(1 / multiplier) * sum * table.table[k].second);
        for (size_t i = 2; i < polinom_deg; i++)
            // i-й коэффицент = (-1)^i * множитель * сумму комбинаций.
            base_function.polinom_coeffs[i] *= std::pow(-1, i) * (1 / multiplier) * getCoeff<V>(grid, polinom_deg, i, k) * table.table[k].second;
        
        for (size_t i = 0; i < polinom_deg; i++)
            Lagrange_polinom.polinom_coeffs[i] += base_function.polinom_coeffs[i];
    }
    for (auto& elem : Lagrange_polinom.polinom_coeffs)
        if (std::abs(elem) < epsilon)
            elem = 0;
}


template<typename V>
void LagrangeInterpolationFake(Grid<V> const& grid, InterpolationTable<V> const& table)
{
    constexpr static size_t ITERS = 1000;
    static std::ofstream fout("FakeLagrangeInfo.txt");
    auto dist = std::abs(grid.getGrid()[0] - grid.getGrid()[grid.getGrid().size() - 1]) / (ITERS - 1);
    for (size_t i = 0; i < ITERS; i++)
    {
        auto x = (grid.getGrid()[0] + i * dist);
        std::vector<V> base_funcktions_results;
        size_t iterations = 0;
        for (auto& node : grid.getGrid())
        {
            V numerator = 1, denumerator = 1;
            for (size_t j = 0; j < grid.getGrid().size(); j++)
            {
                if (grid.getGrid()[j] != node)
                {
                    numerator *= (x - grid.getGrid()[j]);
                    denumerator *= (node - grid.getGrid()[j]);
                }
            }
            base_funcktions_results.push_back((numerator / denumerator) * table.table[iterations].second);
            // std::cout << base_funcktions_results[iterations] << " ";
            ++iterations;
        }
        // std::cout << "\n";

        V y = 0;
        for (auto& elem : base_funcktions_results)
            y += elem;
        
        fout << x << " " << y << "\n";
    }
    fout.close();
}


template<typename V>
void CubicSplineInterpolation(std::vector<Polinom<V>>& spline, InterpolationTable<V> const& table)
{
    std::vector<std::pair<V, V>> coeffs;
    for (size_t i = 2; i < table.table.size(); i++)
        coeffs.push_back(getRunThroughCoeffs(table, coeffs, i));

    short j = -1;
    for (size_t i = table.table.size() - 1; i > 0; i--, j++)
    {
        Polinom<V> cubic_polinom;
        cubic_polinom.polinom_coeffs.resize(4);

        auto g_i = (table.table[i].second - table.table[i - 1].second) / (table.table[i].first - table.table[i - 1].first);
        auto h_i = table.table[i].first - table.table[i - 1].first;

        V a_i = table.table[i - 1].second, c_i, b_i, d_i;
        if (i == table.table.size() - 1)
        {
            c_i = coeffs[i - 2].second;
            b_i = g_i - (2 * c_i * h_i) / 3;
            d_i = (-c_i) / (3 * h_i);
        }
        else
        {
            c_i = coeffs[i - 2].first * spline[j].polinom_coeffs[1] + coeffs[i - 2].second;
            b_i = g_i - ((spline[j].polinom_coeffs[1] + 2 * c_i) * h_i) / 3;
            d_i = (spline[j].polinom_coeffs[1] - c_i) / (3 * h_i);
        }

        // устанавливаем коэффиценты полинома
        cubic_polinom.polinom_coeffs[0] = d_i;
        cubic_polinom.polinom_coeffs[1] = c_i;
        cubic_polinom.polinom_coeffs[2] = b_i;
        cubic_polinom.polinom_coeffs[3] = a_i;

        for (auto& elem : cubic_polinom.polinom_coeffs)
            if (std::abs(elem) < epsilon)
                elem = 0;

        spline.push_back(cubic_polinom);
    }
    std::reverse(spline.begin(), spline.end());
}