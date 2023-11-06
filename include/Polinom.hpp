#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>


enum class GridType {
    UNIFORM = 0,
    CHEBYSHEVSKAYA
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
            auto chebyshev_node = std::cos(((2 * k - 1) * M_PI) / (2 * nodes_count));
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
    explicit Polinom(size_t deg) { polinom_coeffs.reserve(deg + 1); }
    Polinom(Polinom<T> const& other): polinom_coeffs(other.polinom_coeffs) {}
    Polinom(Polinom<T>&& other) noexcept: polinom_coeffs(std::move(other.polinom_coeffs)) {}

    // Деструкторы
    ~Polinom() = default;
    
    // Оператор присваивания
    Polinom<T>& operator=(Polinom<T> const& other);
    // Оператор перемещения
    Polinom<T>& operator=(Polinom<T>&& other) noexcept;


    // Дружественные функции класса

    template<typename V>
    friend Polinom<V> LarangeInterpolation();

    template<typename V>
    friend std::vector<Polinom<V>> SplineInterpolation();

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