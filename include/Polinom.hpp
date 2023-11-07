#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <numeric>
#include <iostream>
#include <algorithm>


constexpr static size_t NODES_COUNT = 5;


template<typename T>
class Grid;


enum class GridType {
    UNIFORM = 0,
    CHEBYSHEVSKAYA
};


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
        for (auto& pair : table)
            table_fout << "x: " << pair.first << "; y: " << pair.second << "\n";
        table_fout << "\n\n";
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


    // Дружественные функции класса

    // Комбинирует индексы от 0 до polinom_deg с помощью ф-лы: C из polinom_deg по number
    template<typename V>
    friend V getCoeff(Grid<V> const& grid, uint16_t polinom_deg, uint16_t number);

    template<typename V>
    friend void LagrangeInterpolation(Polinom<V>& Lagrange_polinom, Grid<V> const& grid, InterpolationTable<V> const& table);

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


template<typename T>
void Polinom<T>::print() const
{
    static std::ofstream fout("InterpolationInfo.txt");

    for (auto& coeff : polinom_coeffs)
        fout << coeff << " ";

    fout << "\n\n";
}


template<typename V>
V getCoeff(Grid<V> const& grid, uint16_t polinom_deg, uint16_t number)
{
    // Хранит всевозможные комбинации произведений длины number
    std::vector<V> combinations;

    // Битмаска
    std::string bitmask(number, 1);
    bitmask.resize(polinom_deg + 1, 0);
    
    do {
        V combination = 1;
        for (size_t i = 0; i <= polinom_deg; ++i)
            if (bitmask[i]) combination *= grid.getGrid()[i];

        combinations.push_back(combination);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    // Возвращаем сумму (то есть коэффицент) при (polinom_deg - number) степени.
    return std::accumulate(combinations.begin(), combinations.end(), 1);
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
        
        V multiplier = 1;
        for (size_t i = 0; ((i != k) && (i < polinom_deg)); i++)
            multiplier *= (grid_copy[k] - grid_copy[i]);
        

        // заполняем ее коэффиценты
        base_function.polinom_coeffs[0] *= multiplier * table.table[polinom_deg - 1].second;
        base_function.polinom_coeffs[1] = ((-multiplier) * std::accumulate(grid_copy.begin(), grid_copy.end(), 1) * table.table[polinom_deg - 2].second);
        for (size_t i = 2; i < polinom_deg; i++)
            // i-й коэффицент = (-1)^i * множитель * сумму комбинаций.
            base_function.polinom_coeffs[i] = std::pow(-1, i) * multiplier * getCoeff<V>(grid, polinom_deg, i) * table.table[polinom_deg - (i + 1)].second;
        
        for (size_t i = 0; i < Lagrange_polinom.polinom_coeffs.size(); i++)
            Lagrange_polinom.polinom_coeffs[i] += base_function.polinom_coeffs[i];
    }
}