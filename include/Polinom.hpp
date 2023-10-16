#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <filesystem>


template<typename T>
class InterpolationPolinom
{
private:
    std::vector<std::pair<T, T>> table;
    std::vector<T> polinom_coeffs;

public:

    // Конструктор
    InterpolationPolinom(std::vector<T> const& nodes, std::function<T(T)> function)
    {
        fillTable<T>(this->table, nodes, function);
        auto supportive_polinoms = makeSupportivePolinoms();
        polinom_coeffs.resize(nodes.size());

        size_t iter = 0;
        for (auto& polinom : supportive_polinoms)
        {
            for (auto coeff : polinom)
            {
                polinom_coeffs[iter] += (coeff * table[iter].second);
                iter++;
            }
            iter = 0;
        }
    }

    std::vector<T> getPolinomCoeffs() { return polinom_coeffs; }

    // для многочленов Лагранжа - создает вспомогательные полиномы C_{k}(x) и возвращает их вектором.
    std::vector<std::vector<T>> makeSupportivePolinoms();

    // Заполнение сетки pair.first - узлы, pair.second - значение сеточной функции в узле.
    template<typename V>
    friend void fillTable(std::vector<std::pair<V, V>>& table, std::vector<V> const& nodes, std::function<V(V)> function);

    // Вывод таблицы значений (x_i, y_i), по которым производится интерполирование.
    template<typename V>
    friend void outputTable(std::vector<std::pair<V, V>>& table);

    // Комбинирует индексы от 0 до polinom_deg с помощью ф-лы: C из polinom_deg по number
    template<typename V>
    friend std::vector<std::vector<size_t>> comb(short polinom_deg, short number);
};


template<typename V>
void fillTable(std::vector<std::pair<V, V>>& table, std::vector<V> const& nodes, std::function<V(V)> function)
{
    for (size_t i = 0; i < nodes.size(); i++)
    {
        std::pair<V, V> pair;
        pair.first = nodes[i];
        pair.second = function(nodes[i]);
        table.push_back(pair);
    }
}


template<typename V>
std::vector<std::vector<size_t>> comb(short polinom_deg, short number)
{
    std::vector<std::vector<size_t>> combinations;
    std::string bitmask(number, 1);
    bitmask.resize(polinom_deg + 1, 0);
    
    do {
        std::vector<size_t> temp;
        for (size_t i = 0; i <= polinom_deg; ++i)
        {
            if (bitmask[i]) temp.push_back(i);
        }
        combinations.push_back(temp);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    return combinations;
}


template<typename T>
std::vector<std::vector<T>> InterpolationPolinom<T>::makeSupportivePolinoms()
{
    std::vector<std::vector<T>> polinoms;

    for (size_t i = 0; i < this->table.size(); i++)
    {
        std::vector<T> polinom(this->table.size(), 1);
        T multiplier = 1;
        
        for (size_t j = 0; j < this->table.size(); j++)
            if (j != i)
                multiplier *= (this->table[i].first - this->table[j].first);

        for (size_t j = 0; j < this->table.size(); j++)
        {
            if (j == i)
                continue;
            else
            {
                auto temp = comb<T>(this->table.size(), j);
                T sum = 0;
                for (auto& indexes : temp)
                {
                    T mul = 1;
                    for (auto index : indexes)
                    {
                        mul *= this->table[index].first; 
                    }
                    sum += mul;
                    mul = 1;
                }
                polinom[j] = sum;
            }
        }

        for (auto& coeff : polinom)
            coeff *= multiplier;

        polinoms.push_back(polinom);
    }
    
    return polinoms;
}









