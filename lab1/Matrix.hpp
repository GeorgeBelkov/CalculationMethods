#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <algorithm>


template<typename T>
class Matrix
{
private:
    std::vector<std::vector<T>> matrix;
    // Заполнение матрицы
    void fill(std::filesystem::path& filename);

public:
    ~Matrix() = default;

    // Конструкторы
    Matrix(size_t n, size_t m, std::filesystem::path filename);
    Matrix(size_t n, size_t m);
    Matrix(Matrix<T>& other);

    // Оператор присваивания
    Matrix<T>& operator=(Matrix<T>& other);
    // Перегрузка операторов
    std::vector<T>& operator[](size_t index);

    // Геттеры и Сеттеры
    std::vector<std::vector<T>>& getMatrix() { return matrix; }

    // Методы класса
    void output();
    std::vector<T> getColumn(size_t index);
    T scalarMul(std::vector<T>& first, std::vector<T>& second);
    void transpoce();
    // Дописать
    // void addRow(int dest_id, int src_id, uint koeff1, uint koeff2);


    // Дружественные фунции
    template<typename V>
    friend Matrix<V> multiply(Matrix<V>& A, Matrix<V>& B);
};


// Полная специализация шаблонного метода fill для типа float
template<>
void Matrix<float>::fill(std::filesystem::path& filename)
{
    float number;
    std::string str;
    size_t row = 0, column = 0;
    std::ifstream file(filename, std::ios_base::in);
    while (!file.eof())
    {
        file >> number;
        matrix[row][column++] = number;
        if (column == matrix[0].size())
        {
            column = 0;
            row++;
        }
    }
}


// Полная специализация шаблонного метода fill для типа double
template<>
void Matrix<double>::fill(std::filesystem::path& filename)
{
    double number;
    size_t row = 0, column = 0;
    std::ifstream file(filename, std::ios_base::in);
    while (!file.eof())
    {
        file >> std::setprecision(7) >> number;
        std::swap(matrix[row][column++], number);
        if (column == matrix[0].size())
        {
            column = 0;
            row++;
        }
    }
}


// Конструктор. Заполняет матрицу значениями из файла.
template<typename T>
Matrix<T>::Matrix(size_t n, size_t m, std::filesystem::path filename)
{
    matrix.reserve(n);
    for (size_t i = 0; i < n; i++)
    {
        std::vector<T> vector(m);
        matrix.push_back(vector);
    }
    fill(filename);
}


// Конструктор "по умолчанию". В матрице нет значений.
template<typename T>
Matrix<T>::Matrix(size_t n, size_t m)
{
    matrix.reserve(n);
    for (size_t i = 0; i < n; i++)
    {
        std::vector<T> vector(m);
        matrix.push_back(vector);
    }
}


// Конструктор копирования
template<typename T>
Matrix<T>::Matrix(Matrix& other)
{
    matrix.reserve(other.matrix.size());
    for (size_t i = 0; i < other.matrix.size(); i++)
    {
        matrix.push_back(std::vector<T>(other.matrix[i]));
    }
}


// Вывод матрицы в консоль.
template<typename T>
void Matrix<T>::output()
{
    for (auto& row : matrix)
    {
        for (T elem : row)
            std::cout << elem << " ";
        std::cout << "\n";
    }
    std::cout << std::setprecision(7) << std::endl;
}


// Дает возможность получать строку матрицы A, как A[i], i - индекс строки.
template<typename T>
std::vector<T>& Matrix<T>::operator[](size_t index)
{
    if (index < matrix.size())
        return matrix[index];
    else
    {
        std::cout << "index of row out of range!\n\n";
        exit(1);
    }
}


// Скалярное умножение
template<typename T>
T Matrix<T>::scalarMul(std::vector<T>& first, std::vector<T>& second)
{
    T result = 0;
    for (size_t i = 0; i < first.size(); i++)
    {
        result += (first[i] * second[i]);
    }
    return result;
}


// Возвращает столбец матрицы с индексом - index.
template<typename T>
std::vector<T> Matrix<T>::getColumn(size_t index)
{
    if (index < matrix[0].size())
    {
        std::vector<T> column;
        column.reserve(matrix.size());
        for (size_t i = 0; i < matrix.size(); i++)
            column.push_back(matrix[i][index]);
        return column;
    }
    else
    {
        std::cout << "index of column out of range!\n\n";
        exit(1);
    }
}


// Транспонирует матрицу. B.transpoce() => B^(T)
template<typename T>
void Matrix<T>::transpoce()
{
    auto temp = matrix;

    matrix.resize(temp[0].size());
    for (auto& row : matrix)
        row.resize(temp.size());
    
    for (size_t i = 0; i < temp.size(); i++)
        for (size_t j = 0; j < temp[0].size(); j++)
            matrix[i][j] = temp[j][i];
}


// Оператор присваивания.
template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix& other)
{
    bool is_different;
    int counter = 0;
    for (size_t i = 0; i < this->matrix.size(); i++)
    {
        if (this->matrix[i] == other.matrix[i])
            counter++;
    }

    if (counter == this->matrix.size())
        is_different = false;
    else
        is_different = true;


    if (is_different &&
        this->matrix.size() == other.matrix.size() &&
        this->matrix[0].size() == other.matrix[0].size())
    {
        for (size_t i = 0; i < this->matrix.size(); i++)
        {
            this->matrix[i].clear();
            this->matrix[i] = other.matrix[i];
        }
    }
    else
        std::cout << "They are not different!!!\n\n";
    
    return *this;
}



// Дружественные функции - реализация



// Умножение матриц
template<typename V>
Matrix<V> multiply(Matrix<V>& A, Matrix<V>& B)
{
    if (A.matrix[0].size() == B.matrix.size())
    {
        Matrix<V> C(A.matrix.size(), B.matrix[0].size());
        for (size_t i = 0; i < C.matrix.size(); i++)
            for (size_t j = 0; j < C.matrix[0].size(); j++)
            {
                auto column = B.getColumn(j);
                C.matrix[i][j] = C.scalarMul(A[i], column);
            }
        return C;
    }
    else
    {
        std::cout << "dimentions of A and B are different!\n\n";
        exit(1);
    }
}