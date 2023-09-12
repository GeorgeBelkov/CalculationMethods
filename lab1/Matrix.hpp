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

    // Вспомогательные методы класса
    void output();
    std::vector<T> getColumn(size_t index, size_t first_row_id = 0);
    T scalarMul(std::vector<T>& first, std::vector<T>& second);
    void transpoce();
    void mulRowAndScalar(size_t row_id, double scalar);
    void addRow(size_t src, size_t dest);
    std::pair<Matrix<T>, Matrix<T>> divideExtendedMatrix();
    // Методы по заданию
    bool forwardGaussStep();
    Matrix<T> backwardGaussStep();

    // Дружественные фунции
    template<typename V>
    friend Matrix<V> multiply(Matrix<V>& A, Matrix<V>& B);

    template<typename V>
    friend Matrix<V> makeE(size_t n);
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
            std::cout << elem << "\t";
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


// Умножение строки на скаляр.
template<typename T>
void Matrix<T>::mulRowAndScalar(size_t row_id, double scalar)
{
    for (auto& elem : matrix[row_id])
        elem *= scalar;
}


// Сложение строк (внутри матрицы)
template<typename T>
void Matrix<T>::addRow(size_t src, size_t dest)
{
    if (src < 0 or dest < 0)
    {
        std::cout << "Matrix hasn`t got below zero row or column!";
        exit(1);
    }
    
    auto& source = matrix[src];
    auto& destination = matrix[dest];
    for (size_t i = 0; i < source.size(); i++)
        destination[i] += source[i];
}


// Прямой ход Метода Гаусса
template<typename T>
bool Matrix<T>::forwardGaussStep()
{
    static double epsilon = 1e-10;
    size_t step = 0;
    size_t matrix_order = this->matrix.size();
    while (step < matrix_order)
    {
        auto column = this->getColumn(step, step);
        // Если "ведущий" элемент 0
        if (std::abs(column[0] - 0) < epsilon)
        {
            // Идем циклом по столбцу от [i][i] до [i][matrix_order]
            // Ищем ненулевой элемент
            for (size_t i = 0; i < column.size(); i++)
            {
                // Нашли
                if (std::abs(column[i] - 0) > epsilon)
                {
                    std::swap(this->operator[](step), this->operator[](step + i));
                    break;
                }
                // Не нашли --> выходим из функции
                if (i == column.size() - 1)
                {
                    std::cout << "Matrix determinate is zero!\n\n";
                    return true;
                }
            }
        }
        // Иначе
        else
        {
            T max = std::abs(column[0]);
            size_t iter = step;
            for (size_t i = 1; i < column.size(); i++)
            {
                if (std::abs(column[i]) > max)
                {
                    max = std::abs(column[i]);
                    iter = step + i;
                }
            }
            if (iter != step)
                std::swap(this->operator[](step), this->operator[](iter));
            
            
            // Складываем каждую нижележащую строку с первой умноженной на [i][step] / [step][step]
            for (size_t i = step + 1; i < matrix_order; i++)
            {

                std::vector<T> row_copy = this->operator[](step);

                if (this->matrix[i][step] == 0)
                    continue;
                T value = this->matrix[i][step] / this->matrix[step][step];

                this->mulRowAndScalar(step, -value);
                this->addRow(step, i);
                std::swap(row_copy, this->operator[](step));
            }
        }
        ++step;
    }
    return false;
}


// Обратный ход Метода Гаусса
template<typename T>
Matrix<T> Matrix<T>::backwardGaussStep()
{
    Matrix<T> Answer(this->matrix.size(), 1);
    for (size_t i = this->matrix.size() - 1; i >= 0; i--)
    {
        T sum = 0, b = this->matrix[i][matrix[0].size() - 1];
        for (size_t j = i + 1; j < this->matrix.size(); j++)
            sum += (this->matrix[i][j] * Answer[j][0]);
        Answer[i][0] = (b - sum) / this->matrix[i][i];

        if (i == 0) break;
    }
    return Answer;
}


// Возвращает пару матриц - матрицу коэффицентов и свободных членов
template<typename T>
std::pair<Matrix<T>, Matrix<T>> Matrix<T>::divideExtendedMatrix()
{
    Matrix<T> left(this->matrix.size(), this->matrix[0].size() - 1);
    Matrix<T> right(this->matrix.size(), 1);
    for (size_t i = 0; i < this->matrix.size(); i++)
    {
        for (size_t j = 0; j < this->matrix[0].size(); j++)
        {
            (j == this->matrix[0].size() - 1) ?
            right[i][0] = this->matrix[i][j] :
            left.matrix[i][j] = this->matrix[i][j];
        }
    }
    return std::make_pair(left, right);
}


// Возвращает столбец матрицы с индексом - index.
template<typename T>
std::vector<T> Matrix<T>::getColumn(size_t index, size_t first_row_id)
{
    if (index < matrix[0].size() && first_row_id < matrix.size() && first_row_id >= 0)
    {
        std::vector<T> column;
        column.reserve(matrix.size());
        for (; first_row_id < matrix.size(); first_row_id++)
            column.push_back(matrix[first_row_id][index]);
        return column;
    }
    else
    {
        std::cout << "index of column or first_row_id out of range!\n\n";
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

template<typename V>
Matrix<V> makeE(size_t n)
{
    Matrix<V> E(n,n);
    for (size_t i = 0; i < E.matrix.size(); i++)
    {
        E.matrix[i][i] = 1;
    }
    return E;
}