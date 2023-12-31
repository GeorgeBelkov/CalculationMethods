#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <cmath>
#include <fstream>



constexpr static double epsilon = 10e-8;
constexpr static double ITER_PARAM = 0.0072;
constexpr static double RELAX_PARAM = 0.18;
constexpr static size_t MATRIX_ORD = 201;


enum class IterativeMethod {
    SIMPLE_ITER = 0,
    JACOBI,
    SEIDEL,
    RELAXATION
};


enum class ExitQRCriteria {
    SIMPLE = 0,
    MODULE,
    SQUARE_MODULE
};


enum class Norm {
    octahedral = 1,
    spherical,
    cubic
};


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
    Matrix(Matrix<T>&& other) noexcept: matrix(std::move(other.matrix)) {}

    // Операторы присваивания
    Matrix<T>& operator=(Matrix<T>& other);
    Matrix<T>& operator=(Matrix<T>&& other) noexcept;

    // Перегрузка операторов
    std::vector<T>& operator[](size_t index);
    Matrix<T> operator+(Matrix<T>& B);
    Matrix<T> operator-(Matrix<T>& B);

    // Геттеры и Сеттеры
    std::vector<std::vector<T>>& getMatrix() { return matrix; }

    // Вспомогательные методы класса
    void output(std::ofstream& fout);
    std::vector<T> getColumn(size_t index, size_t first_row_id = 0);
    T scalarMul(std::vector<T>& first, std::vector<T>& second);
    T scalarMul(Matrix<T>& first, Matrix<T>& second);
    T vectorLength(Matrix<T>& vector);
    void transpoce();
    void mulRowAndScalar(size_t row_id, T scalar);
    void mulMatrixAndScalar(T scalar);
    void addRow(size_t src, size_t dest);
    std::pair<Matrix<T>, Matrix<T>> makeMatrixGoodForIter(IterativeMethod method, std::pair<Matrix<T>, Matrix<T>> divided);
    std::pair<Matrix<T>, Matrix<T>> divideExtendedMatrix();
    Matrix<T> reduction(size_t reduction_deg);
    Matrix<T> Unv();
    std::pair<T,T> Cond();

    // Методы по заданию
    // lab1
    bool forwardGaussStep();
    Matrix<T> backwardGaussStep();
    // lab2
    Matrix<T> simpleIterationsMethod(Matrix<T> initial_approx);
    Matrix<T> methodJacobi(Matrix<T> initial_approx);
    // lab4
    Matrix<T> iterativeQR(Matrix<T> A);
    void makeHessenbergForm(Matrix<T>& A);
    std::pair<Matrix<T>, Matrix<T>> QR(Matrix<T>& A);
    std::pair<Matrix<T>, T> reverseIterationsMethod(bool rayleigh, Matrix<T>& A, Matrix<T> initial_approx, T lambda, size_t* iter = nullptr);

    // Дружественные фунции

    template<typename V>
    friend void exitQRCriteria(Matrix<V>& A, bool& continue_flag, size_t const size, ExitQRCriteria condition);

    template<typename K>
    friend K eqlidNorm(Matrix<K> a);

    template<typename V>
    friend V octahedralNorm(Matrix<V> a);
    
    template<typename V>
    friend V cubicNorm(Matrix<V> a);

    template<typename V>
    friend Matrix<V> makeExtendedMatrix(Matrix<V>& A, Matrix<V>& b);

    template<typename V>
    friend Matrix<V> multiply(Matrix<V>& A, Matrix<V>& B);

    template<typename V>
    friend Matrix<V> makeE(size_t n);

    template<typename V>
    friend Matrix<V> makeRot(size_t num1, size_t num2, Matrix<V>& mat);
    
    template<typename V>
    friend bool QRforSLAU(Matrix<V>& A,Matrix<V>& Q, Matrix<V>& R,Matrix<V>& Uns);

    template<typename K>
    friend void Rotation(Matrix<K>& A, size_t first_id, size_t second_id, Matrix<K>& B, bool transpose);

    template<typename K>
    friend void RotationForHess(Matrix<K>& A, size_t first_id, size_t second_id, Matrix<K>& B, bool transpose);

    template<typename V>
    friend Matrix<V> pgr(V delta,Matrix<V> b);

    template<typename V>
    friend std::pair<Matrix<V>, double> relaxationMethod(Matrix<V> initial_approx, Matrix<V> U, Matrix<V> D, Matrix<V> L, Matrix<V> b, IterativeMethod method);
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
    file.close();
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
    file.close();
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


// Оператор перемещаения
template<typename T>
Matrix<T>&  Matrix<T>::operator=(Matrix<T>&& other) noexcept
{
    this->matrix = std::move(other.matrix);
    return *this;
}


// Вывод матрицы в консоль.
template<typename T>
void Matrix<T>::output(std::ofstream& fout)
{
    for (auto& row : matrix)
    {
        for (T elem : row)
            fout << elem << "\t";
        fout << "\n";
    }
    fout << std::setprecision(7) << std::endl;
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


// Сумма матриц
template<typename T>
Matrix<T> Matrix<T>::operator+(Matrix<T>& B)
{
    Matrix<T> sum(this->matrix.size(), this->matrix[0].size());
    for(size_t i = 0; i < this->matrix.size(); i++)
        for(size_t j = 0; j < this->matrix[0].size(); j++)
            sum.matrix[i][j] = this->matrix[i][j] + B.matrix[i][j];

    return sum;
}



template<typename T>
Matrix<T> Matrix<T>::operator-(Matrix<T>& B)
{
    Matrix<T> diff(this->matrix.size(), this->matrix[0].size());
    for(size_t i = 0; i < this->matrix.size(); i++)
        for(size_t j = 0; j < this->matrix[0].size(); j++)
            diff.matrix[i][j] = this->matrix[i][j] - B.matrix[i][j];

    return diff;
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


template<typename T>
T Matrix<T>::scalarMul(Matrix<T>& first, Matrix<T>& second)
{
    T result = 0;
    for (size_t i = 0; i < first.matrix.size(); i++)
    {
        result += (first[i][0] * second[i][0]);
    }
    return result;
}


template<typename T>
T Matrix<T>::vectorLength(Matrix<T>& vector)
{
    T sum = 0;
    for (size_t i = 0; i < vector.matrix.size(); i++)
        sum += std::pow(vector.matrix[i][0], 2);

    return std::sqrt(sum);
}


// Умножение строки на скаляр.
template<typename T>
void Matrix<T>::mulRowAndScalar(size_t row_id, T scalar)
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
        //std::cout << "Matrix hasn`t got below zero row or column!";
        exit(1);
    }
    
    auto& source = matrix[src];
    auto& destination = matrix[dest];
    for (size_t i = 0; i < source.size(); i++)
        destination[i] += source[i];
}


//
template<typename T>
Matrix<T> Matrix<T>::Unv()
{
 Matrix<T> Unsv(this->matrix.size(),this->matrix.size());
 for(size_t i = 0; i < this->matrix.size();i++)
 {
   Matrix<T> A = *this;
   Matrix<T> e(this->matrix.size(),1);
   e.matrix[i][0] = 1;
   Matrix<T> _A = makeExtendedMatrix(A,e);
   _A.forwardGaussStep();
   Matrix<T> usv = _A.backwardGaussStep();
   for(size_t j = 0 ; j < this->matrix.size(); j++)
   {
    Unsv.matrix[j][i] = usv.matrix[j][0];
   }
 }
 return Unsv;
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
                    //std::cout << "Matrix determinate is zero!\n\n";
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


// Возвращает расширенную матрицу
template<typename T>
Matrix<T> makeExtendedMatrix(Matrix<T>& A, Matrix<T>& b)
{
    Matrix<T> A_(A.matrix.size(),  A.matrix[0].size() + 1);
    for (size_t i = 0; i < A.matrix.size(); i++)
    {
        for (size_t j = 0; j < A.matrix[0].size() + 1; j++)
        {
            (j == A.matrix[0].size()) ?
            A_.matrix[i][j] = b.matrix[i][0] :
            A_.matrix[i][j] = A.matrix[i][j];
        }
    }
    return A_;
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
        //std::cout << "index of column or first_row_id out of range!\n\n";
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
    if (this->matrix.size() == other.matrix.size() &&
        this->matrix[0].size() == other.matrix[0].size())
    {
        for (size_t i = 0; i < this->matrix.size(); i++)
        {
            this->matrix[i].clear();
            this->matrix[i] = other.matrix[i];
        }
    }

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
        //std::cout << "dimentions of A and B are different!\n\n";
        exit(1);
    }
}


// Создает единичную матрицу размера nxn
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


// Возвращает евклидову норму вектора
template<typename K>
K eqlidNorm(Matrix<K> a)
{
    K sum = 0;
    for(int i = 0; i < a.matrix.size(); i++)
        sum += pow(a.matrix[i][0], 2);

    return sqrt(sum);
}


// Первая норма матрицы
template<typename K>
K octahedralNorm(Matrix<K> A)
{
    K max_sum = 0;
    for(size_t j = 0; j < A.matrix[0].size(); j++)
    {  
        K sum = 0;
        for(size_t i = 0; i < A.matrix.size(); i++)
            sum += std::abs(A.matrix[i][j]);
        if (sum > max_sum)
            max_sum = sum;
    }
    return max_sum;
}


// Вторая норма матрицы
template<typename K>
K cubicNorm(Matrix<K> A)
{
    K max_sum = 0;
    for(size_t i = 0; i < A.matrix.size(); i++)
    {  
        K sum = 0;
        for(size_t j = 0; j < A.matrix[0].size(); j++)
            sum += std::abs(A.matrix[i][j]);

        if (sum > max_sum)
            max_sum = sum;
    }
    return max_sum;
}


// Быстрое умножение на матрицу поворота
template<typename K>
void RotationForHess(Matrix<K>& A, size_t first_id, size_t second_id, Matrix<K>& B, bool transpose)
{
    auto first_row = A.matrix[first_id];
    auto second_row = A.matrix[second_id];
    K c1 = B.matrix[first_id][first_id - 1] / std::sqrt(std::pow(B.matrix[first_id][first_id - 1], 2) + std::pow(B.matrix[second_id][first_id - 1], 2)); 
    K c2 = B.matrix[second_id][first_id - 1] / std::sqrt(std::pow(B.matrix[first_id][first_id - 1], 2) + std::pow(B.matrix[second_id][first_id - 1], 2));
    for(size_t i = 0; i < A.matrix[0].size(); i++)
    {
        if (!transpose)
        {
            A.matrix[first_id][i] = c1 * first_row[i] + c2 * second_row[i];
            A.matrix[second_id][i] = (-c2 * first_row[i]) + c1 * second_row[i];
        }
        else
        {
            A.matrix[first_id][i] = c1 * first_row[i] + (-c2 * second_row[i]);
            A.matrix[second_id][i] = c2 * first_row[i] + c1 * second_row[i];
        }
    }
}


template<typename K>
void Rotation(Matrix<K>& A, size_t first_id, size_t second_id, Matrix<K>& B, bool transpose)
{
    auto first_row = A.matrix[first_id];
    auto second_row = A.matrix[second_id];
    K c1 = B.matrix[first_id][first_id] / std::sqrt(std::pow(B.matrix[first_id][first_id], 2) + std::pow(B.matrix[second_id][first_id], 2)); 
    K c2 = B.matrix[second_id][first_id] / std::sqrt(std::pow(B.matrix[first_id][first_id], 2) + std::pow(B.matrix[second_id][first_id], 2));
    for(size_t i = 0; i < A.matrix[0].size(); i++)
    {
        if (!transpose)
        {
            A.matrix[first_id][i] = c1 * first_row[i] + c2 * second_row[i];
            A.matrix[second_id][i] = (-c2 * first_row[i]) + c1 * second_row[i];
        }
        else
        {
            A.matrix[first_id][i] = c1 * first_row[i] + (-c2 * second_row[i]);
            A.matrix[second_id][i] = c2 * first_row[i] + c1 * second_row[i];
        }
    }
}


//QR разложение
template<typename K>
bool QRforSLAU(Matrix<K>& _A, Matrix<K>& Q, Matrix<K>& R, Matrix<K>& Uns)
{ 
    static double epsilon = 1e-10;
    auto Q1 = makeE<K>(_A.matrix.size());

    for(size_t i = 0; i < _A.matrix.size() - 1; i++)
        for(size_t j = i + 1; j < _A.matrix.size(); j++)
        {
            if(_A.matrix[j][i])
            {
                Rotation(Q1, i, j, _A);
                Rotation(_A, i, j, _A);
            }
            else
                continue;
        }

    for(size_t i = 0; i < _A.matrix.size(); i++)
        if(std::abs(_A.matrix[i][i]) < epsilon)
            return true;

    auto Ap =_A.divideExtendedMatrix();

    R = Ap.first;
    Q1.transpoce();
    Q = Q1;
    Matrix<K>  x = Ap.second;
    Matrix<K> R_Ext = makeExtendedMatrix(R,x);
    Matrix<K> Uns_ = R_Ext.backwardGaussStep();
    Uns = Uns_;

    return false;
}


// Возвращает число обусловленности
template<typename T>
std::pair<T,T> Matrix<T>::Cond()
{
    auto ThisUnv = this->Unv();
    T n11 = octahedralNorm<T>(*this);
    T n12 = octahedralNorm<T>(ThisUnv);
    T n21 = cubicNorm<T>(*this);
    T n22 = cubicNorm<T>(ThisUnv);
    return std::pair(n12 * n11, n21 * n22);
}


// Вносимая погрешность 
template<typename T>
Matrix<T> pgr(T delta, Matrix<T> b)
{
    Matrix<T> usv(b.matrix.size(), 1);
    for(size_t i =0; i < b.matrix.size(); i++)
        usv.matrix[i][0] = b.matrix[i][0] + delta;
    return usv;
}



////////////////////////////////////////////////////

//         //\\    ///////    |||||||||
//        //  \\   //    //   ||    ///
//       ////\\\\  ////////       ///
//      //      \\ //    //     ///
/////////        \\///////    //////////

////////////////////////////////////////////////////



// Методы, необходимые для работы Второй лабораторной работы


// Умножение матрицы на число
template<typename T>
void Matrix<T>::mulMatrixAndScalar(T scalar)
{
    for (size_t i = 0; i < this->matrix.size(); i++)
        this->mulRowAndScalar(i, scalar);
}


// Метод преобразования матрицы системы к виду, необходимому для итерационных методов
template<typename T>
std::pair<Matrix<T>, Matrix<T>> Matrix<T>::makeMatrixGoodForIter(IterativeMethod method, std::pair<Matrix<T>, Matrix<T>> divided)
{
    /* 
    В зависимости от метода возвращает пару - матрицу С и вектор-столбец
    правой части умноженный на итерационный параметр - y 
    */

    // divided --> std::pair (first = A, second = b, СЛАУ Ax = b)

    if (method == IterativeMethod::SIMPLE_ITER)
    {
        auto E = makeE<T>(this->matrix.size());
        divided.first.mulMatrixAndScalar(ITER_PARAM);
        divided.second.mulMatrixAndScalar(ITER_PARAM);
        Matrix<T> C = E - divided.first;
        Matrix<T> y = divided.second;

        return std::make_pair(C, y);
    }
    else if (method == IterativeMethod::JACOBI || method == IterativeMethod::SEIDEL)
    {
        Matrix<T> C = divided.first;
        Matrix<T> y = divided.second;
        switch (method)
        {
        case IterativeMethod::JACOBI:
        {
            // c(i,j) = -(a(i, j) / a(i, i)) если (i != j).
            for (size_t i = 0; i < this->matrix.size(); i++)
                for (size_t j = 0; j < this->matrix.size(); j++)
                {
                    if (i == j)
                    {
                        y.matrix[i][0] = divided.second.matrix[i][0] / divided.first.matrix[i][j];
                        C.matrix[i][j] = 0;
                    }
                    else
                        C.matrix[i][j] = -(divided.first.matrix[i][j] / divided.first.matrix[i][i]);
                }
            break;
        }
        case IterativeMethod::SEIDEL:
        {
            for (size_t i = 0; i < MATRIX_ORD; i++)
            {
                C.matrix[i][i] = 0;
                y.matrix[i][0] = divided.second.matrix[i][0] / divided.first.matrix[i][i];
                if (i == 0)
                    C.matrix[i][i + 1] = -(divided.first.matrix[i][i + 1] / divided.first.matrix[i][i]);
                else if (i == MATRIX_ORD - 1)
                    C.matrix[i][i - 1] = -(divided.first.matrix[i][i - 1] / divided.first.matrix[i][i]);
                else
                {
                    C.matrix[i][i + 1] = -(divided.first.matrix[i][i + 1] / divided.first.matrix[i][i]);
                    C.matrix[i][i - 1] = -(divided.first.matrix[i][i - 1] / divided.first.matrix[i][i]);
                }
            }
        }
        default:
            break;
        }
        return std::make_pair(C, y);
    }
    else if (method == IterativeMethod::RELAXATION)
    {
        Matrix<T> C = divided.first;
        Matrix<T> y = divided.second;

        for (size_t i = 0; i < MATRIX_ORD; i++)
        {
            C.matrix[i][i] = 0;
            y.matrix[i][0] = ((divided.second.matrix[i][0] * RELAX_PARAM)/ divided.first.matrix[i][i]);
            if (i == 0)
                C.matrix[i][i + 1] = -((divided.first.matrix[i][i + 1] * RELAX_PARAM) / divided.first.matrix[i][i]);
            else if (i == MATRIX_ORD - 1)
                C.matrix[i][i - 1] = -((divided.first.matrix[i][i - 1] * RELAX_PARAM) / divided.first.matrix[i][i]);
            else
            {
                C.matrix[i][i + 1] = -((divided.first.matrix[i][i + 1] * RELAX_PARAM) / divided.first.matrix[i][i]);
                C.matrix[i][i - 1] = -((divided.first.matrix[i][i - 1] * RELAX_PARAM) / divided.first.matrix[i][i]);
            }
        }
        return std::make_pair(C, y);
    }
    else
    {
        std::cout << "Unknown method!\n";
        exit(1);
    }

    return std::make_pair(divided.first, divided.second);
}


template<typename T>
Matrix<T> Matrix<T>::simpleIterationsMethod(Matrix<T> initial_approx)
{
    size_t iter_counter = 0;
    
    // Получаем вид С = (E - ITER_PARAM * A) & y = iter_param * b
    auto pair = this->makeMatrixGoodForIter(IterativeMethod::SIMPLE_ITER, this->divideExtendedMatrix());

    // Считаем норму + проверка на сходимость метода (|| C || < 1)
    auto norm_C = cubicNorm(pair.first);
    if (norm_C >= 1)
    {
        std::cout << "|| C || >= 1 !!! (|| C || = " << norm_C << ")\n";
        exit(1);
    }
    
    Matrix<T> solution(this->matrix.size(), 1);
    Matrix<T> temp(this->matrix.size(), 1);
    std::ofstream SIM("SIM-iters.txt");
    // Выполняем итерацию. new_solution = x ^ (n + 1) ; temp = x ^ (n).
    do
    {
        iter_counter++;
        // x ^ (n + 1) = C * x ^ (n) + ITER_PARAM * b 
        Matrix<T> new_solution = multiply<T>(pair.first, initial_approx) + pair.second;
        temp = initial_approx;
        solution = new_solution;
        initial_approx = solution;
        solution.output(SIM);
        SIM << "\n";
    } while (cubicNorm(solution - temp) >= (((1 - norm_C) / norm_C) * epsilon));
    SIM.close();
    // Есть возможносоть посмотреть количество итераций метода
    // Чем меньше итерационный параметр -> тем больше итераций.
    std::cout << "Number of iterations in SIM: " << iter_counter << "\n";

    return solution;
}


template<typename T>
Matrix<T> Matrix<T>::methodJacobi(Matrix<T> initial_approx)
{
    size_t iter_counter = 0;
    
    // Получаем вид С = -D ^ (-1) * (U + L) & y = D ^ (-1) * b
    auto pair = this->makeMatrixGoodForIter(IterativeMethod::JACOBI, this->divideExtendedMatrix());

    auto norm_C = cubicNorm(pair.first);
    if (norm_C >= 1)
    {
        std::cout << "|| C || >= 1 !!! (|| C || = " << norm_C << ")\n";
        exit(1);
    }

    Matrix<T> solution(this->matrix.size(), 1);
    Matrix<T> temp(this->matrix.size(), 1);
    // Выполняем итерацию. new_solution = x ^ (n + 1) ; temp = x ^ (n).
    do
    {
        iter_counter++;
        // x ^ (n + 1) = C * x ^ (n) + y
        Matrix<T> new_solution = multiply<T>(pair.first, initial_approx) + pair.second;
        temp = initial_approx;
        solution = new_solution;
        initial_approx = solution;
    } while (cubicNorm(solution - temp) >= (((1 - norm_C) / norm_C) * epsilon));
    
    // Есть возможносоть посмотреть количество итераций метода
    // Чем меньше итерационный параметр -> тем больше итераций.
    std::cout << "Number of iterations in Jacobi: " << iter_counter << "\n";

    return solution;
}


// Дружественные функции



template<typename V>
std::pair<Matrix<V>, double>  relaxationMethod(Matrix<V> initial_approx, Matrix<V> U, Matrix<V> D, Matrix<V> L, Matrix<V> b, IterativeMethod method)
{
    static double eps = 10e-6;
    size_t iter_counter = 0;
    // U - "Верхняя треугольная", L - "нижняя". Формируем С.
    Matrix<V> C(MATRIX_ORD, MATRIX_ORD);
    for (size_t i = 0; i < MATRIX_ORD; i++)
    {
        C.matrix[i][i] = D[i][0];
        if (i != MATRIX_ORD - 1)
        {
            C.matrix[i][i + 1] = U.matrix[i][0];
            C.matrix[i + 1][i] = L.matrix[i][0];
        }
    }
    std::ofstream y("vectorY");
    C.output(y);
    y << "\n\n";
    auto pair = (method == IterativeMethod::SEIDEL) ?
                (C.makeMatrixGoodForIter(IterativeMethod::SEIDEL, std::make_pair(C, b))) :
                (C.makeMatrixGoodForIter(IterativeMethod::RELAXATION, std::make_pair(C, b)));

    auto norm_C = cubicNorm(pair.first);
    if (norm_C >= 1)
    {
        std::cout << "|| C || >= 1 !!! (|| C || = " << norm_C << ")\n";
        exit(1);
    }
    // std::cout << "|| C || = " << norm_C << "\n";
    Matrix<V> solution(MATRIX_ORD, 1);
    Matrix<V> temp(MATRIX_ORD, 1);


    y << "Vector Y: \n";
    pair.second.output(y);
    y << "\n\nMatrix C: \n";
    pair.first.output(y);

    double param = (method == IterativeMethod::SEIDEL) ? 1 : RELAX_PARAM;
    // Выполняем итерацию.
    do
    {
        iter_counter++;
        Matrix<V> new_solution(MATRIX_ORD, 1);
        for (size_t i = 0; i < MATRIX_ORD; i++)
        {
            if (i == 0)
                new_solution[i][0] = (initial_approx[i][0] * (1 - param) + initial_approx[i + 1][0] * pair.first.matrix[i][i + 1] + pair.second.matrix[i][0]);
            else if (i == MATRIX_ORD - 1)
                new_solution[i][0] = (initial_approx[i][0] * (1 - param) + new_solution[i - 1][0] * pair.first.matrix[i][i - 1] + pair.second.matrix[i][0]);
            else
            {
                new_solution[i][0] = ((initial_approx[i][0] * (1 - param)) + (initial_approx[i + 1][0] * pair.first.matrix[i][i + 1]) +
                            (new_solution[i - 1][0] * pair.first.matrix[i][i - 1]) + pair.second.matrix[i][0]);
            }
        }
        temp = initial_approx;
        solution = new_solution;
        initial_approx = solution;
    } while (cubicNorm(solution - temp) / cubicNorm(temp) >= epsilon);
    
    // Есть возможносоть посмотреть количество итераций метода
    // Чем меньше итерационный параметр -> тем больше итераций.
    if (method == IterativeMethod::SEIDEL)
        std::cout << "Number of iterations in Seidel: " << iter_counter << "\n";
    else
        std::cout << "Number of iterations in Relax: " << iter_counter << "\n";
    
    
    auto residual = multiply(C, solution);
    auto norm = cubicNorm(b - residual);
    return std::make_pair(solution, norm);
}



////////////////////////////////////////////////////

//         //\\    ///////     |||||||
//        //  \\   //    //   ||    |||
//       ////\\\\  ////////      |||||
//      //      \\ //    //   ||    |||
/////////        \\///////     |||||||

////////////////////////////////////////////////////


template<typename T>
Matrix<T> Matrix<T>::reduction(size_t reduction_deg)
{
    if (this->matrix.size() <= reduction_deg - 1)
    {
        std::cout << "Невозможно предъявить сужение для данной матрицы!!";
        exit(1);
    }
    auto new_size = this->matrix.size() - reduction_deg;
    Matrix<T> answer(new_size, new_size);
    for (size_t i = 0; i < new_size; i++)
        for (size_t j = 0; j < new_size; j++)
            answer[i][j] = this->matrix[i][j];

    return answer;
}


template<typename T>
std::pair<Matrix<T>, Matrix<T>> Matrix<T>::QR(Matrix<T>& A)
{
    auto Q = makeE<T>(A.matrix.size());
    auto& temp = A;
    for (size_t i = 0; i < A.matrix.size() - 1; i++)
        for (size_t j = i + 1; j < A.matrix.size(); j++)
        {
            if (A.matrix[j][i])
            {
                Rotation(Q, i, j, A, false);
                Rotation(A, i, j, A, false);
            }
            else
                continue;
        }

    Q.transpoce();

    return std::make_pair(A, Q);
}


template<typename V>
void exitQRCriteria(Matrix<V>& A, bool& continue_flag, size_t const size, ExitQRCriteria condition)
{
    switch (condition)
    {
    case ExitQRCriteria::SIMPLE:
    {
        auto last_row = A.operator[](size - 1);

        // Флаг - true, если точность
        // В последней строке еще не достигнута
        for (size_t i = 0; i < size - 1; i++)
            if (std::abs(last_row[i]) > epsilon)
            {
                continue_flag = true;
                return;
            }
        return;
    }
    case ExitQRCriteria::MODULE:
    {
        V sum = 0;
        for (size_t i = 0; i < size; i++)
        {
            for (size_t j = 0; j < i; j++)
            {
                sum += std::abs(A.matrix[i][j]);
            }
        }
        if (sum < epsilon)
            continue_flag = true;
        
        return;
    }
    case ExitQRCriteria::SQUARE_MODULE:
    {
        V sum = 0;
        for (size_t i = 0; i < size; i++)
        {
            for (size_t j = 0; j < i; j++)
            {
                sum += std::pow(std::abs(A.matrix[i][j]), 2);
            }
        }
        if (sum < epsilon)
            continue_flag = true;
        
        return;
    }
    default:
    {
        std::cout << "\nУсловие выхода из итерационного метода QR - некорректно!\n";
        exit(1);
    }
    }
}


template<typename T>
Matrix<T> Matrix<T>::iterativeQR(Matrix<T> A)
{
    // static double eps = 10e-6;
    size_t size = this->matrix.size();
    Matrix<T> answer(size, 1);
    while (size != 1)
    {
        // auto sigma = A[size - 1][size - 1];
        auto sigma = 0;
        auto E = makeE<T>(size);
        E.mulMatrixAndScalar(sigma);
        A = A - E;
        auto QR = A.QR(A);
        auto temp = multiply(QR.first, QR.second) + E;

        bool continue_flag = false;
        
        exitQRCriteria(temp, continue_flag, size, ExitQRCriteria::SIMPLE);

        if (continue_flag)
        {
            A = temp; // A = A^k+1
            continue;
        }
        
        answer.matrix[size - 1][0] = temp[size - 1][size - 1];
        auto newA = temp.reduction(1);

        // Берем угловой минор на порядок меньше
        for (size_t iter = 0; iter < size; iter++)
        {
            if (iter < (size - 1))
                A[iter].resize(size - 1);
            else
            {
                A[iter].clear();
                A.matrix.resize(size - 1);
            }
        }
        A = newA;
        --size;

        if (size == 1)
            answer.matrix[size - 1][0] = newA[0][0];
    }
    return answer;
}


template<typename T>
void Matrix<T>::makeHessenbergForm(Matrix<T>& A)
{
    auto& temp = A;
    for (size_t i = 1; i < A.matrix.size() - 1; i++)
        for (size_t j = i + 1; j < A.matrix.size(); j++)
        {
            A.transpoce();
            if (A.matrix[j][i - 1])
                RotationForHess(A, i, j, A, false);
            else
                continue;
            A.transpoce();
            if (A.matrix[j][i - 1])
                RotationForHess(A, i, j, A, false);
            else
                continue;
        }
}


template<typename T>
std::pair<Matrix<T>, T> Matrix<T>::reverseIterationsMethod(bool rayleigh, Matrix<T>& A, Matrix<T> initial_approx, T lambda, size_t* iter)
{
    auto E = makeE<T>(A.matrix.size());
    if (rayleigh)
    {
        auto image = multiply(A, initial_approx);                                           // image       = Ax
        lambda = (scalarMul(image, initial_approx) / eqlidNorm(initial_approx));            // lambda      = (Ax, x) / (x, x)
        E.mulMatrixAndScalar(lambda);
        auto temp_matrix = A - E;                                                           // temp_matrix = A - lambda * E

        // Создается расширенная матрица системы для метода Гаусса
        auto linear_sys = makeExtendedMatrix<T>(temp_matrix, initial_approx);
        if (!linear_sys.forwardGaussStep())
        {
            auto new_approx = linear_sys.backwardGaussStep();                      // y^k+1
            new_approx.mulMatrixAndScalar(1 / eqlidNorm(new_approx));              // x^k+1 = y^k+1 / || y^k+1 ||
            while ((1.0 - std::abs(scalarMul(new_approx, initial_approx) / (vectorLength(new_approx) * vectorLength(initial_approx)))) > epsilon)
            {
                if (iter)
                    (*iter)++;
                for (int i = 0; i < E.matrix.size(); i++)
                    E.matrix[i][i] = 1;
                
                image = multiply(A, new_approx);                                   // image  = A*x^k
                lambda = (scalarMul(image, new_approx) / eqlidNorm(new_approx));   // lambda = (A*x^k, x^k) / (x^k, x^k)
                E.mulMatrixAndScalar(lambda);
                temp_matrix = A - E;
                initial_approx = new_approx;
                linear_sys = makeExtendedMatrix<T>(temp_matrix, new_approx);
                if (!linear_sys.forwardGaussStep())
                {
                    new_approx = linear_sys.backwardGaussStep();
                    new_approx.mulMatrixAndScalar(1 / eqlidNorm(new_approx));
                }
            }
            initial_approx = new_approx;
        }
        return std::make_pair(initial_approx, lambda);
    }
    else
    {
        E.mulMatrixAndScalar(lambda);
        auto matrix = A - E;
        auto linear_sys = makeExtendedMatrix<T>(matrix, initial_approx);
        if (!linear_sys.forwardGaussStep())
        {
            auto new_approx = linear_sys.backwardGaussStep();
            new_approx.mulMatrixAndScalar(1 / eqlidNorm(new_approx));
            while ((1 - std::abs(scalarMul(new_approx, initial_approx) / (vectorLength(new_approx) * vectorLength(initial_approx)))) > epsilon)
            {
                if (iter)
                    (*iter)++;
                initial_approx = new_approx;
                linear_sys = makeExtendedMatrix<T>(matrix, new_approx);
                if (!linear_sys.forwardGaussStep())
                {
                    new_approx = linear_sys.backwardGaussStep();
                    new_approx.mulMatrixAndScalar(1 / eqlidNorm(new_approx));
                }
            }
            initial_approx = new_approx;
        }
        return std::make_pair(initial_approx, lambda);
    }
}
