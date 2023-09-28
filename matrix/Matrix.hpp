#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <cmath>
#include <fstream>


constexpr static double epsilon = 10e-4;


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
    Matrix<T> operator+(Matrix<T>& B);
    Matrix<T> operator-(Matrix<T>& B);

    // Геттеры и Сеттеры
    std::vector<std::vector<T>>& getMatrix() { return matrix; }

    // Вспомогательные методы класса
    void output(std::ofstream& fout);
    std::vector<T> getColumn(size_t index, size_t first_row_id = 0);
    T scalarMul(std::vector<T>& first, std::vector<T>& second);
    void transpoce();
    void mulRowAndScalar(size_t row_id, T scalar);
    void mulMatrixAndScalar(T scalar);
    void addRow(size_t src, size_t dest);
    std::pair<Matrix<T>, Matrix<T>> makeMatrixGoodForIter(T iter_param, Matrix<T>& E);
    std::pair<Matrix<T>, Matrix<T>> divideExtendedMatrix();
    Matrix<T> Unv();
    std::pair<T,T> Cond();

    // Методы по заданию
    bool forwardGaussStep();
    Matrix<T> backwardGaussStep();
    Matrix<T> simpleIterationsMethod(Matrix<T> initial_approx, Matrix<T>& E, T iter_param);

    // Дружественные фунции

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
    friend bool QR(Matrix<V>& A,Matrix<V>& Q, Matrix<V>& R,Matrix<V>& Uns);

    template<typename K>
    friend void Rotation(Matrix<K>& A, size_t num1, size_t num2,Matrix<K>& B);

    template<typename V>
    friend Matrix<V> pgr(V delta,Matrix<V> b);
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
        //std::cout << "index of row out of range!\n\n";
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
    for(int i = 0;i < a.matrix.size();i++)
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
void Rotation(Matrix<K>& A, size_t num1, size_t num2,Matrix<K>& B)
{
    auto str1 = A.matrix[num1];
    auto str2 = A.matrix[num2];
    K c1 = B.matrix[num1][num1]/sqrt(pow(B.matrix[num1][num1],2) + pow(B.matrix[num2][num1],2)); 
    K c2 = B.matrix[num2][num1]/sqrt(pow(B.matrix[num1][num1],2) + pow(B.matrix[num2][num1],2));
    for(size_t i = 0; i < A.matrix[0].size(); i++)
    {
        A.matrix[num1][i] = c1*str1[i] +c2*str2[i];
        A.matrix[num2][i] = -c2*str1[i] +c1*str2[i];
    }
}


//QR разложение
template<typename K>
bool QR(Matrix<K>& _A, Matrix<K>& Q, Matrix<K>& R, Matrix<K>& Uns)
{ 
    static double epsilon = 1e-10;
    auto Q1 = makeE<K>(_A.matrix.size());

    for(size_t i = 0; i < _A.matrix.size() - 1; i++)
        for(size_t j = i + 1; j < _A.matrix.size(); j++)
        {
            if(_A.matrix[i][j])
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
std::pair<Matrix<T>, Matrix<T>> Matrix<T>::makeMatrixGoodForIter(T iter_param, Matrix<T>& E)
{
    /* Возвращает пару - матрицу С и вектор-столбец
    правой части умноженный на итерационный параметр - y */

    auto divided = this->divideExtendedMatrix();
    divided.first.mulMatrixAndScalar(iter_param);
    divided.second.mulMatrixAndScalar(iter_param);
    Matrix<T> y = divided.second;
    Matrix<T> C = E - divided.first;
    return std::make_pair(C, y);
}


template<typename T>
Matrix<T> Matrix<T>::simpleIterationsMethod(Matrix<T> initial_approx, Matrix<T>& E, T iter_param)
{
    auto pair = this->makeMatrixGoodForIter(iter_param, E);
    auto norm_C = cubicNorm(pair.first);
    if (norm_C >= 1)
    {
        std::cout << "|| C || >= 1 !!! (|| C || = " << norm_C << ")";
        exit(1);
    }
    
    Matrix<T> solution(this->matrix.size(), 1);
    Matrix<T> temp(this->matrix.size(), 1);
    do
    {
        Matrix<T> new_solution = multiply<T>(pair.first, initial_approx) + pair.second;
        temp = initial_approx;
        solution = new_solution;
        initial_approx = solution;
    } while (cubicNorm(solution - temp) >= (((1 - norm_C) / norm_C) * epsilon));
    return solution;
}