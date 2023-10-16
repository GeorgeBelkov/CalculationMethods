#include "../include/Matrix.hpp"


using type_t = double;


int main()
{
    Matrix<type_t> A(4, 4, "/home/george/Desktop/Projects/C++/CalculationMethods/lab4/tests/test.txt");
    auto iterativeA = A;
    auto answer = A.iterativeQR(iterativeA);

    std::ofstream fout("out.txt");
    answer.output(fout);
}