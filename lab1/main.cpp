#include "Matrix.hpp"


int main()
{
    Matrix<double> A(3, 3, "/home/ysan/numer/CalculationMethods/lab1/matrixA.txt");
    Matrix<double> C(3, 3, "/home/ysan/numer/CalculationMethods/lab1/matrixC.txt");
    A.output();
    Matrix<double> B = A;
    B.output();
    B = C;
    B.output();
    B = B;
    B[0] = A[1];
    B.output();
    Matrix<double> D = multiply(A, B);
    D.output();
    D.transpoce();
    D.output();
    Matrix<double> E = makeE<double>(4);
    E.output();
    C.mulRowAndScalar(1, 3);
    C.output();
    C.addRow(1, 2);
    C.output();
}