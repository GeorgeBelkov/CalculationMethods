#include "Matrix.hpp"


int main()
{
    Matrix<double> A(3, 3, "/home/ysan/numer/CalculationMethods/lab1/matrixA.txt");
    //Matrix<double> C(3, 3, "/home/ysan/numer/CalculationMethods/lab1/matrixC.txt");
    
    Matrix<double> Q(3,3),R(3,3);
    QR(A,Q,R);
    Q.output();
    R.output();

    Matrix<double> Q_ = Q;
    Q_.transpoce();
    Matrix<double> E = multiply(Q_,Q);
    E.output();
    //Matrix<double> Temp = makeRot(0,1,A);
    //Matrix<double> M1 = multiply(Temp,A);
    //Matrix<double> Temp1 = makeRot(0,2,M1);
    //Matrix<double> M2 = multiply(Temp1,M1);
    //Temp.output();
    //M2.output();
}
    