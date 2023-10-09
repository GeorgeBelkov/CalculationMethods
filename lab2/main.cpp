#include "../matrix/Matrix.hpp"
#define T double


template<typename V>
void setVectors(Matrix<V>& U, Matrix<V>& D, Matrix<V>& L, Matrix<V>& b)
{
    for (size_t i = 0; i < MATRIX_ORD; i++)
    {
        D.getMatrix()[i][0] = 4;
        if (i == 0)
        {
            U.getMatrix()[i][0] = 1;
            b.getMatrix()[i][0] = 6;
        }
        else if (i == MATRIX_ORD - 1)
        {
            L.getMatrix()[i - 1][0] = 1;
            b.getMatrix()[i][0] = 6;
        }
        else
        {
            U.getMatrix()[i][0] = L.getMatrix()[i - 1][0] = 1;
            b.getMatrix()[i][0] = 10 - 2 * ((i + 1) % 2);
        }
    }
}


int main()
{
    std::ofstream fout("output_lab2.txt");

    std::array<IterativeMethod, 4> methods = {
        IterativeMethod::SIMPLE_ITER,
        IterativeMethod::JACOBI,
        IterativeMethod::SEIDEL,
        IterativeMethod::RELAXATION
    };

    Matrix<T> Test(4, 5, "/home/george/Desktop/Projects/C++/CalculationMethods/lab2/tests/test2.txt");
    Matrix<T> approx_for_test1(4, 1, "/home/george/Desktop/Projects/C++/CalculationMethods/lab2/tests/test1_approx.txt");
    Matrix<T> approx_for_test2(MATRIX_ORD, 1);

    Matrix<T> U(MATRIX_ORD - 1, 1), L(MATRIX_ORD - 1, 1), D(MATRIX_ORD, 1), b(MATRIX_ORD, 1);
    // Заполняем векторы U, D, L, b
    setVectors<T>(U, D, L, b);


    for (auto method : methods)
    {
        switch (method)
        {
        case IterativeMethod::SIMPLE_ITER:
        {
            fout << "ITERATIVE METHOD: simple Iterations Method\n\n";
            auto solution = Test.simpleIterationsMethod(approx_for_test1);
            fout << "SOLUTION:\n";
            solution.output(fout);
            fout << "\nNorm of the residual vector: ";
            auto pair = Test.divideExtendedMatrix();
            auto temp = multiply(pair.first, solution);
            fout << cubicNorm(pair.second - temp) << "\n\n";
            break;
        }
        case IterativeMethod::JACOBI:
        {
            fout << "ITERATIVE METHOD: Jacobi Method\n\n";
            auto solution = Test.methodJacobi(approx_for_test1);
            fout << "SOLUTION:\n";
            solution.output(fout);
            fout << "\nNorm of the residual vector: ";
            auto pair = Test.divideExtendedMatrix();
            auto temp = multiply(pair.first, solution);
            fout << cubicNorm(pair.second - temp) << "\n\n";
            break;
        }
        case IterativeMethod::SEIDEL:
        {
            fout << "ITERATIVE METHOD: Seidel Method\n\n";
            auto solution = relaxationMethod(approx_for_test2, U, D, L, b, method);
            fout << "SOLUTION:\n";
            solution.first.output(fout);
            fout << "\nNorm of the residual vector: " << solution.second << "\n\n";
            break;
        }
        case IterativeMethod::RELAXATION:
        {
            fout << "ITERATIVE METHOD: Relaxation Method\n\n";
            auto solution = relaxationMethod(approx_for_test2, U, D, L, b, method);
            fout << "SOLUTION:\n";
            solution.first.output(fout);
            fout << "\nNorm of the residual vector: " << solution.second << "\n\n";
            break;
        }
        default:
            break;
        }
        
    }
    


}