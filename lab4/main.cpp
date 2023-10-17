#include "../include/Matrix.hpp"


using type_t = double;


int main()
{
    Matrix<type_t> A(4, 4, "/home/george/Desktop/Projects/C++/CalculationMethods/lab4/tests/test.txt");
    Matrix<type_t> x(4, 1, "/home/george/Desktop/Projects/C++/CalculationMethods/lab4/approx/init_approx.txt");
    Matrix<type_t> lambdas(4, 1, "/home/george/Desktop/Projects/C++/CalculationMethods/lab4/approx/lambda_approx.txt");

    std::ofstream fout("out.txt");


    auto iterativeA = A, iterativeHes = A;
    A.makeHessenbergForm(iterativeHes);
    auto Hes = iterativeHes;

    auto iterative_answer = A.iterativeQR(iterativeA);
    auto hessenberg_answer = Hes.iterativeQR(iterativeHes);

    fout << "Reverse Iterations Method (without ralay):\n\n";
    for (size_t i = 0; i < lambdas.getMatrix().size(); i++)
    {
        fout << "e" << i + 1 << ":\n";
        auto reverse_iterative = A.reverseIterationsMethod(false, A, x, lambdas.getMatrix()[i][0]);
        reverse_iterative.first.output(fout);

        fout << "lambda: " << reverse_iterative.second << "\n\n";
    }
    fout << "Reverse Iterations Method (with ralay):\n\n";
    for (size_t i = 0; i < lambdas.getMatrix().size(); i++)
    {
        fout << "e" << i + 1 << ":\n";
        auto reverse_iterative = A.reverseIterationsMethod(true, A, x, lambdas.getMatrix()[i][0]);
        reverse_iterative.first.output(fout);

        fout << "lambda: " << reverse_iterative.second << "\n\n";
    }

    fout << "IterativeQR (without Hessenberg form): \n\n";
    iterative_answer.output(fout);
    fout << "IterativeQR (with Hessenberg form): \n\n";
    hessenberg_answer.output(fout);
}