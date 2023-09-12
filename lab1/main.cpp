#include "Matrix.hpp"


static constexpr int TESTS_COUNT = 5;


int main()
{
    for (size_t i = 1; i <= TESTS_COUNT; i++)
    {
        std::cout << i << " TEST:\n\n";
        Matrix<double> Test(4, 5, "/home/george/Desktop/Projects/C++/CalculationMethods/lab1/tests/test" + std::to_string(i) + ".txt");
        // Matrix<double> Test(4, 5, "/home/george/Desktop/Projects/C++/CalculationMethods/lab1/tests/test2.txt");
        auto pair = Test.divideExtendedMatrix();
        pair.first.output();
        pair.second.output();
        if (!Test.forwardGaussStep())
        {
            Matrix<double> answer = Test.backwardGaussStep();
            answer.output();
        }
    }

    
}