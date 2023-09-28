#include "../matrix/Matrix.hpp"
#define T double

static constexpr int TESTS_COUNT = 5;
static constexpr T INIT_PARAM = 0.01;


int main()
{
    std::ofstream fout("output_lab2.txt");
    Matrix<T> Test(4, 5, "/home/george/Desktop/Projects/C++/CalculationMethods/lab2/tests/test1.txt");
    Matrix<T> approx(4, 1, "/home/george/Desktop/Projects/C++/CalculationMethods/lab2/tests/first_approx.txt");
    auto E = makeE<T>(Test.getMatrix().size());
    auto solution = Test.simpleIterationsMethod(approx, E, INIT_PARAM);
    solution.output(fout);
}