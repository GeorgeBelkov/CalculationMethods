#include "../include/Matrix.hpp"


using type_t = double;

static constexpr size_t PHI_NODES_NUM = 50;
static constexpr size_t ETTA_NODES_NUM = 25;


void makeGrid(std::vector<double>& nodes, std::pair<double, double> segment, size_t nodes_number)
{
    auto dist = std::abs(segment.first - segment.second) / nodes_number;
    int iter = 0;
    for (auto& node : nodes)
    {
        node = segment.first + iter * dist;
        iter++;
    }
}


void sphToDec(Matrix<type_t>& vector, double etta, double phi, double r = 1)
{
    if (vector.getMatrix().size() != 3)
    {
        std::cout << "Для перевода из сферических в декартовы нужен вектор размерности 3!";
        exit(1);
    }
    vector.getMatrix()[0][0] = r * std::sin(etta) * std::cos(phi);
    vector.getMatrix()[1][0] = r * std::sin(etta) * std::sin(phi);
    vector.getMatrix()[2][0] = r * std::sin(etta);
}



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


    // Дополнительное задание
    std::ofstream experiment("experement.txt");
    Matrix<type_t> approx(3, 1);
    Matrix<type_t> matrix(3, 3, "/home/george/Desktop/Projects/C++/CalculationMethods/lab4/experiment/matrix.txt");
    Matrix<type_t> lambda(3, 1, "/home/george/Desktop/Projects/C++/CalculationMethods/lab4/experiment/lambdas.txt");
    Matrix<type_t> e1(3, 1, "/home/george/Desktop/Projects/C++/CalculationMethods/lab4/experiment/e1.txt");
    Matrix<type_t> e2(3, 1, "/home/george/Desktop/Projects/C++/CalculationMethods/lab4/experiment/e2.txt");
    Matrix<type_t> e3(3, 1, "/home/george/Desktop/Projects/C++/CalculationMethods/lab4/experiment/e3.txt");

    std::vector<double> phi_nodes(PHI_NODES_NUM), etta_nodes(ETTA_NODES_NUM);
    auto phi_segment = std::make_pair<double, double>(0, 2 * M_PI);
    auto etta_segment = std::make_pair<double, double>(0, M_PI);

    makeGrid(phi_nodes, phi_segment, PHI_NODES_NUM);
    makeGrid(etta_nodes, etta_segment, ETTA_NODES_NUM);
    auto grid = std::make_pair(phi_nodes, etta_nodes);


    size_t ctr = 1;


    for (size_t k = 0; k < lambda.getMatrix().size(); k++)
    {
        for (size_t i = 0; i < PHI_NODES_NUM; i++)
        {
            for (size_t j = 0; j < ETTA_NODES_NUM; j++)
            {
                sphToDec(approx, etta_nodes[j], phi_nodes[i]);
                size_t* iter_counter = &ctr;
                auto answ = matrix.reverseIterationsMethod(true, matrix, approx, lambda.getMatrix()[k][0], iter_counter);
                std::vector<double> norms = { eqlidNorm(e1 - answ.first), eqlidNorm(e2 - answ.first), eqlidNorm(e3 - answ.first) };
                auto nearest = std::min_element(norms.begin(), norms.end());
                for (int iter = 0; iter < norms.size(); iter++)
                {
                    if (norms[iter] == *nearest)
                    {
                        experiment << (iter + 1) << " " << *iter_counter << "\n";
                        break;
                    }
                }
                ctr = 1;
            }
        }
    }
}