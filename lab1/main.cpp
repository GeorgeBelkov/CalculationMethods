#include "../include/Matrix.hpp"
#define T double

static constexpr int TESTS_COUNT = 5;


int main()
{
    std::ofstream fout("output.txt");
    for (size_t i = 1; i <= TESTS_COUNT; i++)
    {
        fout << i << " TEST:\n\n";
        Matrix<T> Test(4, 5, "/home/george/Desktop/Projects/C++/CalculationMethods/lab1/tests/test" + std::to_string(i) + ".txt");
        // Matrix<double> Test(4, 5, "/home/george/Desktop/Projects/C++/CalculationMethods/lab1/tests/test2.txt");
        Test.output(fout);
        Matrix<T> Q(4,4),R(4,4),usvQR(4,1);
        Matrix<T> TestG = Test;
        fout << "Gauss:\n";
        if (!TestG.forwardGaussStep())
          {
            Matrix<T>  usv = TestG.backwardGaussStep();
            fout << "Result:\n";
            usv.output(fout);
            std::pair<Matrix<T>, Matrix<T>>  pr = Test.divideExtendedMatrix();
            Matrix<T> pr1 = pr.first;
            Matrix<T> pr2 = pr.second;
            Matrix<T> b = multiply(pr1,usv);
            Matrix<T> dif = pr2 - b;
            fout << "|b - b1| = " << eqlidNorm(dif)<<std::endl;
            fout << "1) Cond A = " << pr1.Cond().first<<std::endl;
            fout << "2) Cond A = " << pr1.Cond().second<<std::endl;
            fout << "b := b + 0.01"<<std::endl;
            Matrix<T> _b = pgr(0.01,b);
            Matrix<T> new_sys = makeExtendedMatrix(pr1,_b);
            if(!new_sys.forwardGaussStep())
              {
                Matrix<T>  usv2 = new_sys.backwardGaussStep();
                fout << "Result2:\n";
                usv2.output(fout);
                fout << "Diff:\n";
                Matrix<T> diff2 = usv2 - usv;
                fout <<  eqlidNorm(diff2)<<std::endl;
              }


          }
        else
         fout<< "No solution\n";
        fout << "QR:\n";

        if (!QR(Test,Q,R,usvQR))
        { fout << "Q:\n";
         Q.output(fout);
         fout << "R:\n";
         R.output(fout);
         fout << "Result:\n";
         usvQR.output(fout);

         auto pr = Test.divideExtendedMatrix();
         Matrix<T> pr1 = pr.first;
         Matrix<T> pr2 = pr.second;
         Matrix<T> b = multiply(pr1,usvQR);
         Matrix<T> dif = pr2 - b;
         fout << "|b - b1| = " << eqlidNorm(dif)<<std::endl;
         fout << "b := b + 0.01"<<std::endl;
         Matrix<T> _b = pgr(0.01,b);
         Matrix<T> new_sys = makeExtendedMatrix(pr1,_b);
         Matrix<T> usvQR2(4,1);
         if(!QR(new_sys,Q,R,usvQR2))
              {
                fout << "Result2:\n";
                usvQR2.output(fout);
                fout << "Diff:\n";
                Matrix<T> diff2 = usvQR2 - usvQR;
                fout <<  eqlidNorm(diff2)<<std::endl;
              }
        }
        else
         fout << "No solution\n";
        
        
    } 

    
}