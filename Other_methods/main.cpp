#include <iostream>
#include "Matrix.h"
#include "SlauSolver.h"
#include "SLAU.h"
#include <chrono>

void SaveVec(const std::string &filename, int32_t precession)
{

}


int main()
{   

    /* Истенный вектор x */
    double x[4394];
    double e = -1;
    for(int i = 0; i < 4394; i++)
    {
        if(i % 25 == 0)
            e += 1.0;
        x[i] = e;
    }


    double x_true_norm = SLAUSolvers::IterationSolvers::Norma(x, 4394);

    SLAUSolvers::IterationSolvers::InitDateFile InitStruct;

    InitStruct.type = SLAUSolvers::MatrixType::SPARSE;
    InitStruct.di = "Test3/di.txt";
    InitStruct.ig = "Test3/ig.txt";
    InitStruct.jg = "Test4/jg.txt";
    InitStruct.ggl = "Test4/ggl.txt";
    InitStruct.ggu = "Test4/ggu.txt";
    InitStruct.kuslau = "Test4/kuslau.txt";
    InitStruct.f = "Test4/f.txt";
    
    
    SLAUSolvers::IterationSolvers::SLAU slau;

    SLAUSolvers::IterationSolvers::Load(InitStruct, slau, true);
    SLAUSolvers::IterationSolvers::MultA(slau.matrix, x, slau.f);

    SLAUSolvers::IterationSolvers::SetSolveMode(slau, SLAUSolvers::IterationSolvers::Solvemode::LOS_NOSYMETRIC_LU_FACT);
    
    auto begin = std::chrono::steady_clock::now();
    SLAUSolvers::IterationSolvers::SolveSLAU(slau, true);
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << elapsed_ms.count() << std::endl;


    SLAUSolvers::IterationSolvers::SaveX(slau, "Test4/x.txt", "\n");
    
    
    
     SLAUSolvers::IterationSolvers::ActionVec(1.0, x, -1.0, slau.x, x, slau.N);
     double r = SLAUSolvers::IterationSolvers::Norma(x, slau.N)/x_true_norm;
    
  

    /* Преобразуем внутренний формат матрицы решателя SLAU в формат для конвертации */
    /*Matrix::SparseMatrix SparseMatrix(slau.N, slau.size);
    int N = slau.N;
    int size = slau.size;
    std::copy(slau.matrix.di, slau.matrix.di + N, SparseMatrix.di.begin());
    std::copy(slau.matrix.ggl, slau.matrix.ggl + size, SparseMatrix.ggl.begin());
    std::copy(slau.matrix.ggu,slau.matrix.ggu + size, SparseMatrix.ggu.begin() );
    std::copy(slau.matrix.ig,slau.matrix.ig + N+1, SparseMatrix.ig.begin() );
    std::copy(slau.matrix.jg,slau.matrix.jg + size, SparseMatrix.jg.begin() );
    
    SparseMatrix.PrintDenseMatrix("Sparse.txt");
    // Конвертим в Профильный формат 
    SLAU_ProfileMatrix LUSolver;
    Matrix::SparseMatrix2ProfileMatrix(SparseMatrix, LUSolver.Matr);
    LUSolver.x.resize(N);
    LUSolver.f.resize(N);
    LUSolver.N = N;
    LUSolver.size = LUSolver.Matr.size;
    std::copy(slau.f, slau.f+N, LUSolver.f.begin());*/




    /* LU решатель */
    
    /*auto begin = std::chrono::steady_clock::now();
    SolveSlau(LUSolver, LUSolver.x);
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    //std::cout << elapsed_ms.count() << std::endl;

    SLAUSolvers::IterationSolvers::ActionVec(1.0, x, -1.0, LUSolver.x.data(), x, slau.N);
    double r = SLAUSolvers::IterationSolvers::Norma(x, slau.N)/x_true_norm;

    double Axf[12];
    SLAUSolvers::IterationSolvers::MultA(slau.matrix, LUSolver.x.data(), Axf);
    SLAUSolvers::IterationSolvers::ActionVec(1.0, Axf, -1.0, slau.f, Axf, slau.N);
    double r_Axf = SLAUSolvers::IterationSolvers::Norma(Axf, slau.N)/SLAUSolvers::IterationSolvers::Norma(slau.f, slau.N);
    */
    std::cout << "|  | " << elapsed_ms.count() << " |" << "| " << r << "|\n";

    return 0;
}