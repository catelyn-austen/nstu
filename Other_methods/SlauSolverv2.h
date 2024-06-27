#ifndef SLAU_SOLVER_V2_H_
#define SLAU_SOLVER_V2_H_

#include <vector>
#include <string>

#include "Matrix.h"


namespace SLAUSolvers
{
    enum class MatrixType
	{
		NONE, // Не установленный тип - default
		DENSE,
		SPARSE,
		SPARSE_SYMETRIC,
		PROFILE,
		PROFILE_SYMETRIC
	};

    /* Структуры для инициализации слау */
    namespace InitDataStructs
    {
        struct InitDataSparseMatrix
        {
            // Тип матрицы 
			MatrixType type = MatrixType::NONE;
		
            std::string gg = "";
			std::string ggl = "";
			std::string ggu = "";
			std::string di = "";
			std::string kuslau = "";
			std::string f = ""; // Вектор правой части 
			std::string ig = "";
			std::string jg = "";
        };

        struct InitDataProfileMatrix
        {
            // Тип матрицы 
			MatrixType type = MatrixType::NONE;
            
            std::string gg = "";
			std::string ggl = "";
			std::string ggu = "";
			std::string di = "";
			std::string kuslau = "";
			std::string f = ""; // Вектор правой части 
			std::string ig = "";
        };
    };

    /*  Методы прямых решателей 
        LU - decomposition
    */
	namespace ForwardSolvers
	{
        
	};

    /* Методы Итерационных решателей включает в себя следующие методы:
		ЛОС не симметричных
		ЛОС с факторизацией LU
        BCG не симметричных
        BCG с факторизацией LU
		Формат хранения матрицы разряженный строчно столбцовый  
	*/

    namespace IterationSolvers
	{

		class Solver
		{
			private:
				Matrix::SparseMatrix matrix;
				Matrix::SparseMatrix factor_matrix;
				std::vector<double> f; // Вектор правой части
				std::vector<double> x; // Вектор решенеия 

				void Load(const InitDataStructs::InitDataSparseMatrix &Initdata);

			public:

				Solver(const InitDataStructs::InitDataSparseMatrix &Initdata);
				Solver(const Matrix::SparseMatrix &matr): matrix(matr) {}

				double LOS_LU();
				double LOS();

				double BCG_LU();
				double BCG();
				
		};

    };


    class Solver
    {

    };
};


#endif