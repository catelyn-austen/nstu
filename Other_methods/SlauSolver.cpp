#include "SlauSolver.h"

#include <iostream>
#include <fstream>
#include <exception>
#include <cmath>
#include <iomanip>
#include <stdio.h>

/* Реализация функций для итерационных решателей  */
namespace SLAUSolvers
{
	// Используем пространство имен итерационных решений 
	namespace IterationSolvers
	{
		using namespace std;
		/* Работа с памятью и загрузка из внешних источников */


		void Allocate_Memory_SLAU(SLAU& slau, const int N, const int size, const MatrixType type)
		{
			if (type == MatrixType::SPARSE)
			{
				try
				{
					slau.N = N;
					slau.size = size;
					slau.matrix.N = N;
					slau.matrix.size = size;
	

					if (slau.MainMemoryPool_ != nullptr)
						delete[] slau.MainMemoryPool_;

					slau.MainMemoryPool_ = new double[2 * N + 2 * size];
					double* MainMemoryPool_ = slau.MainMemoryPool_;
					int index = 0; // Индексация 

					/* Распределяем память массив double */
					slau.matrix.di = &MainMemoryPool_[index];
					index += N;
					slau.f = &MainMemoryPool_[index];
					index += N;
					slau.matrix.ggl = &MainMemoryPool_[index];
					index += size;
					slau.matrix.ggu = &MainMemoryPool_[index];
					
					/* Выделяем память под ig, jg */
					if (slau.matrix.ig != nullptr)
						delete[] slau.matrix.ig;

					if (slau.matrix.jg != nullptr)
						delete[] slau.matrix.jg;

					slau.matrix.ig = new int[N + 1];
					slau.matrix.jg = new int[size];
					int* ig = slau.matrix.ig;
					int* jg = slau.matrix.jg;

					/* Инициализация памяти нулями */
					for (int i = 0; i < 2 * N + 2 * size; i++)
						MainMemoryPool_[i] = 0;

					for (int i = 0; i < N + 1; i++)
						ig[i] = 0;

					for (int i = 0; i < size; i++)
						jg[i] = 0;
				}
				catch (bad_alloc &e)
				{
					cout <<  e.what();
					cout << "\n Работа программы прекращена\n";
					exit(-1);
				}
				catch (...)
				{
					cout << "Неизвестная ошибка\n";
					exit(-1);
				}
				
			}
			else if (type == MatrixType::SPARSE_SYMETRIC)
			{
				try
				{
					slau.N = N;
					slau.size = size;
					slau.matrix_s.N = N;
					slau.matrix_s.size = size;

					/* Псевдонимы */

					if (slau.MainMemoryPool_ != nullptr)
						delete[] slau.MainMemoryPool_;

					slau.MainMemoryPool_ = new double[2 * N + 1 * size];
					double* MainMemoryPool_ = slau.MainMemoryPool_;
					int index = 0; // Индексация 

					/* Распределяем память массив double */
					slau.matrix_s.di = &MainMemoryPool_[index];
					index += N;
					slau.f = &MainMemoryPool_[index];
					index += N;
					slau.matrix_s.gg = &MainMemoryPool_[index];
					

					/* Выделяем память под ig, jg */
					if (slau.matrix_s.ig != nullptr)
						delete[] slau.matrix_s.ig;

					if (slau.matrix_s.jg != nullptr)
						delete[] slau.matrix_s.jg;

					slau.matrix_s.ig = new int[N + 1];
					slau.matrix_s.jg = new int[size];

					int* ig = slau.matrix_s.ig;
					int* jg = slau.matrix_s.jg;

					/* Инициализация памяти нулями */
					for (int i = 0; i < 2 * N + 1 * size; i++)
						MainMemoryPool_[i] = 0;

					for (int i = 0; i < N + 1; i++)
						ig[i] = 0;

					for (int i = 0; i < size; i++)
						jg[i] = 0;
				}
				catch (bad_alloc& e)
				{
					cout << e.what();
					cout << "\n Работа программы прекращена\n";
					exit(-1);
				}
				catch (...)
				{
					cout << "Неизвестная ошибка\n";
					exit(-1);
				}
			}
			else
			{
				cout << "Invalid matrix format\n";
				exit(-1);
			}
		}

		void Allocate_Memory_Fact(SLAU& slau)
		{
			/* В зависимости от типа факторизации вводим вспомогательные данные */
			if (slau.mode == Solvemode::NONE)
			{
				cout << "Способ решения не задан\n";
				exit(-1);
			}


			/* Симметричные матрицы */
			if (slau.mode == Solvemode::LOS_SYMETRIC_DIAG_FACT || slau.mode == Solvemode::LOS_NOSYMETRIC_DIAG_FACT)
			{
				try
				{
					if (slau.FactMemoryPool_ != nullptr)
						delete[] slau.FactMemoryPool_;
					/* Нужен всего один вектор под диагональ */

					if (slau.mode == Solvemode::LOS_SYMETRIC_DIAG_FACT)
					{
						slau.Fmatr_s.N = slau.N; // Даем размер диагонали 
						slau.Fmatr_s.di = new double[slau.N];
						double* di = slau.Fmatr_s.di;
						slau.FactMemoryPool_ = di;
						for (int i = 0; i < slau.N; i++)
							di[i] = 0;
					}
					else if (slau.mode == Solvemode::LOS_NOSYMETRIC_DIAG_FACT)
					{
						slau.Fmatr.N = slau.N; // Даем размер диагонали 
						slau.Fmatr.di = new double[slau.N];
						double* di = slau.Fmatr.di;
						slau.FactMemoryPool_ = di;
						for (int i = 0; i < slau.N; i++)
							di[i] = 0;
					}
				
				}
				catch (bad_alloc& e)
				{
					cout << e.what();
					cout << "\n Работа программы прекращена\n";
					exit(-1);
				}
				catch (...)
				{
					cout << "Неизвестная ошибка\n";
					exit(-1);
				}
			} 
			else if (slau.mode == Solvemode::MSG_SYMETRIC_DIAG_FACT)
			{
				/* Нужен всего один вектор под диагональ причем что это полная копия*/
				slau.Fmatr_s.N = slau.N; // Даем размер диагонали 

			}
			else if (slau.mode == Solvemode::LOS_SYMETRIC_LLT_FACT || slau.mode == Solvemode::MSG_SYMETRIC_LLT_FACT)
			{
				try
				{
					if (slau.FactMemoryPool_ != nullptr)
						delete[]  slau.FactMemoryPool_;

					slau.Fmatr_s.N = slau.matrix_s.N;
					slau.Fmatr_s.size = slau.matrix_s.size;
					slau.Fmatr_s.ig = slau.matrix_s.ig;
					slau.Fmatr_s.jg = slau.matrix_s.jg;
					
					// нужно только под gg и di сумарно N + size
					slau.FactMemoryPool_ = new double[slau.N + slau.size];
					double* FactMemoryPool_ = slau.FactMemoryPool_;

					int index = 0;
					slau.Fmatr_s.di = &FactMemoryPool_[index];
					index += slau.N;
					slau.Fmatr_s.gg = &FactMemoryPool_[index];
					/* Обнуляем память */
					index += slau.size;

					for (int i = 0; i < index; i++)
						FactMemoryPool_[i] = 0;
				}
				catch (bad_alloc& e)
				{
					cout << e.what();
					cout << "\n Работа программы прекращена\n";
					exit(-1);
				}
				catch (...)
				{
					cout << "Неизвестная ошибка\n";
					exit(-1);
				}
			}

			/* Не симметричные матрицы */
			if (slau.mode == Solvemode::LOS_NOSYMETRIC_DIAG_FACT)
			{
				try {
					if (slau.FactMemoryPool_ != nullptr)
						delete[] slau.FactMemoryPool_;

					int N = slau.N;
					slau.FactMemoryPool_ = new double[N];
					double* FactMemoryPool_ = slau.FactMemoryPool_;
					slau.Fmatr.di = FactMemoryPool_;

					slau.Fmatr.N = slau.matrix.N;
					slau.Fmatr.size = slau.matrix.size;
					slau.Fmatr.ig = slau.matrix.ig;
					slau.Fmatr.jg = slau.matrix.jg;

					/* Обнуление памяти */
					for (int i = 0; i < N; i++)
						FactMemoryPool_[i] = 0;
				}
				catch (bad_alloc& e)
				{
					cout << e.what();
					cout << "\n Работа программы прекращена\n";
					exit(-1);
				}
				catch (...)
				{
					cout << "Неизвестная ошибка\n";
					exit(-1);
				}
			}
			else if ( 	slau.mode == Solvemode::LOS_NOSYMETRIC_LUsq_FACT ||
						slau.mode == Solvemode::LOS_NOSYMETRIC_LU_FACT   ||
					 	slau.mode == Solvemode::BCG_NOSYMETRIC_LU_FACT    
					)
			{
				try
				{
					if (slau.FactMemoryPool_ != nullptr)
						delete[] slau.FactMemoryPool_;

					slau.Fmatr.N = slau.matrix.N;
					slau.Fmatr.size = slau.matrix.size;
					slau.Fmatr.ig = slau.matrix.ig;
					slau.Fmatr.jg = slau.matrix.jg;
					
					int N = slau.matrix.N;
					int size = slau.matrix.size;
					/* di(N) ggl(size) ggu(size) => N + 2*size */
					slau.FactMemoryPool_ = new double[N + 2 * size];
					double* FactMemoryPool_ = slau.FactMemoryPool_;

					int index = 0;
					slau.Fmatr.di = &FactMemoryPool_[index];
					index += N;
					slau.Fmatr.ggl = &FactMemoryPool_[index];
					index += size;
					slau.Fmatr.ggu = &FactMemoryPool_[index];

					/* Зануляем память */
					for (int i = 0; i < N + 2 * size; i++)
						FactMemoryPool_[i] = 0;
				}
				catch (bad_alloc& e)
				{
					cout << e.what();
					cout << "\n Работа программы прекращена\n";
					exit(-1);
				}
				catch (...)
				{
					cout << "Неизвестная ошибка\n";
					exit(-1);
				}

			}
			
		}

		void Allocate_Memory_Support(SLAU& slau)
		{
			/* Распишем в файле необходимое количество векторов */
			/* Раздача векторов будет происходить непосредственно в своей функции  */

			try
			{
				if (slau.SuportMemoryPool_ != nullptr)
					delete[] slau.SuportMemoryPool_;


				/* Симметричные матрицы */
				if (slau.mode == Solvemode::MSG_SYMETRIC_CLASSIC)
				{
					/* Память под следующие вектора 
						rk - 
						zk - 
						Azk - A*zk
						r_true - истинная невязка f - A*xk - не нужно это дубликат rk
						итого 3*N
					*/
					slau.SuportMemoryPool_ = new double[3 * slau.N];
					double* SuportMemoryPool_ = slau.SuportMemoryPool_;
					// Зануляем память 
					for (int i = 0; i < 3 * slau.N; i++)
						SuportMemoryPool_[i] = 0;

				}
				else if (	slau.mode == Solvemode::MSG_SYMETRIC_DIAG_FACT ||
							slau.mode == Solvemode::MSG_SYMETRIC_LLT_FACT  ||
						 	slau.mode == Solvemode::LOS_SYMETRIC_CLASSIC   ||
							slau.mode == Solvemode::LOS_NOSYMETRIC_CLASSIC
						)
				{
					/* Память под следующие вектора МСГ
						rk -
						zk -
						Azk - A*zk
						Mrk - Произведение вектора M^-1*rk
						r_true - истинная невязка f - A*xk - не нужно это дубликат rk
						итого 4*N
					*/

				    /* Память под следующие вектора для классического ЛОС
						rk -
						zk -
						Ark - A*rk
						pk - 
						r_true - истинная невязка f - A*xk - не нужно это дубликат rk
						итого 4*N
					*/

					slau.SuportMemoryPool_ = new double[4 * slau.N];
					double* SuportMemoryPool_ = slau.SuportMemoryPool_;
					// Зануляем память 
					for (int i = 0; i < 4 * slau.N; i++)
						SuportMemoryPool_[i] = 0;
				}
				else if (	slau.mode == Solvemode::LOS_SYMETRIC_DIAG_FACT   ||
							slau.mode == Solvemode::LOS_SYMETRIC_LLT_FACT    ||
							slau.mode == Solvemode::LOS_NOSYMETRIC_LUsq_FACT ||
							slau.mode == Solvemode::LOS_NOSYMETRIC_DIAG_FACT ||
							slau.mode == Solvemode::LOS_NOSYMETRIC_LU_FACT
						)
				{
					/* Память под следующие вектора для классического ЛОС
						rk -
						zk -
						Ark - A*rk
						pk -
						Qrk - Q^-1*rk
						SAQrk = S^-1*A*Q^-1*rk
						r_true - истинная невязка f - A*xk 
						итого 7*N
					*/
					slau.SuportMemoryPool_ = new double[7 * slau.N];
					double* SuportMemoryPool_ = slau.SuportMemoryPool_;
					// Зануляем память 
					for (int i = 0; i < 7 * slau.N; i++)
						SuportMemoryPool_[i] = 0;
				}
				else if (	slau.mode == Solvemode::BCG_NOSYMETRIC_CLASSIC ||
						 	slau.mode == Solvemode::BCG_NOSYMETRIC_LU_FACT
						)
				{
					/* Память под следующие вектора для классического ЛОС
						rk -
						pk -
						zk -
						sk -
						Azk - A*z(k-1)
						ATsk - A(T)*s(k-1)
						r_true - истинная невязка f - A*xk 
						итого 7*N
					*/
					slau.SuportMemoryPool_ = new double[7 * slau.N];
					double* SuportMemoryPool_ = slau.SuportMemoryPool_;
					// Зануляем память 
					for (int i = 0; i < 7 * slau.N; i++)
						SuportMemoryPool_[i] = 0;

				}
				/* Не симметричные матрицы */
				
			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}
		}

		void DeAllocateMemory(SLAU& slau)
		{
			if (slau.MainMemoryPool_ != nullptr)
				delete[] slau.MainMemoryPool_;

			if (slau.FactMemoryPool_ != nullptr)
				delete[] slau.FactMemoryPool_;
			
			if (slau.SuportMemoryPool_ != nullptr)
				delete[] slau.SuportMemoryPool_;

			if (slau.matrix.N != -1)
			{
				if (slau.matrix.ig != nullptr)
					delete[] slau.matrix.ig;

				if (slau.matrix.jg != nullptr)
					delete[] slau.matrix.jg;
			}

			if (slau.matrix_s.N != -1)
			{
				if (slau.matrix_s.ig != nullptr)
					delete[] slau.matrix_s.ig;

				if (slau.matrix_s.jg != nullptr)
					delete[] slau.matrix_s.jg;
			}
		}

		void Load(InitDateFile& initfile, SLAU& slau, bool PrintProgress)
		{

			auto condPrint = [](const string& str, bool flag) { if (flag) cout << str << "\n"; };
		
			/* Произведем загрузку с файлов */
			ifstream in;
			in.open(initfile.kuslau);
			in >> slau.N >> slau.eps >> slau.maxiter;
			in.close();
			
			if (initfile.type == MatrixType::SPARSE)
			{
				try {
					bool is_CFmt = false; // Какой формат входных данных true - С формат false - fortran format
					/* Псевдонимы */
					
					slau.matrix.N = slau.N;
					int N = slau.N;
					
					/* ig */
					in.open(initfile.ig);

					slau.matrix.ig = new int[N + 1];
					int* ig = slau.matrix.ig;
					in >> ig[0];
					if (ig[0] == 1)
					{
						is_CFmt = false;
						ig[0]--;
					}
					else is_CFmt = true;

					if (!is_CFmt)
					{
						// Делаем декремент на 1
						for (int i = 1; i <= N; i++)
						{
							in >> ig[i];
							ig[i]--;
						}
					}
					else
					{
						// Без декремента 
						for (int i = 1; i <= N; i++)
							in >> ig[i];
					}

					
					in.close();
					condPrint("ig Загружен", PrintProgress);
					/******************/

					/* Вычислим размер size for gg, jg */
					int size = ig[N];
					slau.size = size;
					slau.matrix.size = size;

					/* jg */
					slau.matrix.jg = new int[size];
					int* jg = slau.matrix.jg;
					in.open(initfile.jg);
					if (!is_CFmt)
					{
						for (int i = 0; i < size; i++)
						{
							in >> jg[i];
							jg[i]--;
						}
					}
					else
					{
						for (int i = 0; i < size; i++)
							in >> jg[i];
					}
					in.close();
					condPrint("jg Загружен", PrintProgress);
					/******************/
					/* Выделяем память под di(N), gg(size), f(N)  */
					slau.MainMemoryPool_ = new double[2 * N + 2*size];
					double* MainMemoryPool_ = slau.MainMemoryPool_;
					int index = 0;
					slau.matrix.di = &MainMemoryPool_[index];
					index += N;
					slau.f = &MainMemoryPool_[index];
					index += N;
					slau.matrix.ggl = &MainMemoryPool_[index];
					index += size;
					slau.matrix.ggu = &MainMemoryPool_[index];

					double* di = slau.matrix.di;
					double* ggl = slau.matrix.ggl;
					double* ggu = slau.matrix.ggu;
					double* f = slau.f;

					/* di */
					in.open(initfile.di);
					for (int i = 0; i < N; i++)
						in >> di[i];
					in.close();
					condPrint("di Загружен", PrintProgress);


	
					/****************/

					/* f */
					in.open(initfile.f);
					for (int i = 0; i < N; i++)
						in >> f[i];
					in.close();
					condPrint("f Загружен", PrintProgress);
					/***************/

					/* ggl */
					in.open(initfile.ggl);
					for (int i = 0; i < size; i++)
						in >> ggl[i];
					in.close();
					condPrint("ggl Загружен", PrintProgress);
					/**************/

					/* ggu */
					in.open(initfile.ggu);
					for (int i = 0; i < size; i++)
						in >> ggu[i];
					in.close();
					condPrint("ggu Загружен", PrintProgress);
					/**************/
				}
				catch (bad_alloc& e)
				{
					cout << e.what();
					cout << "\n Работа программы прекращена\n";
					exit(-1);
				}
				catch (...)
				{
					cout << "Неизвестная ошибка\n";
					exit(-1);
				}
			}
			else if (initfile.type == MatrixType::SPARSE_SYMETRIC)
			{
				try {
					bool is_CFmt = false; // Какой формат входных данных true - С формат false - fortran format
					/* Псевдонимы */
					slau.matrix_s.N = slau.N;
					int N = slau.N;

					/* ig */
					in.open(initfile.ig);

					slau.matrix_s.ig = new int[N + 1];
					int* ig = slau.matrix_s.ig;
					in >> ig[0];
					if (ig[0] == 1)
					{
						is_CFmt = false;
						ig[0]--;
					}
					else is_CFmt = true;

					if (!is_CFmt)
					{
						// Делаем декремент на 1
						for (int i = 1; i <= N; i++)
						{
							in >> ig[i];
							ig[i]--;
						}
					}
					else
					{
						// Без декремента 
						for (int i = 1; i <= N; i++)
							in >> ig[i];
					}
					in.close();
					condPrint("ig Загружен", PrintProgress);
					/******************/

					/* Вычислим размер size for gg, jg */
					int size = ig[N];
					slau.size = size;
					slau.matrix_s.size = size;

					/* jg */
					slau.matrix_s.jg = new int[size];
					int* jg = slau.matrix_s.jg;
					in.open(initfile.jg);
					if (!is_CFmt)
					{
						for (int i = 0; i < size; i++)
						{
							in >> jg[i];
							jg[i]--;
						}
					}
					else
					{
						for (int i = 0; i < size; i++)
							in >> jg[i];
					}
					in.close();
					condPrint("jg Загружен", PrintProgress);
					/******************/
					/* Выделяем память под di(N), gg(size), f(N)  */
					slau.MainMemoryPool_ = new double[2 * N + size];
					double* MainMemoryPool_ = slau.MainMemoryPool_;


					int index = 0;
					slau.matrix_s.di = &MainMemoryPool_[index];
					index += N;
					slau.f = &MainMemoryPool_[index];
					index += N;
					slau.matrix_s.gg = &MainMemoryPool_[index];
					

					double* di = slau.matrix_s.di;
					double* gg = slau.matrix_s.gg;
					double* f = slau.f;
					/* di */
					in.open(initfile.di);
					for (int i = 0; i < N; i++)
						in >> di[i];
					in.close();
					condPrint("di Загружен", PrintProgress);
					/****************/

					/* f */
					in.open(initfile.f);
					for (int i = 0; i < N; i++)
						in >> f[i];
					in.close();
					condPrint("f Загружен", PrintProgress);
					/***************/

					/* gg */
					in.open(initfile.gg);
					for (int i = 0; i < size; i++)
						in >> gg[i];
					in.close();
					condPrint("gg Загружен", PrintProgress);
					/**************/
				}
				catch (bad_alloc& e)
				{
					cout << e.what();
					cout << "\n Работа программы прекращена\n";
					exit(-1);
				}
				catch (...)
				{
					cout << "Неизвестная ошибка\n";
					exit(-1);
				}

			}
			else
			{
				cout << "Некорректный формат матрицы \n";
				exit(-1);
			}
		}

		void SetSolveMode(SLAU& slau, Solvemode mode)
		{
			slau.mode = mode;
			Allocate_Memory_Fact(slau); // Память для факторизаций 
			Allocate_Memory_Support(slau); // Память для вспомогательных векторов 
		}

		/*****************************************/

		/* Работа С матрицами и векторами */

		void MultA(const Sparse_matrix_symetric& matr, const double* x, double* res)
		{
			double* gg = matr.gg;
			double* di = matr.di;
			int* ig = matr.ig;
			int* jg = matr.jg;
			int N = matr.N;
			int size = matr.size;

			for (int i = 0; i < N; i++)
			{
				res[i] = di[i] * x[i];
				int igi0 = ig[i];
				int igi1 = ig[i + 1];
				for (int j = igi0; j < igi1; j++)
				{
					int k = jg[j];
					res[i] += gg[j] * x[k];
					res[k] += gg[j] * x[i];
				}
			}
		}

		void MultA(const Sparse_matrix& matr, const double* x, double* res)
		{
			double* ggl = matr.ggl;
			double* ggu = matr.ggu;
			double* di = matr.di;
			int* ig = matr.ig;
			int* jg = matr.jg;
			int N = matr.N;
			int size = matr.size;

			for (int i = 0; i < N; i++)
			{
				res[i] = di[i] * x[i];
				int igi0 = ig[i];
				int igi1 = ig[i + 1];
				for (int j = igi0; j < igi1; j++)
				{
					int k = jg[j];
					res[i] += ggl[j] * x[k];
					res[k] += ggu[j] * x[i];
				}
			}
		}

		void MultAT(const Sparse_matrix& matr, const double* x, double* res)
		{
			double* ggl = matr.ggl;
			double* ggu = matr.ggu;
			double* di = matr.di;
			int* ig = matr.ig;
			int* jg = matr.jg;
			int N = matr.N;
			int size = matr.size;

			for (int i = 0; i < N; i++)
			{
				res[i] = di[i] * x[i];
				int igi0 = ig[i];
				int igi1 = ig[i + 1];
				for (int j = igi0; j < igi1; j++)
				{
					int k = jg[j];
					res[i] += ggu[j] * x[k];
					res[k] += ggl[j] * x[i];
				}
			}
		}

		void MultFactU(const Fact_matrix &matr, const double *x, double *res)
		{
			double *ggu = matr.ggu;
			double *di = matr.di;
			int* ig = matr.ig;
			int* jg = matr.jg;
			int N = matr.N;
			int size = matr.size;

			for (int i = 0; i < N; i++)
			{
				res[i] = 1.0 * x[i];
				int igi0 = ig[i];
				int igi1 = ig[i + 1];
				for (int j = igi0; j < igi1; j++)
				{
					int k = jg[j];
					res[k] += ggu[j] * x[i];
				}
			}
		}

		void ActionVec(const double a, const double* x1, const double b, const double* x2, double* y, const int N)
		{
			for (int i = 0; i < N; i++)
				y[i] = a * x1[i] + b * x2[i];
		}


		double ScalarMult(const double* x1, const double* x2, const int N)
		{
			double res = 0.0;
			for (int i = 0; i < N; i++)
				res += x1[i] * x2[i];
			return res;
		}

		double Norma(const double* x, const int N)
		{
			return sqrt(ScalarMult(x, x, N));
		}

		void CopyVec(const double* from, double* to, const int N)
		{
			for (int i = 0; i < N; i++)
				to[i] = from[i];
		}

		/* Прямой ход Гауса */
		void normal(Fact_matrix_symetric& matr, const double* b, double* res)
		{
			int N = matr.N;
			double* di_fact = matr.di;
			double* gg_fact = matr.gg;
			int* ig = matr.ig;
			int* jg = matr.jg;

			if(b != res)
				CopyVec(b, res, N);

			for (int i = 0; i < N; i++) {
				int igi0 = ig[i];
				int igi1 = ig[i + 1];
				for (int j = igi0; j < igi1; j++)
					res[i] -= gg_fact[j] * res[jg[j]];

				res[i] /= di_fact[i];
			}
		}



		void normal(Fact_matrix& matr, const double* b, double* res)
		{
			int N = matr.N;
			double* di_fact = matr.di;
			double* ggl_fact = matr.ggl;
			int* ig = matr.ig;
			int* jg = matr.jg;

			if (b != res)
				CopyVec(b, res, N);

			for (int i = 0; i < N; i++) {
				int igi0 = ig[i];
				int igi1 = ig[i + 1];
				for (int j = igi0; j < igi1; j++)
					res[i] -= ggl_fact[j] * res[jg[j]];

				res[i] = res[i] / di_fact[i];
			}
		}

		void normalUT(Fact_matrix& matr, const double* b, double* res)
		{
			int N = matr.N;
			double* di_fact = matr.di;
			double* ggu_fact = matr.ggu;
			int* ig = matr.ig;
			int* jg = matr.jg;

			if (b != res)
				CopyVec(b, res, N);

			for (int i = 0; i < N; i++) {
				int igi0 = ig[i];
				int igi1 = ig[i + 1];
				for (int j = igi0; j < igi1; j++)
					res[i] -= ggu_fact[j] * res[jg[j]];

				res[i] = res[i] /1.0; //di_fact[i];
			}
		}
		/* Обратный ход */
		void reverse(Fact_matrix_symetric& matr, const double* x, double* res)
		{
			int N = matr.N;
			double* di_fact = matr.di;
			double* gg_fact = matr.gg;
			int* ig = matr.ig;
			int* jg = matr.jg;

			if (x != res)
				CopyVec(x, res, N);

			for (int j = N - 1; j >= 0; j--) {
				res[j] = res[j] /1.0; //di_fact[j]; // LU modification

				int igj0 = ig[j];
				int igj1 = ig[j + 1];
				for (int i = igj0; i < igj1; i++)
				{
					int k = jg[i];
					res[k] -= gg_fact[i] * res[j];
				}
			}
		}

		/* 
			@warning  Используется модификация для работы с LU разложением. нужно пофиксить и сделать функции для работы с другими LU факторизациями
		*/
		void reverse(Fact_matrix& matr, const double* x, double* res)
		{
			int N = matr.N;
			double* di_fact = matr.di;
			double* ggu_fact = matr.ggu;
			int* ig = matr.ig;
			int* jg = matr.jg;

			if (x != res)
				CopyVec(x, res, N);

			for (int j = N - 1; j >= 0; j--) {
				res[j] = res[j] /1.0; //di_fact[j];  // LU modification 
				int igj0 = ig[j];
				int igj1 = ig[j + 1];
				for (int i = igj0; i < igj1; i++)
				{
					int k = jg[i];
					res[k] -= ggu_fact[i] * res[j];
				}
			}
		}

		void reverseLT(Fact_matrix &matr, const double *x, double *res)
		{
			int N = matr.N;
			double* di_fact = matr.di;
			double* ggl_fact = matr.ggl;
			int* ig = matr.ig;
			int* jg = matr.jg;

			if (x != res)
				CopyVec(x, res, N);

			for (int j = N - 1; j >= 0; j--) {
				res[j] = res[j] / di_fact[j];
				int igj0 = ig[j];
				int igj1 = ig[j + 1];
				for (int i = igj0; i < igj1; i++)
				{
					int k = jg[i];
					res[k] -= ggl_fact[i] * res[j];
				}
			}
		}

		/*****************************************/

		/* Факторизации */
		void DiagFactor(SLAU& slau)
		{
			/* Память уже выделена поэтому этим вопросом не занимаемся если мы здесь то точно все хорошо */
			if (slau.mode == Solvemode::MSG_SYMETRIC_DIAG_FACT)
			{
				// M = D
				//CopyVec(slau.matrix_s.di, slau.Fmatr_s.di, slau.N);
				slau.Fmatr_s.di = slau.matrix_s.di;
			}
			else if (slau.mode == Solvemode::LOS_NOSYMETRIC_DIAG_FACT)
			{
				// M = D = S*Q, где S = Q = sqrt(Dii) - матрицы предсталяют собой диагональ на которой стоят корни квадратные из элементов 
				int N = slau.N;
				double* di = slau.matrix.di;
				double* di_fact = slau.Fmatr.di;
				for (int i = 0; i < N; i++)
					di_fact[i] = sqrt(di[i]);
			}
			else if (slau.mode == Solvemode::LOS_SYMETRIC_DIAG_FACT)
			{
				// Ничем не отличается от варианта выше 
				int N = slau.N;
				double* di = slau.matrix_s.di;
				double* di_fact = slau.Fmatr_s.di;
				for (int i = 0; i < N; i++)
					di_fact[i] = sqrt(di[i]);
			}
			else
			{
				cout << "Неподходящий режим\n";
				exit(-1);
			}
		}

		void LLTFactor(SLAU& slau)
		{
			/* Проверка, что режим тот выбран */
			if (slau.mode == Solvemode::MSG_SYMETRIC_LLT_FACT || slau.mode == Solvemode::LOS_SYMETRIC_LLT_FACT)
			{
				/* Псевдонимы  */
				int N = slau.matrix_s.N;
				int size = slau.matrix_s.size;
				double* di = slau.matrix_s.di;
				double* gg = slau.matrix_s.gg;
				int* ig = slau.matrix_s.ig;
				int* jg = slau.matrix_s.jg;

				double* di_fact = slau.Fmatr_s.di;
				double* gg_fact = slau.Fmatr_s.gg;

				// Потом пересмотрю эту операцию 
				CopyVec(di, di_fact, N);
				CopyVec(gg, gg_fact, size);

				for (int i = 0; i < N; i++)
				{
					double sum_diag = 0;

					int igi0 = ig[i];
					int igi1 = ig[i + 1];
					for (int j = igi0; j < igi1; j++)
					{
						double sum = 0;
						int jgj0 = jg[j];
						int jk = ig[jgj0];
						int ik = ig[i];

						int igjgj1 = ig[jgj0 + 1];
						while ((ik < j) && (jk < igjgj1))
						{
							int l = jg[jk] - jg[ik];
							if (l == 0) {
								sum += gg_fact[jk] * gg_fact[ik];
								ik++; jk++;
							}
							jk += (l < 0);
							ik += (l > 0);
						}
						gg_fact[j] -= sum;
						gg_fact[j] /= di_fact[jg[j]];
						sum_diag += gg_fact[j] * gg_fact[j];
					}
					di_fact[i] -= sum_diag;
					di_fact[i] = sqrt(abs(di_fact[i]));
				}
			}
			else
			{
				cout << "Неподходящий режим\n";
				exit(-1);
			}
		}

		void LUsqFactor(SLAU& slau)
		{

			if (slau.mode == Solvemode::LOS_NOSYMETRIC_LUsq_FACT)
			{
				double sum_u, sum_l, sum_d;
				
				/* Псевдонимы  */
				int N = slau.matrix.N;
				int size = slau.matrix.size;
				double* di = slau.matrix.di;
				double* ggu = slau.matrix.ggu;
				double* ggl = slau.matrix.ggl;
				int* ig = slau.matrix.ig;
				int* jg = slau.matrix.jg;

				double* di_fact = slau.Fmatr.di;
				double* ggu_fact = slau.Fmatr.ggu;
				double* ggl_fact = slau.Fmatr.ggl;

				// Потом пересмотрю эту операцию 
				CopyVec(di, di_fact, N);
				CopyVec(ggu, ggu_fact, size);
				CopyVec(ggl, ggl_fact, size);

				for (int i = 0; i < N; ++i) {

					int i0 = ig[i];
					int i1 = ig[i + 1];


					for (int k = i0; k < i1; ++k) {

						int j = jg[k];
						int j0 = ig[j];
						int j1 = ig[j + 1];
						sum_l = 0;
						sum_u = 0;
						int ki = i0;
						int kj = j0;

						while (ki < k && kj < j1) {

							if (jg[ki] == jg[kj]) {
								sum_l += ggl_fact[ki] * ggu_fact[kj];
								sum_u += ggu_fact[ki] * ggl_fact[kj];
								ki++;
								kj++;
							}
							else {
								if (jg[ki] > jg[kj]) kj++;
								else ki++;
							}
						}

						ggl_fact[k] = (ggl_fact[k] - sum_l) / di_fact[j];
						ggu_fact[k] = (ggu_fact[k] - sum_u) / di_fact[j];
					}



					sum_d = 0.0;
					for (int k = i0; k < i1; ++k)
						sum_d += ggl_fact[k] * ggu_fact[k];
					di_fact[i] = sqrt(di_fact[i] - sum_d);
				}

			}
		}

		void LUFactor(SLAU& slau)
		{
			if (slau.mode == Solvemode::LOS_NOSYMETRIC_LU_FACT || slau.mode == Solvemode::BCG_NOSYMETRIC_LU_FACT)
			{
				double sum_u, sum_l, sum_d;
				
				/* Псевдонимы  */
				int N = slau.matrix.N;
				int size = slau.matrix.size;
				double* di = slau.matrix.di;
				double* ggu = slau.matrix.ggu;
				double* ggl = slau.matrix.ggl;
				int* ig = slau.matrix.ig;
				int* jg = slau.matrix.jg;

				double* di_fact = slau.Fmatr.di;
				double* ggu_fact = slau.Fmatr.ggu;
				double* ggl_fact = slau.Fmatr.ggl;

				// Потом пересмотрю эту операцию 
				CopyVec(di, di_fact, N);
				CopyVec(ggu, ggu_fact, size);
				CopyVec(ggl, ggl_fact, size);

				for (int i = 0; i < N; ++i) {

					int i0 = ig[i];
					int i1 = ig[i + 1];


					for (int k = i0; k < i1; ++k) {

						int j = jg[k];
						int j0 = ig[j];
						int j1 = ig[j + 1];
						sum_l = 0;
						sum_u = 0;
						int ki = i0;
						int kj = j0;

						while (ki < k && kj < j1) {

							if (jg[ki] == jg[kj]) {
								sum_l += ggl_fact[ki] * ggu_fact[kj];
								sum_u += ggu_fact[ki] * ggl_fact[kj];
								ki++;
								kj++;
							}
							else {
								if (jg[ki] > jg[kj]) kj++;
								else ki++;
							}
						}

						ggl_fact[k] = (ggl_fact[k] - sum_l);
						ggu_fact[k] = (ggu_fact[k] - sum_u) / di_fact[j];
					}



					sum_d = 0.0;
					for (int k = i0; k < i1; ++k)
						sum_d += ggl_fact[k] * ggu_fact[k];
 					di_fact[i] = di_fact[i] - sum_d;
				}

			}
		}
		/********************************/


		/* MSG симетричные матрицы  */

		/* Частные решатели для каждого типа решателя */
		double MSG_Symetric_Classic(SLAU& slau, bool printIteration)
		{
			auto printIterationProcess = [printIteration](double res, int k ) 
			{
				if(printIteration)
					cout << "Iteration k = " << k+1 << "  || f-A*xk ||/|| f || = " << res << "\n";
			};
			double res = 0.0; // Возвращаемое значение невязки 
			int index = 0;
			int N = slau.N;

			/* Вспомогательные вектора */
			double* rk = &slau.SuportMemoryPool_[index];
			index += N;
			double* zk = &slau.SuportMemoryPool_[index];
			index += N;
			double *Azk = &slau.SuportMemoryPool_[index];
			double* r_true = rk;
			/* Псевдонимы для используемых векторов */
			double* f = slau.f;
			Sparse_matrix_symetric& matr = slau.matrix_s;

			/* Вспомогательные величины */
			int maxiter = slau.maxiter;
			double eps = slau.eps;
			double f_norm = Norma(f, N);

			/* Если x нет то производим инициализацию */
			try
			{
				// Если начального приближения не было то 
				double* xk;
				if (slau.x == nullptr)
				{

					slau.x = new double[N];
					xk = slau.x;
					/* Зануление */
					for (int i = 0; i < N; i++)
						xk[i] = 0;

					/* Проводим инициализацию векторов */
					CopyVec(f, rk, N); // r0 = b
					CopyVec(rk, zk, N); // z0 = r0
				}
				else
				{
					xk = slau.x;

					MultA(matr, xk, rk); // A*x0
					ActionVec(1, f, -1, rk, rk, N); // r0 = f - A*x0

					CopyVec(rk, zk, N); // z0 = r0
				}

				CopyVec(rk, r_true, N);
				double r_true_norm = Norma(r_true, N);
				res = r_true_norm / f_norm;
				printIterationProcess(res, 0);
				if (res < eps) return res;

				/* Итерационный процесс */
				for (int k = 1; k < maxiter; k++)
				{
					MultA(matr, zk, Azk); // A*z(k-1)
					double dot_rk = ScalarMult(rk, rk, N); // (r(k-1); r(k-1) )
					double alpha_k = dot_rk / ScalarMult(Azk, zk, N); // (r(k-1); r(k-1) ) / (A*z(k-1); z(k-1))
					ActionVec(1, xk, alpha_k, zk, xk, N); // xk = x(k-1) + alpha_k*z(k-1)
					ActionVec(1, rk, -alpha_k, Azk, rk, N); // rk = r(k-1) - alpha_k*A*z(k-1)
					double betta_k = ScalarMult(rk, rk, N) / dot_rk; // (rk; rk)/ dot_rk
					ActionVec(1, rk, betta_k, zk, zk, N); // zk = rk + betta_k*z(k-1)

					// Расчет невязки можно подумать как это более логично делать 
					MultA(matr, xk, r_true);
					ActionVec(1, f, -1, r_true, r_true, N); // r_true = f - A*xk
					double r_true_norm = Norma(r_true, N);
					
					/* Распечатка итерационного процесса с выводом невязки*/
					res = r_true_norm / f_norm;
					printIterationProcess(res, k);

					// Выход по невязке 
					if (res < eps) break;
				}

			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}
			return res;
		}

		double MSG_Symetric_DiagFact(SLAU& slau, bool printIteration)
		{
			auto printIterationProcess = [printIteration](double res, int k)
			{
				if (printIteration)
					cout << "Iteration k = " << k + 1 << "  || f-A*xk ||/|| f || = " << res << "\n";
			};

			double res = 0.0; // Возвращаемое значение невязки 
			int N = slau.N;

			/* Расчет M^-1*rk  для диагональной матрицы это всего лишь умножение */
			auto MultD = [N](Fact_matrix_symetric& Fmatr, const double *x, double*res) -> void 
			{
				double* di = Fmatr.di;
				for (int i = 0; i < N; i++)
					res[i] = x[i] / di[i];
			};

			/* Вспомогательные вектора */
			int index = 0;
			double* rk = &slau.SuportMemoryPool_[index];
			index += N;
			double* zk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Azk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Mrk = &slau.SuportMemoryPool_[index];
			double* r_true = rk;

			/* Псевдонимы для используемых векторов */
			double* f = slau.f;
			Sparse_matrix_symetric& matr = slau.matrix_s;
			Fact_matrix_symetric& Fmatr_s = slau.Fmatr_s;
			/* Вспомогательные величины */
			int maxiter = slau.maxiter;
			double eps = slau.eps;
			double f_norm = Norma(f, N);

			/* Если x нет то производим инициализацию */
			try
			{
				// Если начального приближения не было то 
				double* xk;
				if (slau.x == nullptr)
				{

					slau.x = new double[N];
					xk = slau.x;
					/* Зануление */
					for (int i = 0; i < N; i++)
						xk[i] = 0;

					/* Проводим инициализацию векторов */
					CopyVec(f, rk, N); // r0 = b
				}
				else
				{
					xk = slau.x;
					MultA(matr, xk, rk); // A*x0
					ActionVec(1, f, -1, rk, rk, N); // r0 = f - A*x0
				}
				MultD(Fmatr_s, rk, zk); // z0 = M^-1*r0
				CopyVec(zk, Mrk, N);
				CopyVec(rk, r_true, N);
				double r_true_norm = Norma(r_true, N);
				res = r_true_norm / f_norm;
				printIterationProcess(res, 0);
				if (res < eps) return res;

				for (int k = 1; k < maxiter; k++)
				{
					MultA(matr, zk, Azk); // A*z(k-1)
					double dot_Mrk = ScalarMult(Mrk, rk, N); // (M^-1*r(k-1); r(k-1))
					double alphak = dot_Mrk / ScalarMult(Azk, zk, N); // 

					ActionVec(1, xk, alphak, zk, xk, N); // xk = xk-1 + alpahak*zk-1
					ActionVec(1, rk, -alphak, Azk, rk, N); // rk = rk-1 - alphak*Azk-1

					MultD(Fmatr_s, rk, Mrk);
					double bettak = ScalarMult(Mrk, rk, N) / dot_Mrk; //

					ActionVec(1, Mrk, bettak, zk, zk, N);


					// Расчет невязки 
					MultA(matr, xk, r_true);
					ActionVec(1, f, -1, r_true, r_true, N); // r_true = f - A*xk
					double r_true_norm = Norma(r_true, N);

					/* Распечатка итерационного процесса с выводом невязки*/
					res = r_true_norm / f_norm;
					printIterationProcess(res, k);

					// Выход по невязке 
					if (res < eps) break;
				}

			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}

			return res;
		}

		double MSG_Symetric_LLTFact(SLAU& slau, bool printIteration)
		{
			auto printIterationProcess = [printIteration](double res, int k)
			{
				if (printIteration)
					cout << "Iteration k = " << k + 1 << "  || f-A*xk ||/|| f || = " << res << "\n";
			};

			double res = 0.0; // Возвращаемое значение невязки 
			int N = slau.N;

			/* Расчет M^-1*rk  Эквивалентно решению слау по гаусу */
			auto MultD = [N](Fact_matrix_symetric& Fmatr, const double* x, double* res) -> void
			{
				normal(Fmatr, x, res);
				reverse(Fmatr, res, res);
			};

			LLTFactor(slau); // Факторизация 

			/* Вспомогательные вектора */
			int index = 0;
			double* rk = &slau.SuportMemoryPool_[index];
			index += N;
			double* zk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Azk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Mrk = &slau.SuportMemoryPool_[index];
			double* r_true = rk;
			

			/* Псевдонимы для используемых векторов */
			double* f = slau.f;
			Sparse_matrix_symetric& matr = slau.matrix_s;
			Fact_matrix_symetric& Fmatr_s = slau.Fmatr_s;
			/* Вспомогательные величины */
			int maxiter = slau.maxiter;
			double eps = slau.eps;
			double f_norm = Norma(f, N);

			/* Если x нет то производим инициализацию */
			try
			{
				// Если начального приближения не было то 
				double* xk;
				if (slau.x == nullptr)
				{

					slau.x = new double[N];
					xk = slau.x;
					/* Зануление */
					for (int i = 0; i < N; i++)
						xk[i] = 0;

					/* Проводим инициализацию векторов */
					CopyVec(f, rk, N); // r0 = b
				}
				else
				{
					xk = slau.x;
					MultA(matr, xk, rk); // A*x0
					ActionVec(1, f, -1, rk, rk, N); // r0 = f - A*x0
					
				}
				MultD(Fmatr_s, rk, zk); // z0 = M^-1*r0
				CopyVec(zk, Mrk, N);
				CopyVec(rk, r_true, N);
				double r_true_norm = Norma(r_true, N);
				res = r_true_norm / f_norm;
				printIterationProcess(res, 0);
				if (res < eps) return res;

				for (int k = 1; k < maxiter; k++)
				{
					MultA(matr, zk, Azk); // A*z(k-1)
					double dot_Mrk = ScalarMult(Mrk, rk, N); // (M^-1*r(k-1); r(k-1))
					double alphak = dot_Mrk / ScalarMult(Azk, zk, N); // 

					ActionVec(1, xk, alphak, zk, xk, N); // xk = xk-1 + alpahak*zk-1
					ActionVec(1, rk, -alphak, Azk, rk, N); // rk = rk-1 - alphak*Azk-1

					MultD(Fmatr_s, rk, Mrk);
					double bettak = ScalarMult(Mrk, rk, N) / dot_Mrk; //

					ActionVec(1, Mrk, bettak, zk, zk, N);


					// Расчет невязки 
					MultA(matr, xk, r_true);
					ActionVec(1, f, -1, r_true, r_true, N); // r_true = f - A*xk
					double r_true_norm = Norma(r_true, N);

					/* Распечатка итерационного процесса с выводом невязки*/
					res = r_true_norm / f_norm;
					printIterationProcess(res, k);

					// Выход по невязке 
					if (res < eps) break;
				}

			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}

			return res;
		}

		/********************************/

		/* MSG  НЕ симетричные матрицы - еще хороший вопрос а нужно ли делать MSG для не симетричных куда лучше сделать би-сопряженные градиенты */
		/********************************/

		/* LOS  симетричные матрицы */

		double LOS_Symetric_Classic(SLAU& slau, bool printIteration)
		{
			auto printIterationProcess = [printIteration](double res, int k)
			{
				if (printIteration)
					cout << "Iteration k = " << k + 1 << "  || f-A*xk ||/|| f || = " << res << "\n";
			};
			double res = 0.0; // Возвращаемое значение невязки 
			int index = 0;
			int N = slau.N;

			/* Вспомогательные вектора */
			double* rk = &slau.SuportMemoryPool_[index];
			index += N;
			double* zk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Ark = &slau.SuportMemoryPool_[index];
			index += N;
			double* pk = &slau.SuportMemoryPool_[index];
			double* r_true = rk; // Псевдоним 
			/* Псевдонимы для используемых векторов */
			double* f = slau.f;
			Sparse_matrix_symetric& matr = slau.matrix_s;

			/* Вспомогательные величины */
			int maxiter = slau.maxiter;
			double eps = slau.eps;
			double f_norm = Norma(f, N);
			
			/* Если x нет то производим инициализацию */
			try
			{
				// k = 0
				// Если начального приближения не было то 
				double* xk;
				if (slau.x == nullptr)
				{

					slau.x = new double[N];
					xk = slau.x;
					/* Зануление */
					for (int i = 0; i < N; i++)
						xk[i] = 0;

					/* Проводим инициализацию векторов */
					CopyVec(f, rk, N); // r0 = b
				}
				else
				{
					xk = slau.x;
					MultA(matr, xk, rk); // A*x0
					ActionVec(1, f, -1, rk, rk, N); // r0 = f - A*x0
				}
				MultA(matr, rk, pk); // p0 = A*r0
				CopyVec(rk, zk, N);

				double r_true_norm = Norma(r_true, N);
				res = r_true_norm / f_norm;
				printIterationProcess(res, 0);
				if (res < eps) return res;

				for (int k = 1; k < maxiter; k++)
				{
					double dot_pk = ScalarMult(pk, pk, N); // (pk-1; pk-1)
					double alphak = ScalarMult(pk, rk, N) / dot_pk;
					ActionVec(1.0, xk, alphak, zk, xk, N);
					ActionVec(1.0, rk, -alphak, pk, rk ,N);
					MultA(matr, rk, Ark);
					double bettak = -1.0 * (ScalarMult(pk, Ark, N))/dot_pk;
					ActionVec(1.0, rk, bettak, zk, zk, N);
					ActionVec(1.0, Ark, bettak, pk, pk, N);

					// Расчет невязки 
					MultA(matr, xk, r_true);
					ActionVec(1, f, -1, r_true, r_true, N); // r_true = f - A*xk
					double r_true_norm = Norma(r_true, N);

					/* Распечатка итерационного процесса с выводом невязки*/
					res = r_true_norm / f_norm;
					printIterationProcess(res, k);

					// Выход по невязке 
					if (res < eps) break;
				}
			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}

			return res;
		}

		double LOS_Symetric_DiagFact(SLAU& slau, bool printIteration)
		{
			auto printIterationProcess = [printIteration](double res, int k)
			{
				if (printIteration)
					cout << "Iteration k = " << k + 1 << "  || f-A*xk ||/|| f || = " << res << "\n";
			};

			double res = 0.0; // Возвращаемое значение невязки 
			int index = 0;
			int N = slau.N;

			/* Произведение Q^-1*vec диагональная матрица умножается на вектор  */
			auto MultD = [N](const double* D, const double* x, double* res)
			{
				for (int i = 0; i < N; i++)
					res[i] = x[i] / D[i];
			};

			/* Вспомогательные вектора */
			double* rk = &slau.SuportMemoryPool_[index];
			index += N;
			double* zk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Ark = &slau.SuportMemoryPool_[index];
			index += N;
			double* pk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Qrk = &slau.SuportMemoryPool_[index];
			index += N;
			double* SAQrk = &slau.SuportMemoryPool_[index];
			index += N;
			double* r_true = &slau.SuportMemoryPool_[index]; 

			/* Псевдонимы для используемых векторов */
			double* f = slau.f;
			DiagFactor(slau); // Факторизация 
			Sparse_matrix_symetric& matr = slau.matrix_s;
			double* FactDiag = slau.Fmatr_s.di; // Факторизованная диагональ 

			/* Вспомогательные величины */
			int maxiter = slau.maxiter;
			double eps = slau.eps;
			double f_norm = Norma(f, N);

			try
			{
				// k = 0
				// Если начального приближения не было то 
				double* xk;
				if (slau.x == nullptr)
				{

					slau.x = new double[N];
					xk = slau.x;
					/* Зануление */
					for (int i = 0; i < N; i++)
						xk[i] = 0;

					/* Проводим инициализацию векторов */
					CopyVec(f, rk, N); // r0 = b
					CopyVec(rk, r_true, N); // Копирование истинной невязки 
					MultD(FactDiag, rk, rk);
				}
				else
				{
					xk = slau.x;
					MultA(matr, xk, rk); // A*x0
					ActionVec(1, f, -1, rk, rk, N); // r0 = f - A*x0
					CopyVec(rk, r_true, N); // Копирование истинной невязки 
					MultD(FactDiag, rk, rk);
				}

				MultD(FactDiag, rk, zk); // zk = Q^-1*rk
				MultA(matr, zk, pk); // p0 = S^-1*A*r0
				MultD(FactDiag, pk, pk);


				double r_true_norm = Norma(r_true, N);
				res = r_true_norm / f_norm;
				printIterationProcess(res, 0);
				if (res < eps) return res;

				for (int k = 1; k < maxiter; k++)
				{
					double dot_pk = ScalarMult(pk, pk, N); // (pk-1; pk-1)
					double alphak = ScalarMult(pk, rk, N) / dot_pk;
					ActionVec(1.0, xk, alphak, zk, xk, N);
					ActionVec(1.0, rk, -alphak, pk, rk, N);

					MultD(FactDiag, rk, Qrk);
					MultA(matr, Qrk, SAQrk);
					MultD(FactDiag, SAQrk, SAQrk);
					double bettak = -1.0 * (ScalarMult(pk, SAQrk, N))/dot_pk;
					ActionVec(1.0, Qrk, bettak, zk, zk, N);
					ActionVec(1.0, SAQrk, bettak, pk, pk, N);

					// Расчет невязки 
					MultA(matr, xk, r_true);
					ActionVec(1, f, -1, r_true, r_true, N); // r_true = f - A*xk
					double r_true_norm = Norma(r_true, N);

					/* Распечатка итерационного процесса с выводом невязки*/
					res = r_true_norm / f_norm;
					printIterationProcess(res, k);

					// Выход по невязке 
					if (res < eps) break;
				}

			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}

			return res;
		}


		double LOS_Symetric_LLTFact(SLAU& slau, bool printIteration)
		{
			auto printIterationProcess = [printIteration](double res, int k)
			{
				if (printIteration)
					cout << "Iteration k = " << k + 1 << "  || f-A*xk ||/|| f || = " << res << "\n";
			};
			double res = 0.0; // Возвращаемое значение невязки 
			int index = 0;
			int N = slau.N;

			LLTFactor(slau); // Факторизация 

			/* Вспомогательные вектора */
			double* rk = &slau.SuportMemoryPool_[index];
			index += N;
			double* zk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Ark = &slau.SuportMemoryPool_[index];
			index += N;
			double* pk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Qrk = &slau.SuportMemoryPool_[index];
			index += N;
			double* SAQrk = &slau.SuportMemoryPool_[index];
			index += N;
			double* r_true = &slau.SuportMemoryPool_[index];

			/* Псевдонимы для используемых векторов */
			double* f = slau.f;
			Sparse_matrix_symetric& matr = slau.matrix_s;
			Fact_matrix_symetric& Fmatr_s = slau.Fmatr_s;
			/* Вспомогательные величины */
			int maxiter = slau.maxiter;
			double eps = slau.eps;
			double f_norm = Norma(f, N);

			try
			{
				// Если начального приближения не было то 
				double* xk;
				if (slau.x == nullptr)
				{

					slau.x = new double[N];
					xk = slau.x;
					/* Зануление */
					for (int i = 0; i < N; i++)
						xk[i] = 0;

					/* Проводим инициализацию векторов */
					CopyVec(f, rk, N); // r0 = b
				}
				else
				{
					xk = slau.x;
					MultA(matr, xk, rk); // A*x0
					ActionVec(1, f, -1, rk, rk, N); // r0 = f - A*x0

				}
				CopyVec(rk, r_true, N); // Истинная невязка еще до перехода к другим матрицам 

				normal(Fmatr_s, rk, rk); // S^-1*rk
				reverse(Fmatr_s, rk, zk); // zk = Q^-1*rk
				MultA(matr, zk, pk); // A*z0
				normal(Fmatr_s, pk, pk); // pk = S^-1*A*z0 

				
				double r_true_norm = Norma(r_true, N);
				res = r_true_norm / f_norm;
				printIterationProcess(res, 0);
				if (res < eps) return res;

				for (int k = 1; k < maxiter; k++)
				{
					double dot_pk = ScalarMult(pk, pk, N); // (pk-1; pk-1)
					double alphak = ScalarMult(pk, rk, N) / dot_pk;
					ActionVec(1.0, xk, alphak, zk, xk, N);
					ActionVec(1.0, rk, -alphak, pk, rk, N);

					/* Вычесляем произведения обратных матриц на вектор */
					reverse(Fmatr_s, rk, Qrk); // Q^-1*rk
					MultA(matr, Qrk, SAQrk);
					normal(Fmatr_s, SAQrk, SAQrk); // S^-1*A*Q^-1*rk


					double bettak = -1.0 * (ScalarMult(pk, SAQrk, N)) / dot_pk;
					ActionVec(1.0, Qrk, bettak, zk, zk, N);
					ActionVec(1.0, SAQrk, bettak, pk, pk, N);

					// Расчет невязки 
					MultA(matr, xk, r_true);
					ActionVec(1, f, -1, r_true, r_true, N); // r_true = f - A*xk
					double r_true_norm = Norma(r_true, N);

					/* Распечатка итерационного процесса с выводом невязки*/
					res = r_true_norm / f_norm;
					printIterationProcess(res, k);

					// Выход по невязке 
					if (res < eps) break;
				}
			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}


			return res;
		}
		/********************************/


		/* LOS  НЕ симетричные матрицы */
		double LOS_Classic(SLAU& slau, bool printIteration)
		{
			auto printIterationProcess = [printIteration](double res, int k)
			{
				if (printIteration)
					cout << "Iteration k = " << k + 1 << "  || f-A*xk ||/|| f || = " << res << "\n";
			};
			double res = 0.0; // Возвращаемое значение невязки 
			int index = 0;
			int N = slau.N;

			/* Вспомогательные вектора */
			double* rk = &slau.SuportMemoryPool_[index];
			index += N;
			double* zk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Ark = &slau.SuportMemoryPool_[index];
			index += N;
			double* pk = &slau.SuportMemoryPool_[index];
			double* r_true = rk; // Псевдоним 
			/* Псевдонимы для используемых векторов */
			double* f = slau.f;
			Sparse_matrix& matr = slau.matrix;

			/* Вспомогательные величины */
			int maxiter = slau.maxiter;
			double eps = slau.eps;
			double f_norm = Norma(f, N);

			/* Если x нет то производим инициализацию */
			try
			{
				// k = 0
				// Если начального приближения не было то 
				double* xk;
				if (slau.x == nullptr)
				{

					slau.x = new double[N];
					xk = slau.x;
					/* Зануление */
					for (int i = 0; i < N; i++)
						xk[i] = 0;

					/* Проводим инициализацию векторов */
					CopyVec(f, rk, N); // r0 = b
				}
				else
				{
					xk = slau.x;
					MultA(matr, xk, rk); // A*x0
					ActionVec(1, f, -1, rk, rk, N); // r0 = f - A*x0
				}
				MultA(matr, rk, pk); // p0 = A*r0
				CopyVec(rk, zk, N);

				double r_true_norm = Norma(r_true, N);
				res = r_true_norm / f_norm;
				printIterationProcess(res, 0);
				if (res < eps) return res;

				for (int k = 1; k < maxiter; k++)
				{
					double dot_pk = ScalarMult(pk, pk, N); // (pk-1; pk-1)
					double alphak = ScalarMult(pk, rk, N) / dot_pk;
					ActionVec(1.0, xk, alphak, zk, xk, N);
					ActionVec(1.0, rk, -alphak, pk, rk, N);
					MultA(matr, rk, Ark);
					double bettak = -1.0 * (ScalarMult(pk, Ark, N)) / dot_pk;
					ActionVec(1.0, rk, bettak, zk, zk, N);
					ActionVec(1.0, Ark, bettak, pk, pk, N);

					// Расчет невязки 
					MultA(matr, xk, r_true);
					ActionVec(1, f, -1, r_true, r_true, N); // r_true = f - A*xk
					double r_true_norm = Norma(r_true, N);

					/* Распечатка итерационного процесса с выводом невязки*/
					res = r_true_norm / f_norm;
					printIterationProcess(res, k);

					// Выход по невязке 
					if (res < eps) break;
				}
			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}

			return res;
		}

		double LOS_DiagFact(SLAU& slau, bool printIteration)
		{
			auto printIterationProcess = [printIteration](double res, int k)
			{
				if (printIteration)
					cout << "Iteration k = " << k + 1 << "  || f-A*xk ||/|| f || = " << res << "\n";
			};

			double res = 0.0; // Возвращаемое значение невязки 
			int index = 0;
			int N = slau.N;

			/* Произведение Q^-1*vec диагональная матрица умножается на вектор  */
			auto MultD = [N](const double* D, const double* x, double* res)
			{
				for (int i = 0; i < N; i++)
					res[i] = x[i] / D[i];
			};

			/* Вспомогательные вектора */
			double* rk = &slau.SuportMemoryPool_[index];
			index += N;
			double* zk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Ark = &slau.SuportMemoryPool_[index];
			index += N;
			double* pk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Qrk = &slau.SuportMemoryPool_[index];
			index += N;
			double* SAQrk = &slau.SuportMemoryPool_[index];
			index += N;
			double* r_true = &slau.SuportMemoryPool_[index];

			/* Псевдонимы для используемых векторов */
			double* f = slau.f;
			DiagFactor(slau); // Факторизация 
			Sparse_matrix& matr = slau.matrix;
			double* FactDiag = slau.Fmatr.di; // Факторизованная диагональ 

		
			/* Вспомогательные величины */
			int maxiter = slau.maxiter;
			double eps = slau.eps;
			double f_norm = Norma(f, N);

			try
			{
				// k = 0
				// Если начального приближения не было то 
				double* xk;
				if (slau.x == nullptr)
				{

					slau.x = new double[N];
					xk = slau.x;
					/* Зануление */
					for (int i = 0; i < N; i++)
						xk[i] = 0;

					/* Проводим инициализацию векторов */
					CopyVec(f, rk, N); // r0 = b
					CopyVec(rk, r_true, N); // Копирование истинной невязки 
					MultD(FactDiag, rk, rk);
				}
				else
				{
					xk = slau.x;
					MultA(matr, xk, rk); // A*x0
					ActionVec(1, f, -1, rk, rk, N); // r0 = f - A*x0
					CopyVec(rk, r_true, N); // Копирование истинной невязки 
					MultD(FactDiag, rk, rk);
				}

				MultD(FactDiag, rk, zk); // zk = Q^-1*rk
				MultA(matr, zk, pk); // p0 = S^-1*A*r0
				MultD(FactDiag, pk, pk);


				double r_true_norm = Norma(r_true, N);
				res = r_true_norm / f_norm;
				printIterationProcess(res, 0);
				if (res < eps) return res;

				for (int k = 1; k < maxiter; k++)
				{
					double dot_pk = ScalarMult(pk, pk, N); // (pk-1; pk-1)
					double alphak = ScalarMult(pk, rk, N) / dot_pk;
					ActionVec(1.0, xk, alphak, zk, xk, N);
					ActionVec(1.0, rk, -alphak, pk, rk, N);

					MultD(FactDiag, rk, Qrk);
					MultA(matr, Qrk, SAQrk);
					MultD(FactDiag, SAQrk, SAQrk);
					double bettak = -1.0 * (ScalarMult(pk, SAQrk, N)) / dot_pk;
					ActionVec(1.0, Qrk, bettak, zk, zk, N);
					ActionVec(1.0, SAQrk, bettak, pk, pk, N);

					// Расчет невязки 
					MultA(matr, xk, r_true);
					ActionVec(1, f, -1, r_true, r_true, N); // r_true = f - A*xk
					double r_true_norm = Norma(r_true, N);

					/* Распечатка итерационного процесса с выводом невязки*/
					res = r_true_norm / f_norm;
					printIterationProcess(res, k);

					// Выход по невязке 
					if (res < eps) break;
				}

			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}

			return res;
		}

		double LOS_LUsqFact(SLAU& slau, bool printIteration)
		{
			auto printIterationProcess = [printIteration](double res, int k)
			{
				if (printIteration)
					cout << "Iteration k = " << k + 1 << "  || f-A*xk ||/|| f || = " << res << "\n";
			};
			double res = 0.0; // Возвращаемое значение невязки 
			int index = 0;
			int N = slau.N;

			LUsqFactor(slau); // Факторизация 

			/* Вспомогательные вектора */
			double* rk = &slau.SuportMemoryPool_[index];
			index += N;
			double* zk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Ark = &slau.SuportMemoryPool_[index];
			index += N;
			double* pk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Qrk = &slau.SuportMemoryPool_[index];
			index += N;
			double* SAQrk = &slau.SuportMemoryPool_[index];
			index += N;
			double* r_true = &slau.SuportMemoryPool_[index];

			/* Псевдонимы для используемых векторов */
			double* f = slau.f;
			Sparse_matrix& matr = slau.matrix;
			Fact_matrix& Fmatr_s = slau.Fmatr;
			/* Вспомогательные величины */
			int maxiter = slau.maxiter;
			double eps = slau.eps;
			double f_norm = Norma(f, N);

			try
			{
				// Если начального приближения не было то 
				double* xk;
				if (slau.x == nullptr)
				{

					slau.x = new double[N];
					xk = slau.x;
					/* Зануление */
					for (int i = 0; i < N; i++)
						xk[i] = 0;

					/* Проводим инициализацию векторов */
					CopyVec(f, rk, N); // r0 = b
				}
				else
				{
					xk = slau.x;
					MultA(matr, xk, rk); // A*x0
					ActionVec(1, f, -1, rk, rk, N); // r0 = f - A*x0

				}
				CopyVec(rk, r_true, N); // Истинная невязка еще до перехода к другим матрицам 

				normal(Fmatr_s, rk, rk); // S^-1*rk
				reverse(Fmatr_s, rk, zk); // zk = Q^-1*rk
				MultA(matr, zk, pk); // A*z0
				normal(Fmatr_s, pk, pk); // pk = S^-1*A*z0 


				double r_true_norm = Norma(r_true, N);
				res = r_true_norm / f_norm;
				printIterationProcess(res, 0);
				if (res < eps) return res;

				for (int k = 1; k < maxiter; k++)
				{
					double dot_pk = ScalarMult(pk, pk, N); // (pk-1; pk-1)
					double alphak = ScalarMult(pk, rk, N) / dot_pk;
					ActionVec(1.0, xk, alphak, zk, xk, N);
					ActionVec(1.0, rk, -alphak, pk, rk, N);

					/* Вычесляем произведения обратных матриц на вектор */
					reverse(Fmatr_s, rk, Qrk); // Q^-1*rk
					MultA(matr, Qrk, SAQrk);
					normal(Fmatr_s, SAQrk, SAQrk); // S^-1*A*Q^-1*rk


					double bettak = -1.0 * (ScalarMult(pk, SAQrk, N)) / dot_pk;
					ActionVec(1.0, Qrk, bettak, zk, zk, N);
					ActionVec(1.0, SAQrk, bettak, pk, pk, N);

					// Расчет невязки 
					MultA(matr, xk, r_true);
					ActionVec(1, f, -1, r_true, r_true, N); // r_true = f - A*xk
					double r_true_norm = Norma(r_true, N);

					/* Распечатка итерационного процесса с выводом невязки*/
					res = r_true_norm / f_norm;
					printIterationProcess(res, k);

					// Выход по невязке 
					if (res < eps) break;
				}
			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}


			return res;
		}

		
		double LOS_LUFact(SLAU& slau, bool printIteration)
		{
			auto printIterationProcess = [printIteration](double res, int k)
			{
				if (printIteration)
					cout << "Iteration k = " << k + 1 << "  || f-A*xk ||/|| f || = " << res << "\n";
			};
			double res = 0.0; // Возвращаемое значение невязки 
			int index = 0;
			int N = slau.N;

			LUFactor(slau); // Факторизация 

			/* Вспомогательные вектора */
			double* rk = &slau.SuportMemoryPool_[index];
			index += N;
			double* zk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Ark = &slau.SuportMemoryPool_[index];
			index += N;
			double* pk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Qrk = &slau.SuportMemoryPool_[index];
			index += N;
			double* SAQrk = &slau.SuportMemoryPool_[index];
			index += N;
			double* r_true = &slau.SuportMemoryPool_[index];

			/* Псевдонимы для используемых векторов */
			double* f = slau.f;
			Sparse_matrix& matr = slau.matrix;
			Fact_matrix& Fmatr_s = slau.Fmatr;
			/* Вспомогательные величины */
			int maxiter = slau.maxiter;
			double eps = slau.eps;
			double f_norm = Norma(f, N);

			try
			{
				// Если начального приближения не было то 
				double* xk;
				if (slau.x == nullptr)
				{

					slau.x = new double[N];
					xk = slau.x;
					/* Зануление */
					for (int i = 0; i < N; i++)
						xk[i] = 0;

					/* Проводим инициализацию векторов */
					CopyVec(f, rk, N); // r0 = b
				}
				else
				{
					xk = slau.x;
					MultA(matr, xk, rk); // A*x0
					ActionVec(1, f, -1, rk, rk, N); // r0 = f - A*x0

				}
				CopyVec(rk, r_true, N); // Истинная невязка еще до перехода к другим матрицам 

				normal(Fmatr_s, rk, rk); // S^-1*rk
				reverse(Fmatr_s, rk, zk); // zk = Q^-1*rk
				MultA(matr, zk, pk); // A*z0
				normal(Fmatr_s, pk, pk); // pk = S^-1*A*z0 


				double r_true_norm = Norma(r_true, N);
				res = r_true_norm / f_norm;
				printIterationProcess(res, 0);
				if (res < eps) return res;

				for (int k = 1; k < maxiter; k++)
				{
					double dot_pk = ScalarMult(pk, pk, N); // (pk-1; pk-1)
					double alphak = ScalarMult(pk, rk, N) / dot_pk;
					ActionVec(1.0, xk, alphak, zk, xk, N);
					ActionVec(1.0, rk, -alphak, pk, rk, N);

					/* Вычесляем произведения обратных матриц на вектор */
					reverse(Fmatr_s, rk, Qrk); // Q^-1*rk
					MultA(matr, Qrk, SAQrk);
					normal(Fmatr_s, SAQrk, SAQrk); // S^-1*A*Q^-1*rk


					double bettak = -1.0 * (ScalarMult(pk, SAQrk, N)) / dot_pk;
					ActionVec(1.0, Qrk, bettak, zk, zk, N);
					ActionVec(1.0, SAQrk, bettak, pk, pk, N);

					// Расчет невязки 
					MultA(matr, xk, r_true);
					ActionVec(1, f, -1, r_true, r_true, N); // r_true = f - A*xk
					double r_true_norm = Norma(r_true, N);

					/* Распечатка итерационного процесса с выводом невязки*/
					res = r_true_norm / f_norm;
					printIterationProcess(res, k);

					// Выход по невязке 
					if (res < eps) break;
				}
			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}


			return res;

		}
		
		double BCG_Classic(SLAU& slau, bool printIteration)
		{
			auto printIterationProcess = [printIteration](double res, int k)
			{
				if (printIteration)
					cout << "Iteration k = " << k + 1 << "  || f-A*xk ||/|| f || = " << res << "\n";
			};
			double res = 0.0; // Возвращаемое значение невязки 
			int index = 0;
			int N = slau.N;

			/* Вспомогательные вектора */
			double* rk = &slau.SuportMemoryPool_[index];
			index += N;
			double* zk = &slau.SuportMemoryPool_[index];
			index += N;
			double* Azk = &slau.SuportMemoryPool_[index];
			index += N;
			double* pk = &slau.SuportMemoryPool_[index];
			index += N;
			double *sk = &slau.SuportMemoryPool_[index];
			index += N;
			double *ATsk = &slau.SuportMemoryPool_[index];

			double* r_true = rk; // Псевдоним 
			/* Псевдонимы для используемых векторов */
			double* f = slau.f;
			Sparse_matrix& matr = slau.matrix;

			/* Вспомогательные величины */
			int maxiter = slau.maxiter;
			double eps = slau.eps;
			double f_norm = Norma(f, N);

			/* Если x нет то производим инициализацию */
			try
			{
				// k = 0
				// Если начального приближения не было то 
				double* xk;
				if (slau.x == nullptr)
				{

					slau.x = new double[N];
					xk = slau.x;
					/* Зануление */
					for (int i = 0; i < N; i++)
						xk[i] = 0;

					/* Проводим инициализацию векторов */
					CopyVec(f, rk, N); // r0 = b
				}
				else
				{
					xk = slau.x;
					MultA(matr, xk, rk); // A*x0
					ActionVec(1, f, -1, rk, rk, N); // r0 = f - A*x0
				}

				CopyVec(rk, pk, N);
				CopyVec(rk, zk, N);
				CopyVec(rk, sk, N);

				double r_true_norm = Norma(r_true, N);
				res = r_true_norm / f_norm;
				printIterationProcess(res, 0);
				if (res < eps) return res;

				for (int k = 1; k < maxiter; k++)
				{
					/* Алгоритм */
					double prk = ScalarMult(pk, rk, N);  // (pk-1; rk-1)
					MultA(matr, zk, Azk);  // A*zk-1
					double alphak = prk/ScalarMult(sk, Azk, N); // alphak

					ActionVec(1.0, xk, alphak, zk, xk, N); // xk = xk-1 + alphak*zk-1
					ActionVec(1.0,rk, -alphak, Azk, rk, N); // rk = rk-1 - alphak*A*zk-1
					MultAT(matr, sk, ATsk); // A(T)*sk-1
					ActionVec(1.0, pk, -alphak, ATsk, pk, N); // pk = pk-1 - alphak*A(T)*sk-1
					double betak = ScalarMult(pk, rk, N)/prk;
					ActionVec(1.0, rk, betak, zk, zk, N);
					ActionVec(1.0, pk, betak, sk, sk, N);

					// Расчет невязки 
					MultA(matr, xk, r_true);
					ActionVec(1, f, -1, r_true, r_true, N); // r_true = f - A*xk
					double r_true_norm = Norma(r_true, N);

					/* Распечатка итерационного процесса с выводом невязки*/
					res = r_true_norm / f_norm;
					printIterationProcess(res, k);

					// Выход по невязке 
					if (res < eps) break;
				}

			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}

			return res;

		}

		double BCG_LUFact(SLAU& slau, bool printIteration)
		{
			auto printIterationProcess = [printIteration](double res, int k)
			{
				if (printIteration)
					cout << "Iteration k = " << k + 1 << "  || f-A*xk ||/|| f || = " << res << "\n";
			};
			double res = 0.0; // Возвращаемое значение невязки 
			int index = 0;
			int N = slau.N;

			LUFactor(slau); // Факторизация 


			/* Вспомогательные вектора */
			double* rk = &slau.SuportMemoryPool_[index];
			index += N;
			double* zk = &slau.SuportMemoryPool_[index];
			index += N;
			double* LAUzk = &slau.SuportMemoryPool_[index];
			index += N;
			double* pk = &slau.SuportMemoryPool_[index];
			index += N;
			double *sk = &slau.SuportMemoryPool_[index];
			index += N;
			double *UT_AT_LT_sk = &slau.SuportMemoryPool_[index];
			index += N;
			double *tmp_vec = &slau.SuportMemoryPool_[index]; // Используется для промежуточных вычислений 

			double* r_true = rk; // Псевдоним 
			/* Псевдонимы для используемых векторов */
			double* f = slau.f;
			Sparse_matrix& matr = slau.matrix;
			Fact_matrix& Fmatr_s = slau.Fmatr;

			/* Вспомогательные величины */
			int maxiter = slau.maxiter;
			double eps = slau.eps;
			double f_norm = Norma(f, N);

			/* Если x нет то производим инициализацию */

			try
			{
				// k = 0
				// Если начального приближения не было то 
				double* xk;
				if (slau.x == nullptr)
				{

					slau.x = new double[N];
					xk = slau.x;
					/* Зануление */
					for (int i = 0; i < N; i++)
						xk[i] = 0;

					/* Проводим инициализацию векторов */
					CopyVec(f, rk, N); // r0 = b
				}
				else
				{
					xk = slau.x;
					MultA(matr, xk, rk); // A*x0
					ActionVec(1, f, -1, rk, rk, N); // r0 = f - A*x0
					
				}

				normal(Fmatr_s, rk, rk); // L^-1*rk
				
				CopyVec(rk, pk, N);
				CopyVec(rk, zk, N);
				CopyVec(rk, sk, N);

				double r_true_norm = Norma(r_true, N);
				res = r_true_norm / f_norm;
				printIterationProcess(res, 0);
				if (res < eps) return res;

				for (int k = 1; k < maxiter; k++)
				{
					/* Алгоритм */
					double prk = ScalarMult(pk, rk, N);  // (pk-1; rk-1)
					/* Вычисляем произведение L^-1*A*U^-1*zk-1 */
					reverse(Fmatr_s, zk, LAUzk); // U^-1*zk-1
					MultA(matr, LAUzk, tmp_vec); // A*U^-1*zk-1
					normal(Fmatr_s, tmp_vec, LAUzk); // L^-1*A*U^-1*zk-1

					double alphak = prk/ScalarMult(sk, LAUzk, N);
					ActionVec(1.0, xk, alphak, zk, xk, N); // xk = xk-1 + alphak*zk-1
					ActionVec(1.0,rk, -alphak, LAUzk, rk, N); // rk = rk-1 - alphak*LAUzk*zk-1

					/*  Вычисляем произведение U^-T*A(T)*L^-T*sk-1 */
					reverseLT(Fmatr_s, sk, UT_AT_LT_sk); // L^-T*sk-1
					MultAT(matr, UT_AT_LT_sk, tmp_vec); // A(T)*L^-T*sk-1
					normalUT(Fmatr_s, tmp_vec, UT_AT_LT_sk); // U^-T*A(T)*L^-T*sk-1

					ActionVec(1.0, pk, -alphak, UT_AT_LT_sk, pk ,N); // pk = pk-1 - alphak*UT_AT_LT*sk-1
					double betak = ScalarMult(pk, rk, N)/prk;
					ActionVec(1.0, rk, betak, zk, zk, N);
					ActionVec(1.0, pk, betak, sk, sk, N);

					// Расчет невязки 
					reverse(Fmatr_s, slau.x, tmp_vec); // Получили вектор x

					MultA(matr, tmp_vec, r_true); // A*xk
					ActionVec(1, f, -1, r_true, r_true, N); // r_true = f - A*xk
					double r_true_norm = Norma(r_true, N); // Посчитали норму 

					normal(Fmatr_s, rk, rk); // L^-1*rk - сновва невязка предобусловленной сичтемы 
					/* Распечатка итерационного процесса с выводом невязки*/
					
					res = r_true_norm / f_norm;
					printIterationProcess(res, k);

					// Выход по невязке 
					if (res < eps) break;
				}	
			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}

			/* x* = U*xm */
			reverse(Fmatr_s, slau.x, slau.x);
			//CopyVec(tmp_vec, slau.x, N);

			return res;
		}

		/********************************/

		/* Общий решатель СЛАУ выбраной модификации  */
		double SolveSLAU(SLAU& slau, bool printIteration)
		{
			double r = 0.0;
			switch (slau.mode)
			{
			case Solvemode::MSG_SYMETRIC_CLASSIC:
				r = MSG_Symetric_Classic(slau, printIteration);
				break;
			case Solvemode::MSG_SYMETRIC_DIAG_FACT:
				r = MSG_Symetric_Classic(slau, printIteration);
				break;
			case Solvemode::MSG_SYMETRIC_LLT_FACT:
				r = MSG_Symetric_LLTFact(slau, printIteration);
				break;
			case Solvemode::LOS_SYMETRIC_CLASSIC:
				r = LOS_Symetric_Classic(slau, printIteration);
				break;
			case Solvemode::LOS_SYMETRIC_DIAG_FACT:
				r = LOS_Symetric_DiagFact(slau, printIteration);
				break;
			case Solvemode::LOS_SYMETRIC_LLT_FACT:
				r = LOS_Symetric_LLTFact(slau, printIteration);
				break;
			case Solvemode::LOS_NOSYMETRIC_CLASSIC:
				r = LOS_Classic(slau, printIteration);
				break;
			case Solvemode::LOS_NOSYMETRIC_DIAG_FACT:
				r = LOS_DiagFact(slau, printIteration);
				break;
			case Solvemode::LOS_NOSYMETRIC_LUsq_FACT:
				r = LOS_LUsqFact(slau, printIteration);
				break;
			case Solvemode::LOS_NOSYMETRIC_LU_FACT:
				r = LOS_LUFact(slau, printIteration);
				break;
			case Solvemode::BCG_NOSYMETRIC_CLASSIC:
				r = BCG_Classic(slau, printIteration);
				break;
			case Solvemode::BCG_NOSYMETRIC_LU_FACT:
				r = BCG_LUFact(slau, printIteration);
				break;
			default:
				cout << "Неизвестный способ решения СЛАУ";
				break;
			}

			return r;
		}

		/********************************/


		/* Вспомогательные функции для отладки  */
		ostream& operator <<(ostream& cnt, const SLAU& slau)
		{
			auto print = [&cnt](auto* arr, const char* name, int N)
			{
				cnt << name << " = [";
				for (int i = 0; i < N; i++)
					cnt << arr[i] << " ";
				cnt << "]\n";
			};


			cnt << "СЛАУ данные\n";
			cnt << "N = " << slau.N  << " size = " << slau.size << " eps = " << slau.eps << " maxiter = " << slau.maxiter << "\n";
			if (slau.matrix.N != -1)
			{
				print(slau.matrix.ig, "ig", slau.matrix.N+1);
				print(slau.matrix.jg, "jg", slau.matrix.size);
				print(slau.matrix.di, "di", slau.matrix.N);
				print(slau.matrix.ggl, "ggl", slau.matrix.size);
				print(slau.matrix.ggu, "ggu", slau.matrix.size);
				print(slau.f, "f", slau.matrix.N);

				/* Матрица факторизации */
				if (slau.Fmatr.N != -1)
				{
					print(slau.Fmatr.di, "di fact", slau.Fmatr.N);

					if (slau.mode != Solvemode::LOS_NOSYMETRIC_DIAG_FACT)
					{
						print(slau.Fmatr.ggu, "ggu fact", slau.Fmatr.size);
						print(slau.Fmatr.ggl, "ggl fact", slau.Fmatr.size);
					}
				}

			}
			else if (slau.matrix_s.N != -1)
			{
				print(slau.matrix_s.ig, "ig", slau.matrix_s.N+1);
				print(slau.matrix_s.jg, "jg", slau.matrix_s.size);
				print(slau.matrix_s.di, "di", slau.matrix_s.N);
				print(slau.matrix_s.gg, "gg", slau.matrix_s.size);
				print(slau.f, "f", slau.matrix_s.N);

				/* Матрица факторизации */
				if (slau.Fmatr_s.N != -1)
				{
					print(slau.Fmatr_s.di, "di fact", slau.Fmatr_s.N);
					print(slau.Fmatr_s.gg, "gg fact", slau.Fmatr_s.size);
				}
			}
			return cnt;
		}
	
		/* Работа с вектором X */


		void LoadX(SLAU& slau, const string filename)
		{
			try 
			{
				if (slau.x != nullptr)
					delete[] slau.x;

				slau.x = new double[slau.N];
				ifstream in;
				in.open(filename);
				for (int i = 0; i < slau.N; i++)
					in >> slau.x[i];
				in.close();
			}
			catch (bad_alloc& e)
			{
				cout << e.what();
				cout << "\n Работа программы прекращена\n";
				exit(-1);
			}
			catch (...)
			{
				cout << "Неизвестная ошибка\n";
				exit(-1);
			}
		}


		void SaveX(SLAU& slau, string filename, const string delimetr)
		{
			/* Потом переделаю так что бы можно было выставлять формат из вызова функции */
			ofstream out;
			out.open(filename);
			out << fixed;
			for (int i = 0; i < slau.N; i++)
				out  <<setprecision(15) << slau.x[i] << delimetr;
			out.close();
		}

		void PrintVec(double* vec, int N)
		{
			for (int i = 0; i < N; i++)
				cout << vec[i] << " ";
			cout << "\n";
		}

		void PrintX(SLAU& slau, const char *fmt)
		{
			for (int i = 0; i < slau.N; i++)
				printf(fmt, slau.x[i]);
		}

		/********************************/
	};
};