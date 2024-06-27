#include <iostream>
#include <vector>
#include <mpi.h>

#define MATRIX_SIZE 1000

using namespace std;

void printM(const vector<double>& matrix) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j)
            cout << matrix[i * MATRIX_SIZE + j] << "\t";
        cout << endl;
    }
    cout << endl;
}

void fillMatrix(vector<double>& A, vector<double>& B) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            A[i * MATRIX_SIZE + j] = i <= j ? 1 : 0;
            B[i * MATRIX_SIZE + j] = i >= j ? 1 : 0;
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv); // инициализация среды MPI

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // узнаем ранг (номер) процессора в коммуникаторе
    MPI_Comm_size(MPI_COMM_WORLD, &size); // узнаем кол-во процессоров в коммуникаторе

    vector<double> A(MATRIX_SIZE * MATRIX_SIZE);
    vector<double> B(MATRIX_SIZE * MATRIX_SIZE);
    vector<double> C(MATRIX_SIZE * MATRIX_SIZE, 0.0);

    if (rank == 0)
    {
        fillMatrix(A, B);
        /*cout << "Matrix A:" << endl;
        printM(A);
        cout << "Matrix B:" << endl;
        printM(B);*/
    }
        
    clock_t start_time = clock();
    // Рассылка А и В всем процессам
    // параметры: указатель на эл. буффера (в данном случае матрицы), размер буффера, тип эл., номер процесса
    // отправителя, коммуникатор
    MPI_Bcast(A.data(), MATRIX_SIZE * MATRIX_SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B.data(), MATRIX_SIZE * MATRIX_SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    vector<double> localC(MATRIX_SIZE * MATRIX_SIZE, 0.0); // вектор для рассылки процессам

    // Вычисление только той части матричного умножения, которая присвоена данному процессу
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        if (i % size == rank) { // если строка матрицы принадлежит процессу rank
            for (int j = 0; j < MATRIX_SIZE; ++j) {
                for (int k = 0; k < MATRIX_SIZE; ++k)
                {
                    localC[i * MATRIX_SIZE + j] += A[i * MATRIX_SIZE + k] * B[k * MATRIX_SIZE + j];
                }
            }
        }
    }

    // Объединение буффера каждого процесса локальных матриц localC в глобальную матрицу C с операцией сложения
    MPI_Reduce(localC.data(), C.data(), MATRIX_SIZE * MATRIX_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    clock_t end_time = clock();
    if (rank == 0) {
        /*cout << "Matrix C:" << endl;
        printM(C);*/
        cout << "Number of processes: " << size << endl;
    }
    double search_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    cout << search_time << endl;

    // норма С
    double norm = 0;
    for (int i = 0; i < MATRIX_SIZE; i++)
        for (int j = 0; j < MATRIX_SIZE; j++)
            norm = norm + pow(C[i * MATRIX_SIZE + j], 2);
    norm = sqrt(norm);
    cout << norm;

    MPI_Finalize();

    return 0;
}