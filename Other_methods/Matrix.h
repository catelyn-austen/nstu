#ifndef MATRIX_H_
#define MATRIX_H_

#include <cctype>
#include <vector>
#include <iostream>

namespace Matrix
{

    class ProfileMatrix;

    class SparseMatrix
    {
    private:
            /* Friend func */
        friend void SparseMatrix2ProfileMatrix(SparseMatrix &src_Matr, ProfileMatrix &dst_Matr);

    public:
        int32_t N = -1;
        int32_t size = -1;
        std::vector<double> di;
        std::vector<double> ggu;
        std::vector<double> ggl;
        std::vector<int32_t> ig;
        std::vector<int32_t> jg;

    public:
        SparseMatrix() = default;

        SparseMatrix(int32_t N, int32_t size);

        void AllocMemory(int32_t N, int32_t size);

        void InsertElem(double elem, int32_t i, int32_t j);

        double GetElement(int32_t i, int32_t j);

        void PrintDenseMatrix(const std::string& filename = "", int32_t precession = 3);


        ~SparseMatrix() = default;
    
    };

    class ProfileMatrix
    {

    private:

        /* Friend func */
        friend void SparseMatrix2ProfileMatrix(SparseMatrix &src_Matr, ProfileMatrix &dst_Matr);

    public:
        int32_t N = -1;
        int32_t size = -1;

        std::vector<int> ia;
        std::vector<double> di;
        std::vector<double> al;
        std::vector<double> au;

    

    public:
        ProfileMatrix() = default;

        ProfileMatrix(int32_t N, int32_t size);
        void AllocateMemory(int32_t N, int32_t size);

        void InsertElem(double elem, int32_t i, int32_t j);

        double GetElement(int32_t i, int32_t j);

        void PrintDenseMatrix(const std::string& filename = "", int32_t precession = 3);


        ~ProfileMatrix() = default;

        
    };

    void SparseMatrix2ProfileMatrix(SparseMatrix &src_Matr, ProfileMatrix &dst_Matr);
}


#endif