#include "Matrix.h"
#include <iomanip>
#include <fstream>

namespace Matrix
{

    SparseMatrix::SparseMatrix(int32_t N, int32_t size)
    {
        this->N = N;
        this->size = size;
        ggu.resize(size);
        ggl.resize(size);
        di.resize(N);
        ig.resize(N + 1);
        jg.resize(size);

        // ig = {0,0,1,1,3,4,6};
        // jg = {0,0,1,3,1,2};
        // di = {1,5,8,12,15,19};
        // ggl = {4,10,11,14,16,17};
        // ggu = {2,3,6,13,7,9};
    }

    void SparseMatrix::AllocMemory(int32_t N, int32_t size)
    {
        this->N = N;
        this->size = size;
        ggu.resize(size);
        ggl.resize(size);
        di.resize(N);
        ig.resize(N + 1);
        jg.resize(size);
    }

    void SparseMatrix::InsertElem(double elem, int32_t i, int32_t j)
    {
        if (i == j)
        {
            di[i] = elem;
        }

        if (i < j)
        {
            int ind;
            int igj0 = ig[j];
            int igj1 = ig[j + 1];
            for (ind = igj0; ind < igj1; ind++)
                if (jg[ind] == i)
                    ggu[ind] = elem;
        }
        else
        {
            int ind;
            int igi0 = ig[i];
            int igi1 = ig[i + 1];
            for (ind = igi0; ind < igi1; ind++)
                if (jg[ind] == j)
                    ggl[ind] = elem;
        }
    }

    double SparseMatrix::GetElement(int32_t i, int32_t j)
    {

        if (i == j)
        {
            return di[i];
        }

        if (i < j)
        {
            int ind;
            int igj0 = ig[j];
            int igj1 = ig[j + 1];
            for (ind = igj0; ind < igj1; ind++)
                if (jg[ind] == i)
                    return ggu[ind];
        }
        else
        {
            int ind;
            int igi0 = ig[i];
            int igi1 = ig[i + 1];
            for (ind = igi0; ind < igi1; ind++)
                if (jg[ind] == j)
                    return ggl[ind];
        }

        return 0;
    }

    void SparseMatrix::PrintDenseMatrix(const std::string &filename, int32_t precession)
    {
        std::ofstream out(filename);
        out << "Dense Sparse Matrix\n";
        for (int32_t i = 0; i < N; i++)
        {

            for (int32_t j = 0; j < N; j++)
            {
                double elem = GetElement(i, j);
                out << std::left << std::setprecision(precession) << elem << " ";
            }
            out << "\n";
        }
        out.close();
    }

    ProfileMatrix::ProfileMatrix(int32_t N, int32_t size)
    {
        this->N = N;
        this->size = size;

        ia.resize(N + 1);
        au.resize(size);
        al.resize(size);
        di.resize(N);

        // di = {1, 5, 8, 12, 15, 19};
        // ia = {0, 0, 1, 1, 4, 5, 9};
        // al = {4, 10, 11, 0, 14, 16, 17, 0, 0};
        // au = {2, 3, 6, 0, 13, 7, 9, 0, 0};
    }

    void ProfileMatrix::AllocateMemory(int32_t N, int32_t size)
    {
        this->N = N;
        this->size = size;

        ia.resize(N + 1);
        au.resize(size);
        al.resize(size);
        di.resize(N);
    }

    void ProfileMatrix::InsertElem(double elem, int32_t i, int32_t j)
    {
        if (i == j)
        {
            di[i] = elem;
        }

        if (i < j)
        {
            int32_t ind = ia[j + 1] - (j - i);
            if (ind >= ia[j] && ind < ia[j + 1])
                au[ind] = elem;
        }
        else
        {
            int32_t ind = ia[i + 1] - (i - j);
            if (ind >= ia[i] && ind < ia[i + 1])
                al[ind] = elem;
        }
    }

    double ProfileMatrix::GetElement(int32_t i, int32_t j)
    {
        if (i == j)
        {
            return di[i];
        }

        if (i < j)
        {
            int32_t ind = ia[j + 1] - (j - i);
            if (ind >= ia[j] && ind < ia[j + 1])
                return au[ind];
        }
        else
        {
            int32_t ind = ia[i + 1] - (i - j);
            if (ind >= ia[i] && ind < ia[i + 1])
                return al[ind];
        }

        return 0;
    }

    void ProfileMatrix::PrintDenseMatrix(const std::string &filename, int32_t precession)
    {
        std::ofstream out(filename);
        out << "Dense Sparse Matrix\n";
        for (int32_t i = 0; i < N; i++)
        {

            for (int32_t j = 0; j < N; j++)
            {
                double elem = GetElement(i, j);
                out << std::left << std::setprecision(precession) << elem << " ";
            }
            out << "\n";
        }
        out.close();
    }

    void SparseMatrix2ProfileMatrix(SparseMatrix &src_Matr, ProfileMatrix &dst_Matr)
    {
        /* В цикле по элементам */
        dst_Matr.N = src_Matr.N;

        dst_Matr.ia.resize(src_Matr.N + 1);
        dst_Matr.ia[0] = 0;
        // dst_Matr.ia[1] = 0;

        for (int32_t i = 0; i < src_Matr.N; i++)
        {
            int32_t cnt_elem = 0;
            for (int32_t j = 0; j < i; j++)
            {
                double elem = src_Matr.GetElement(i, j);
                if (std::abs(elem) >= 0)
                {
                    cnt_elem = i - j;
                    break;
                }
            }
            dst_Matr.ia[i + 1] = dst_Matr.ia[i] + cnt_elem;
        }

        // Выделение памяти
        dst_Matr.size = dst_Matr.ia[dst_Matr.N];
        dst_Matr.au.resize(dst_Matr.size);
        dst_Matr.al.resize(dst_Matr.size);
        dst_Matr.di.resize(dst_Matr.N);

        for (int32_t i = 0; i < src_Matr.N; i++)
        {
            int32_t cnt_elem = 0;
            for (int32_t j = 0; j < src_Matr.N; j++)
            {
                double elem = src_Matr.GetElement(i, j);

                if (std::abs(elem) > 0)
                    dst_Matr.InsertElem(elem, i, j);
            }
        }
    }

};
