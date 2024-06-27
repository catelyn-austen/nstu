#include <iostream>
#include <fstream>

class FEM
{
public:
    // инициализация локальных векторов и матриц
    double G[3][3]{ };
    double M[3][3]{ };
    double local_A[3][3]{ };
    double local_F[3]{ };
    double board_F[2]{ };
    double board_A[2][2]{ };

    // переменные и матрицы для хранения информации о конечных элементах
    int points_count = 0;
    int finit_el_count = 0;
    int areas_count = 0;
    int border_cond_count = 0;
    double** points = NULL;
    int** finits = NULL;
    int** border_cond = NULL;

    // инициализация векторов для разложения глобальной матрицы А и работы с ней
    int* ig;
    int* jg;
    double* ggl;
    double* ggu;
    double* di;
    // глобальный вектор правой части
    double* F;

    // векторы для решения СЛАУ
    double* LUl;
    double* LUu;
    double* LUd;
    double* r, * z, * p;
    double* Ar, * Ax;
    double* tmp1, * tmp2;
    double* q;

    // функции для вычисления локальных матриц и векторов и для занесения их в глобальную А и F
    double Det(double* p1, double* p2, double* p3);
    double mesG(double* p1, double* p2);
    int G_matrix(double* p1, double* p2, double* p3, int area);
    int M_matrix(double* p1, double* p2, double* p3, int area);
    int local_matrix(int num_of_triangle);
    int cond_23(int* current_border, double beta);
    int cond_1(int* current_border);
    int Global_matrix();

    // функции для решения СЛАУ
    int MultMatrByVec(double* vec, double* res);
    void doForwardLU(double* b, double* res, bool flag);
    void doBackwardLU(double* b, double* res, bool flag);
    int ComputeLU();
    int LosLU();
    int Los();
    int AX();
    void PrintQ();

    // инициализация и деструктор
    FEM(std::string points_file, std::string triangles_file, std::string border_cond_file, std::string areas);
    FEM();
    ~FEM();

    int Init(std::string coords_file, std::string finits_file, std::string border_file, std::string areas);

    int Portrait();

    // функции для вычисления коэф-тов
    double f(double* point, int area);
    double coef_lambda(double* point, int area);
    double coef_gamma(int area);
    double coef_beta(int area);
};

double f_cond1(double* x, int k)
{
    switch (k)
    {
        case 0: return x[0];
            break;
        /*case 1: return 6 * x[0] + 5;
            break;*/
        /*case 2: return 6 * x[0] + 5;
            break;
        case 3: return x[1];*/
            break;
    }
    return 0;
}
double f_cond2(double* x, int k)
{
    switch (k)
    {
        case 0: return x[1] + 1;
            break;
        case 1: return - x[1] - 1;
            break;
        /*case 2: return 2.;
            break;*/
    }
    return 0;
}
double f_cond3(double* x, int k)
{
    switch (k)
    {
        case 0: return x[0] * x[1] + 2 * x[0] + x[1] + 1;
            break;
    }
    return 0;
}

FEM::FEM(std::string coords, std::string finits, std::string board, std::string areas)
{
    Init(coords, finits, board, areas);
}

FEM::FEM()
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            local_A[i][j] = 0;
            M[i][j] = 0;
            G[i][j] = 0;
        }
        local_F[i] = 0;
    }

    points_count = 0;
    finit_el_count = 0;
    areas_count = 0;
    border_cond_count = 0;

    points = NULL;
    finits = NULL;
    border_cond = NULL;

    F = NULL;
    di = NULL;
    ggl = NULL;
    ggu = NULL;
    ig = NULL;
    jg = NULL;

    LUl = NULL;
    LUu = NULL;
    LUd = NULL;
    q = NULL;

    r = z = p = NULL;
    Ar = Ax = NULL;
    tmp1 = tmp2 = NULL;
}

FEM::~FEM()
{
    for (int i = 0; i < points_count; i++)
        delete[] points[i];
    delete[] points;

    for (int i = 0; i < finit_el_count; i++)
        delete[] finits[i];
    delete[] finits;

    for (int i = 0; i < border_cond_count; i++)
        delete[] border_cond[i];
    delete[] border_cond;
}

int FEM::Init(std::string points_file, std::string triangles_file, std::string border_cond_file, std::string areas)
{
    std::ifstream in;

    // чтение из файла с координатами точек
    in.open(points_file);
    if (!in.is_open())
    {
        std::cout << points_file << " file opening error" << std::endl;
        return 1;
    }
    in >> points_count;
    points = new double* [points_count];
    for (int i = 0; i < points_count; i++)
    {
        points[i] = new double[2];
        in >> points[i][0] >> points[i][1];
    }
    in.close();

    // чтение из файла с конечными элементами
    in.open(triangles_file);
    if (!in.is_open())
    {
        std::cout << triangles_file << " file opening error" << std::endl;
        return 1;
    }
    in >> finit_el_count;
    finits = new int* [finit_el_count];
    for (int i = 0; i < finit_el_count; i++)
    {
        finits[i] = new int[4];
        in >> finits[i][0] >> finits[i][1] >> finits[i][2] >> finits[i][3];
    }
    in.close();

    // чтение из файла с краевыми условиями 
    in.open(border_cond_file);
    if (!in.is_open())
    {
        std::cout << border_cond_file << " file opening error" << std::endl;
        return 1;
    }
    in >> border_cond_count;
    border_cond = new int* [border_cond_count];
    for (int i = 0; i < border_cond_count; i++)
    {
        border_cond[i] = new int[5];
        in >> border_cond[i][0] >> border_cond[i][1] >> border_cond[i][2] >> border_cond[i][3] >> border_cond[i][4];
    }
    in.close();

    return 0;
}

void sort(int* arr, int k)
{
    int l = 0;
    for (int i = k - 1; i >= 0; i--)
        for (int j = 0; j <= i; j++)
            if (arr[j] > arr[j + 1])
            {
                l = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = l;
            }
}

int FEM::Portrait()
{
    int k = 0, kk = 0;
    int key = 0;
    int pos = 0;
    int* arr = new int[points_count];
    struct List
    {
        int num;
        List* next;
    };
    List* list = new List[points_count];
    List* p = NULL;
    ig = new int[points_count + 1];
    for (int i = 0; i < points_count; i++)
        list[i].next = NULL;

    for (int i = 0; i < finit_el_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            key = 0;
            k = finits[i][(j == 2)];
            kk = finits[i][((j != 0) + 1)];
            if (k < kk)
            {
                k += kk;
                kk = k - kk;
                k -= kk;
            }
            p = &list[k];
            while (p->next)
            {
                if (p->next->num == kk)
                {
                    key = 1;
                    break;
                }
                p = p->next;
            }
            if (!key)
            {
                p->next = new List;
                p->next->num = kk;
                p->next->next = NULL;
            }
        }
    }

    ig[0] = 0;
    for (int i = 0; i < points_count; i++) {
        k = 0;
        p = &list[i];
        while (p = p->next)
            k++;
        ig[i + 1] = ig[i] + k;
    }

    F = new double[points_count] {};
    di = new double[points_count] {};
    q = new double[points_count] {};
    jg = new int[ig[points_count] - 1]{ };
    ggu = new double[ig[points_count]]{};
    ggl = new double[ig[points_count]]{};

    for (int i = 0; i < points_count; i++) {
        k = 0;
        key = 0;
        p = &list[i];
        while (p = p->next) {
            arr[k] = p->num;
            k++;
            key = 1;
        }
        if (key) {
            sort(arr, --k);
            int ii = 0;
            int jj = 0;
            for (ii = pos, jj = 0; ii <= k + pos; ii++, jj++)
                jg[ii] = arr[jj];
            pos += k + 1;
        }
    }
    return 0;
}
double FEM::Det(double* p1, double* p2, double* p3)
{
    return (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]);
}
double FEM::mesG(double* p1, double* p2)
{
    return (sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1])));
}
int FEM::M_matrix(double* p1, double* p2, double* p3, int area)
{
    double coef = std::abs(Det(p1, p2, p3)) / 24.;
    double* fnew = new double[3];

    fnew[0] = coef * f(p1, area);
    fnew[1] = coef * f(p2, area);
    fnew[2] = coef * f(p3, area);

    local_F[0] = 2 * fnew[0] + fnew[1] + fnew[2];
    local_F[1] = fnew[0] + 2 * fnew[1] + fnew[2];
    local_F[2] = fnew[0] + fnew[1] + 2 * fnew[2];

    coef *= coef_gamma(area);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (j != i)
                M[i][j] = coef;
            else
                M[i][j] = 2 * coef;
    return 1;
}
int FEM::G_matrix(double* p1, double* p2, double* p3, int area)
{
    double det = Det(p1, p2, p3);
    double aver_lambda = (coef_lambda(p1, area) + coef_lambda(p2, area) + coef_lambda(p3, area)) / 3.;
    double coef = aver_lambda * abs(det) / 2.;

    G[0][0] = coef * ((p2[1] - p3[1]) * (p2[1] - p3[1]) / (det * det) + (p3[0] - p2[0]) * (p3[0] -
        p2[0]) / (det * det));
    G[0][1] = coef * ((p2[1] - p3[1]) * (p3[1] - p1[1]) / (det * det) + (p3[0] - p2[0]) * (p1[0] -
        p3[0]) / (det * det));
    G[0][2] = coef * ((p2[1] - p3[1]) * (p1[1] - p2[1]) / (det * det) + (p3[0] - p2[0]) * (p2[0] -
        p1[0]) / (det * det));
    G[1][0] = coef * ((p3[1] - p1[1]) * (p2[1] - p3[1]) / (det * det) + (p1[0] - p3[0]) * (p3[0] -
        p2[0]) / (det * det));
    G[1][1] = coef * ((p3[1] - p1[1]) * (p3[1] - p1[1]) / (det * det) + (p1[0] - p3[0]) * (p1[0] -
        p3[0]) / (det * det));
    G[1][2] = coef * ((p3[1] - p1[1]) * (p1[1] - p2[1]) / (det * det) + (p1[0] - p3[0]) * (p2[0] -
        p1[0]) / (det * det));
    G[2][0] = coef * ((p2[1] - p3[1]) * (p1[1] - p2[1]) / (det * det) + (p3[0] - p2[0]) * (p2[0] -
        p1[0]) / (det * det));
    G[2][1] = coef * ((p3[1] - p1[1]) * (p1[1] - p2[1]) / (det * det) + (p1[0] - p3[0]) * (p2[0] -
        p1[0]) / (det * det));
    G[2][2] = coef * ((p1[1] - p2[1]) * (p1[1] - p2[1]) / (det * det) + (p2[0] - p1[0]) * (p2[0] -
        p1[0]) / (det * det));

    return 1;
}
int FEM::local_matrix(int num_of_triangle)
{
    int a = finits[num_of_triangle][0];
    int b = finits[num_of_triangle][1];
    int c = finits[num_of_triangle][2];
    int area = finits[num_of_triangle][3];

    M_matrix(points[a], points[b], points[c], area);
    G_matrix(points[a], points[b], points[c], area);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            local_A[i][j] = G[i][j] + M[i][j];
    return 0;
}

int FEM::cond_23(int* current_border, double beta)
{
    double coef = 0;
    if (current_border[2] == 1)
        return 0;
    else if (current_border[2] == 2)
    {
        coef = mesG(points[current_border[0]], points[current_border[1]]) / 6.;
        double teta1 = f_cond2(points[current_border[0]], current_border[3]);
        double teta2 = f_cond2(points[current_border[1]], current_border[3]);

        board_F[0] = coef * (2 * teta1 + teta2);
        board_F[1] = coef * (teta1 + 2 * teta2);
        return 1;
    }
    else if (current_border[2] == 3)
    {
        coef = beta * mesG(points[current_border[0]], points[current_border[1]]) / 6.;
        double ubeta1 = f_cond3(points[current_border[0]], current_border[3]);
        double ubeta2 = f_cond3(points[current_border[1]], current_border[3]);
        board_A[0][0] = board_A[1][1] = 2 * coef;
        board_A[0][1] = board_A[1][0] = coef;
        board_F[0] = coef * (2 * ubeta1 + ubeta2);
        board_F[1] = coef * (ubeta1 + 2 * ubeta2);
        return 2;
    }
    return 0;
}
int FEM::cond_1(int* current_border)
{
    int kol = 0, m = 0;
    int i0, i1;

    di[current_border[0]] = 1;
    di[current_border[1]] = 1;
    F[current_border[0]] = f_cond1(points[current_border[0]], current_border[3]);
    F[current_border[1]] = f_cond1(points[current_border[1]], current_border[3]);
    kol = ig[current_border[0] + 1] - ig[current_border[0]];

    for (int i = 0; i < kol; i++)
        ggl[ig[current_border[0]] + i] = 0;
    kol = ig[current_border[1] + 1] - ig[current_border[1]];

    for (int i = 0; i < kol; i++)
        ggl[ig[current_border[1]] + i] = 0;

    for (int i = current_border[0] + 1; i < points_count; i++)
    {
        i0 = ig[i];
        i1 = ig[i + 1];
        for (int p = i0; p < i1; p++)
            if (jg[p] == current_border[0]) {
                ggu[p] = 0;
                continue;
            }
    }

    for (int i = current_border[1] + 1; i < points_count; i++)
    {
        i0 = ig[i];
        i1 = ig[i + 1];
        for (int p = i0; p < i1; p++)
            if (jg[p] == current_border[1])
            {
                ggu[p] = 0;
                continue;
            }
    }
    return 0;
}

int FEM::Global_matrix()
{
    int i = 0;
    int j = 0;
    int k = 0;
    int p = 0;
    int ibeg = 0;
    int iend = 0;
    int h = 0;
    int key = 0;
    int kol = 0;
    int* L = new int[3];
    int* L2 = new int[2];
    int* K = new int[border_cond_count];

    Portrait();

    for (k = 0; k < finit_el_count; k++)
    {
        local_matrix(k);
        memcpy(L, finits[k], 3 * sizeof(int));
        for (i = 0; i < 3; i++)
        {
            ibeg = L[i];
            for (j = i + 1; j < 3; j++)
            {
                iend = L[j];
                if (ibeg < iend)
                {
                    h = ig[iend];
                    while (jg[h++] - ibeg);
                    h--;
                    ggl[h] += local_A[i][j];
                    ggu[h] += local_A[j][i];
                }
                else
                {
                    h = ig[ibeg];
                    while (jg[h++] - iend);
                    h--;
                    ggl[h] += local_A[i][j];
                    ggu[h] += local_A[j][i];
                }
            }
            di[ibeg] += local_A[i][i];
        }
        for (i = 0; i < 3; i++)
            F[L[i]] += local_F[i];
    }

    for (k = 0; k < border_cond_count; k++) {
        key = cond_23(border_cond[k], coef_beta(border_cond[k][4]));
        if (!key) {
            K[p] = k;
            p++;
            continue;
        }
        L2[0] = border_cond[k][0];
        L2[1] = border_cond[k][1];
        if (key == 2) {
            for (i = 0; i < 2; i++) {
                ibeg = L2[i];
                for (j = i + 1; j < 2; j++) {
                    iend = L2[j];
                    if (ibeg < iend) {
                        h = ig[iend];
                        while (jg[h++] - ibeg);
                        h--;
                        ggl[h] += board_A[i][j];
                        ggu[h] += board_A[j][i];
                    }
                    else {
                        h = ig[ibeg];
                        while (jg[h++] - iend);
                        h--;
                        ggl[h] += board_A[i][j];
                        ggu[h] += board_A[j][i];
                    }
                }
                di[ibeg] += board_A[i][i];
            }
        }
        for (i = 0; i < 2; i++)
            F[L2[i]] += board_F[i];
        for (i = 0; i < 2; i++) {
            board_F[i] = 0;
            for (j = 0; j < 2; j++)
                board_A[i][j] = 0;
        }
    }
    for (i = 0; i < p; i++)
        cond_1(border_cond[K[i]]);
    return 0;
}

void NullVec(double* vec, int size)
{
    for (int i = 0; i < size; i++)
        vec[i] = 0;
}
void SubVec(double* vec1, double* vec2, int size, double* res)
{
    for (int i = 0; i < size; i++)
        res[i] = vec1[i] - vec2[i];
}
void SumVec(double* vec1, double* vec2, int size, double* res)
{
    for (int i = 0; i < size; i++)
        res[i] = vec1[i] + vec2[i];
}
double MultScalar(double* vec1, double* vec2, int size)
{
    double res = 0;
    for (int i = 0; i < size; i++)
        res += vec1[i] * vec2[i];
    return res;
}
void MultScalar(double a, double* vec, int size, double* res)
{
    NullVec(res, size);
    for (int i = 0; i < size; i++)
        res[i] = a * vec[i];
}
int FEM::MultMatrByVec(double* vec, double* res)
{
    long st;
    for (int i = 0; i < points_count; i++) {
        res[i] = di[i] * vec[i];
        for (int j = ig[i]; j < ig[i + 1]; j++) {
            st = jg[j];
            res[i] += ggl[j] * vec[st];
            res[st] += ggu[j] * vec[i];
        }
    }
    return 0;
}
void FEM::doForwardLU(double* b, double* res, bool flag)
{
    if (!flag)
    {
        for (int i0, i1, i = 0; i < points_count; i++)
        {
            double sum = 0;
            i0 = ig[i];
            i1 = ig[i + 1];
            for (int k = i0; k < i1; k++)
            {
                int j = jg[k];
                sum += LUl[k] * res[j];
            }
            res[i] = (b[i] - sum) / LUd[i];
        }
    }
    else
        for (int i0, i1, i = 0; i < points_count; i++)
        {
            double sum = 0;
            i0 = ig[i];
            i1 = ig[i + 1];
            for (int k = i0; k < i1; k++)
            {
                int j = jg[k];
                sum += LUu[k] * res[j];
            }
            res[i] = (b[i] - sum);
        }
}
void FEM::doBackwardLU(double* b, double* res, bool flag)
{
    if (!flag)
    {
        for (int i0, i1, i = points_count - 1; i >= 0; i--)
        {
            i0 = ig[i];
            i1 = ig[i + 1];
            res[i] = b[i];
            for (int j, k = i0; k < i1; k++)
            {
                j = jg[k];
                b[j] -= LUu[k] * res[i];
            }
        }
    }
    else
    {
        for (int i0, i1, i = points_count - 1; i >= 0; i--)
        {
            i0 = ig[i];
            i1 = ig[i + 1];
            res[i] = b[i] / LUd[i];
            for (int j, k = i0; k < i1; k++)
            {
                j = jg[k];
                b[j] -= LUl[k] * res[i];
            }
        }
    }
}
int FEM::ComputeLU()
{
    LUl = new double[ig[points_count]]{ };
    LUu = new double[ig[points_count]]{ };
    LUd = new double[points_count] {};

    for (int i = 0; i < points_count; i++)
    {
        int i0 = ig[i];
        int i1 = ig[i + 1];
        double sumD = 0;

        for (int k = i0; k < i1; k++)
        {
            double sumL = 0, sumU = 0;
            int j = jg[k];

            int j0 = ig[j];
            int j1 = ig[j + 1];

            int ik = i0;
            int kj = j0;

            for (; ik < k and kj < j1; )
            {
                if (jg[ik] < jg[kj])
                    ik++;
                if (jg[ik] > jg[kj])
                    kj++;
                if (jg[ik] == jg[kj]) {
                    sumL += LUl[ik] * LUu[kj];
                    sumU += LUl[kj] * LUu[ik];
                    ik++; kj++;
                }

            }
            LUl[k] = ggl[k] - sumL;
            LUu[k] = (ggu[k] - sumU) / LUd[j];
            sumD += LUl[k] * LUu[k];
        }
        LUd[i] = di[i] - sumD;
    }
    return 0;
}
int FEM::LosLU()
{
    r = new double[points_count] {};
    z = new double[points_count] {};
    p = new double[points_count] {};
    Ar = new double[points_count] {};
    Ax = new double[points_count] {};
    tmp1 = new double[points_count] {};
    tmp2 = new double[points_count] {};

    ComputeLU();

    double residual = 0;
    int n = points_count;
    int maxiter = 20000;
    double eps = 1e-15;
    double a, b;

    MultMatrByVec(q, Ax);
    SubVec(F, Ax, points_count, tmp1);
    doForwardLU(tmp1, r, 0);

    memcpy(tmp1, r, n * sizeof(double));
    doBackwardLU(tmp1, z, 0);

    MultMatrByVec(z, tmp1);    //// Az = p
    doForwardLU(tmp1, p, 0);

    residual = MultScalar(r, r, n);

    int k;
    for (k = 0; k < maxiter && residual > eps; k++)
    {
        double pp0 = MultScalar(p, p, n);

        a = MultScalar(p, r, n) / pp0;

        for (int i = 0; i < n; i++)
            q[i] = q[i] + a * z[i];

        for (int i = 0; i < n; i++)
            r[i] = r[i] - a * p[i];

        doBackwardLU(r, tmp1, 0);
        MultMatrByVec(tmp1, tmp2);
        doForwardLU(tmp2, Ar, 0);
        b = -MultScalar(p, Ar, n) / pp0;


        doBackwardLU(r, tmp1, 0);
        for (int i = 0; i < n; i++)
            z[i] = tmp1[i] + b * z[i];

        for (int i = 0; i < n; i++)
            p[i] = Ar[i] + b * p[i];

        residual = MultScalar(r, r, n);
        //std::cout << "Iteration " << k << " Residual norm in square: " << residual << std::endl;
    }
    //std::cout << "Iteration " << k - 1 << " Residual norm in square: " << residual << std::endl;
    return 0;
}

int FEM::Los()
{
    r = new double[points_count];
    z = new double[points_count];
    p = new double[points_count];
    Ar = new double[points_count];
    Ax = new double[points_count];
    double residual = 0;
    int n = points_count;
    int maxiter = 20000;
    double eps = 1e-15;
    double a, b;

    MultMatrByVec(q, Ax);
    SubVec(F, Ax, points_count, r);
    memcpy(z, r, n * sizeof(double));
    MultMatrByVec(z, p);    //// Az = p

    residual = MultScalar(r, r, n);

    int k;
    for (k = 0; k < maxiter && residual > eps; k++)
    {
        double pp0 = MultScalar(p, p, n);

        a = MultScalar(p, r, n) / pp0;

        for (int i = 0; i < n; i++)
            q[i] = q[i] + a * z[i];

        for (int i = 0; i < n; i++)
            r[i] = r[i] - a * p[i];

        MultMatrByVec(r, Ar);
        b = -MultScalar(p, Ar, n) / pp0;

        for (int i = 0; i < n; i++)
            z[i] = r[i] + b * z[i];

        for (int i = 0; i < n; i++)
            p[i] = Ar[i] + b * p[i];

        residual = MultScalar(r, r, n);
        //std::cout << "Iteration " << k << " Residual norm in square: " << residual << std::endl;
    }
    //std::cout << "Iteration " << k - 1 << " Residual norm in square: " << residual << std::endl;
    return 0;
}

int FEM::AX()
{
    double* vectorq = new double[5];
    vectorq[0] = 6.;
    vectorq[1] = 8.;
    vectorq[2] = 4.;
    vectorq[3] = 6.;
    vectorq[4] = 6.;
    double* vectorb = new double[5]{};
    MultMatrByVec(vectorq, vectorb);
    for (int i = 0; i < points_count; i++)
        std::cout << vectorb[i] << std::endl;
    return 0;

}

double FEM::f(double* x, int area)
{
    switch (area)
    {
        case 0: return 0;
            break;
        case 1: return 0;
            break;
        default: return 0;
    }
}
double FEM::coef_lambda(double* point, int area)
{
    switch (area)
    {
        case 0: return 1.;
            break;
        case 1: return 0.;
            break;
        default: return 0;
    }
}
double FEM::coef_gamma(int area)
{
    switch (area)
    {
        case 0: return 0.;
            break;
        case 1: return 1.;
            break;
        default: return 0;
    }
}
double FEM::coef_beta(int area)
{
    switch (area)
    {
        case 0: return 1.;

        default: return 0;
    }
}

void FEM::PrintQ()
{
    for (int i = 0; i < points_count; i++)
        printf_s("%.14lf\n", q[i]);
}
int main()
{
    FEM fem("coords.txt", "finit_elements.txt", "board.txt", "areas.txt");

    fem.Global_matrix();
    fem.Los();
    //fem.AX();
    fem.PrintQ();
}
