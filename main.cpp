#include <vector>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <cstdlib>

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

    void dgetrs_(char* C, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
    void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N,
                const int *K, double *ALPHA, double *A, const int *LDA, double *B,
                const int *LDB, double *BETA, double *C, const int *LDC);

}

void inverse(double* A, int N)
{
    int *IPIV = new int[N];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete[] IPIV;
    delete[] WORK;
}

//~ void mul_A_B_in_C(const std::vector<double> &A, const int A_m, const int A_n,
              //~ const std::vector<double> &B, const int B_m, const int B_n,
              //~ std::vector<double> &C, const int C_m, const int C_n) {
void mul_A_B_in_C(std::vector<double> &A, const int A_m, const int A_n,
                  std::vector<double> &B, const int B_m, const int B_n,
                  std::vector<double> &C, const int C_m, const int C_n) {
    const char N = 'N';
    double alpha = 1, beta = 1;
    C.assign(C_m*C_n, 0);
    dgemm_(&N, &N, &C_m, &C_n, &A_n, &alpha, A.data(), &A_m, B.data(),
           &B_m, &beta, C.data(), &C_m);
}

void full_A(std::vector<double> &A, const int m, const int n, double h_x, double h_y) {
    A.clear();
    A.resize(m*n*m*n, 0);
    double h_x__2 = 1/(h_x*h_x), h_y__2 = 1/(h_y*h_y);
    for(auto i = 0; i < m*n; ++i) {
        A[i*m*n+i] = -2*h_x__2-2*h_y__2;

        if(i%n) {
            A[i*m*n+i-1] = 1*h_x__2;
        }

        if((i+1)%n) {
            A[i*m*n+i+1] = 1*h_x__2;
        }

        if(i >= n) {
            A[i*m*n+i-n] = 1*h_y__2;
        }

        if(i < m*n - n) {
            A[i*m*n+i+n] = 1*h_y__2;
        }
    }
}

void mul_diag_A(std::vector<double> &a, std::vector<double> &A, int m, int n){

    for(auto i = 0; i < m*n; ++i) {
        A[i*m*n+i] *= a[i];

        if(i%n) {
            A[i*m*n+i-1] *= a[i-1];
        }

        if((i+1)%n) {
            A[i*m*n+i+1] *= a[i+1];
        }

        if(i >= n) {
            A[i*m*n+i-n] *= a[i-n];
        }

        if(i < m*n - n) {
            A[i*m*n+i+n] *= a[i+n];
        }
    }
}
void print_A(std::ostream &out, const std::vector<double> &A,const int m, const int n) {

    for(auto i = 0; i < m*n; ++i) {
        for(auto j = 0; j < m*n; ++j) {
            out << A[i*m*n+j] << " ";
        }
        out<<std::endl;
    }
}

void E_plus_A_in_A(std::vector<double> &A, int m, int n){

    for(auto i = 0; i < m*n; ++i) {
        A[i*m*n+i] = 1+A[i*m*n+i];

        //~ if(i%n) {
            //~ A[i*m*n+i-1] *= -1;
        //~ }

        //~ if((i+1)%n) {
            //~ A[i*m*n+i+1] *= -1;
        //~ }

        //~ if(i >= n) {
            //~ A[i*m*n+i-n] *= -1;
        //~ }

        //~ if(i < m*n - n) {
            //~ A[i*m*n+i+n] *= -1;
        //~ }
    }
}

//~ int main(){

    //~ std::vector<double> A = {
        //~ 1,2,
        //~ 3,4
    //~ };
    //~ std::vector<double> B = {
        //~ 1,2,
        //~ 3,4
    //~ };

    //~ std::vector<double> C(4);

    //~ mul_A_B_in_C(A, 2, 2, B, 2, 2, C, 2, 2);

    //~ printf("%f %f\n", C[0], C[2]);
    //~ printf("%f %f\n", C[1], C[3]);

    //~ return 0;
//~ }

double function(double x, double y, double t) {
    //~ return 0;
    //~ return sin(x)*cos(y)+cos(t);
    return sin(x)*cos(y)+cos(t)+1/(t+1);
}

void full_f(std::vector<double> &f, double (&source_function)(double x, double y, double t),
            std::vector<double> &x_axis_grid, std::vector<double> &y_axis_grid, double t ) {
    f.clear();

    for(auto y: y_axis_grid) {
        for(auto x: x_axis_grid) {
            f.push_back(source_function(x, y, t));
        }
    }
}

void full_axis_grid(std::vector<double> &axis_grid, double start, double finish, int number_of_nodes) {
    axis_grid.clear();
    axis_grid.resize(number_of_nodes);
    axis_grid.front() = start;
    double interval = (finish - start)/(number_of_nodes-1);
    for(auto it = axis_grid.begin()+1; it != axis_grid.end()-1; it++ ) {
        *it = *(it-1)+interval;
    }
    axis_grid.back() = finish;
}

void mul_vector_scalar(std::vector<double> &vector, double scalar) {
    for(auto &element: vector) {
        element *= scalar;
    }
}

double bound_1(double x, double t) {
    //~ return 0;
    //~ return 1;
    //~ return sin(x)*cos(t);
    return sin(x)-cos(t);
}
double bound_2(double x, double t) {
    //~ return 0;
    //~ return 1;
    return sin(x)*cos(t);
}
double bound_3(double x, double t) {
    //~ return 0;
    //~ return 1;
    //~ return sin(x)*cos(t);
    return sin(x)+cos(t);
}
double bound_4(double x, double t) {
    //~ return 0;
    //~ return 1;
    //~ return sin(x)*cos(t);
    return sin(x)*1/(t+1);
}

void vector_plus_vector(std::vector<double> &v1,const std::vector<double> &v2) {
    for(int i = 0; i < v1.size(); ++i) {
        v1[i] += v2[i];
    }
}

template <class T>
void print_vector(std::ostream &out, const std::vector<T> &vector) {
    for (auto element: vector) {
        out << element << " ";
    }
    out << std::endl;
}

void boundary_conditions(std::vector<double> &b,
                         std::vector<double> &a,
                         double tau, double h_x, double h_y,
                         double (&left_bound)(double, double),
                         double (&up_bound)(double, double),
                         double (&right_bound)(double, double),
                         double (&down_bound)(double, double),
                         std::vector<double> &x_axis_grid,
                         std::vector<double> &y_axis_grid,
                         double t) {
    for(auto i = 0; i < x_axis_grid.size()-2; ++i) {
        b[i] += down_bound(x_axis_grid[i+1], t+tau) * tau /(h_y*h_y) * a[i];
        b[b.size()-i-1] += up_bound(x_axis_grid[x_axis_grid.size()-i-2], t+tau) * tau /(h_y*h_y)* a[b.size()-i-1];
    }
    for(auto i = 1; i < y_axis_grid.size()-2; ++i) {
        b[i*(x_axis_grid.size()-2)] += left_bound(y_axis_grid[i], t+tau) * tau /(h_x*h_x) * a[i*(x_axis_grid.size()-2)];
        b[(i+1)*(x_axis_grid.size()-2)-1] += right_bound(y_axis_grid[i], t+tau) * tau /(h_x*h_x) * a[(i+1)*(x_axis_grid.size()-2)-1];
    }

}

void calculate_coefficients(std::vector<double> &A,
                            std::vector<double> &x,
                            std::vector<double> &b,
                            std::vector<double> &left_x,
                            std::vector<double> &up_x,
                            std::vector<double> &right_x,
                            std::vector<double> &down_x,
                            double h_x, double h_y, double tau) {

    std::vector<double> b_1;
    auto M = x.size();
    mul_A_B_in_C(A, M, M, x, M, 1, b_1, M, 1);

    for(auto i = 0; i < down_x.size()-2; ++i) {
        b_1[i] += down_x[i+1] /(h_y*h_y);
        b_1[b_1.size()-i-1] += up_x[up_x.size()-i-2] / (h_y*h_y);
    }
    for(auto i = 1; i < left_x.size()-2; ++i) {
        b_1[i*(down_x.size()-2)] += left_x[i] / (h_x*h_x);
        b_1[(i+1)*(down_x.size()-2)-1] += right_x[i] /(h_x*h_x);
    }
    //~ for(auto i = 0; i < x_axis_grid.size()-2; ++i) {
        //~ b[i] += down_bound(x_axis_grid[i+1], t+tau) * tau /(h_y*h_y) * a[i];
        //~ b[b.size()-i-1] += up_bound(x_axis_grid[x_axis_grid.size()-i-2], t+tau) * tau /(h_y*h_y)* a[b.size()-i-1];
    //~ }
    //~ for(auto i = 1; i < y_axis_grid.size()-2; ++i) {
        //~ b[i*(x_axis_grid.size()-2)] += left_bound(y_axis_grid[i], t+tau) * tau /(h_x*h_x) * a[i*(x_axis_grid.size()-2)];
        //~ b[(i+1)*(x_axis_grid.size()-2)-1] += right_bound(y_axis_grid[i], t+tau) * tau /(h_x*h_x) * a[(i+1)*(x_axis_grid.size()-2)-1];
    //~ }

    mul_vector_scalar(b, -1);
    vector_plus_vector(x, b);
    mul_vector_scalar(x, 1./tau);

    for(int i = 0; i < b_1.size(); ++i) {
        x[i] /= b_1[i];
    }
}

void full_bound(std::vector<double> &x, std::vector<double> &x_axis_grid,
                double (&bound_function)(double, double), double t) {
    x.clear();
    for(auto elem: x_axis_grid) {
        x.push_back(bound_function(elem, t));
    }
}

void full_coef_vector(std::vector<double> &a, double (&function)(double x, double y, double t),
                      std::vector<double> &x_axis_grid, std::vector<double> &y_axis_grid) {

    //~ for(auto &element: a) {
        //~ element = 1;
    //~ }
    //~ for(auto &element: a) {
        //~ element = 2.3;
    //~ }
    //~ for(auto i = 0; i < a.size(); ++i)  {
        //~ a[i] = i+1;
    //~ }
    full_f(a, function, x_axis_grid, y_axis_grid, 0);
}

double calc_euc_norm_2(const std::vector<double> &x) {
    double euc_norm_2 = 0;
    for(auto elem: x) {
        euc_norm_2 += elem*elem;
    }
    return euc_norm_2;
}

int args_parsing(int argc, char *argv[], std::vector <double> &args_values) {
    args_values.clear();
    for(int i = 1; i < argc && i <= 8; ++i) {
        args_values.push_back(atof(argv[i]));
    }

    const auto number_of_arguments = 8;
    if(args_values.size() < number_of_arguments) {
        std::cout << "Not enough command line arguments. 8 arguments are expected. " << std::endl;
        std::cout << "8 arguments are expected: x_min, x_max, y_min, y_max, m, k, tau, T. " << std::endl;
        return 1;
    }
    if(args_values[0] >= args_values[1] || args_values[2] >= args_values[3]) {
        std::cout << "Expected that second argument is greater than first, and fourth is greater than third." << std::endl;
        return 1;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    std::vector<double> A;
    std::vector<double> args_values;
    if(args_parsing(argc, argv, args_values)) {
        return 0;
    }

    double x_min = args_values[0], x_max = args_values[1], y_min = args_values[2], y_max = args_values[3];
    int m = args_values[4], n = args_values[5];
    double tau = args_values[6];
    double T = args_values[7];
    double h_x = (x_max-x_min)/(n-1), h_y = (y_max-y_min)/(m-1);
    full_A(A, m-2, n-2, h_x, h_y);
    std::vector<double> A_without_coef(A);
    int M = (m-2)*(n-2);

    std::vector<double> a;
    std::vector<double> x_axis_grid, y_axis_grid;
    full_axis_grid(x_axis_grid, x_min, x_max, n);
    full_axis_grid(y_axis_grid, y_min, y_max, m);

    a.resize(M,1);
    std::vector<double> y_axis_grid_in;
    std::vector<double> x_axis_grid_in;
    y_axis_grid_in.assign(y_axis_grid.begin()+1, y_axis_grid.end()-1);
    x_axis_grid_in.assign(x_axis_grid.begin()+1, x_axis_grid.end()-1);

    full_coef_vector(a, function, x_axis_grid_in, y_axis_grid_in);

    mul_diag_A(a, A, m-2, n-2);

    mul_vector_scalar(A, -tau);
    E_plus_A_in_A(A, m-2, n-2);

    std::vector<int> IPIV(M);
    int INFO = 0;

    // LU decomoposition of a general matrix
    dgetrf_(&M, &M, A.data(), &M, IPIV.data(), &INFO);

    std::vector<double> inverse_A(A);
    int LWORK = M*M;
    std::vector<double>WORK(LWORK);

    dgetri_(&M, inverse_A.data(), &M, IPIV.data(), WORK.data(), &LWORK, &INFO);

    std::vector<double> b(M,0);
    std::vector<double> c;
    char N = 'N';
    int NRHS = 1;
    double t = 0;

    std::vector<double> f;
    std::vector<double> x;
    std::vector<double> left_x, up_x, right_x, down_x;

    std::ofstream out;
    out.open("out.txt");
    out << x_min << std::endl << x_max << std::endl;
    out << y_min << std::endl << y_max << std::endl;
    out << m-2 << std::endl << n-2 << std::endl;

    for(t = 0; t < T; t += tau) {
        //vector b in Ax = b
        full_f(f, function, x_axis_grid, y_axis_grid, t+tau/3);
        mul_vector_scalar(f, tau);
        vector_plus_vector(b, f);
        mul_A_B_in_C(inverse_A, M, M, b, M, 1, c, M, 1);
        full_f(f, function, x_axis_grid, y_axis_grid, t+tau*2/3);
        mul_vector_scalar(f, tau);
        vector_plus_vector(c, f);
        mul_A_B_in_C(inverse_A, M, M, c, M, 1, b, M, 1);
        full_f(f, function, x_axis_grid, y_axis_grid, t+tau);
        mul_vector_scalar(f, tau);
        vector_plus_vector(b, f);

        c = b;
        boundary_conditions(b, a, tau, h_x, h_y, bound_1, bound_2, bound_3, bound_4, x_axis_grid, y_axis_grid, t);
        //vector x in Ax = b
        dgetrs_(&N, &M, &NRHS, A.data(), &M, IPIV.data(), b.data(), &M, &INFO);
        out << t+tau << std::endl;
        print_vector(out, b);

        x = b;
    }

    full_bound(left_x, y_axis_grid, bound_1, t);
    full_bound(up_x, x_axis_grid, bound_2, t);
    full_bound(right_x, y_axis_grid, bound_3, t);
    full_bound(down_x, x_axis_grid, bound_4, t);

    calculate_coefficients(A_without_coef, x, c, left_x, up_x, right_x, down_x, h_x, h_y, tau);
//~ print_vector(std::cout, x);
//~ print_vector(std::cout, a);
    mul_vector_scalar(x, -1);

    vector_plus_vector(x, a);
    std::cout << "Square of norm \"calculated coefficients - coefficients\": " << calc_euc_norm_2(x);
    std::cout << std::endl;

    return 0;
}
