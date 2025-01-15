#include <iostream>

#include "CSRMatrix.h"
#include "Jacobi.h"
#include "GaussSeidel.h"

template<typename IdType, typename ScalarType>
void benchmark(int mesh_size)
{
    auto A = createLaplacian2D(IdType(9), ScalarType(1.0));

    std::vector<ScalarType> b(A.n, 1.0);
    std::vector<ScalarType> x(A.n, 0.0);

    IdType max_iter = 1000;
    ScalarType tolerance = 1e-6;

    {
        // Solve the system using Jacobi
        auto iter = jacobi(A.n, A, b, x, tolerance, max_iter);
        std::cout << "Jacobi, Number of iterations: " << iter << std::endl;
        // std::cout << "Solution:" << std::endl;
        // for (IdType i = 0; i < A.n; ++i)
        // {
        //     std::cout << x[i] << std::endl;
        // }
    }

    {
        // Solve the system using Gauss-Seidel
        std::fill(x.begin(), x.end(), ScalarType(0)); // reset x
        auto iter = gaussSeidel(A.n, A, b, x, tolerance, max_iter);
        std::cout << "Gauss-Seidel, Number of iterations: " << iter << std::endl;
        // std::cout << "Solution:" << std::endl;
        // for (IdType i = 0; i < A.n; ++i)
        // {
        //     std::cout << x[i] << std::endl;
        // }
    }

}

int main()
{
    benchmark<int, double>(9);
    benchmark<int, float>(9);
    
    return 0;
}