#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
#include "SparseMatrix.hh"

int main ()
{
  auto start = std::chrono::high_resolution_clock::now();
  // Gauss-Seidel Testing
  int N = 100; // Size of matrix A
  double TOL = 1e-6; // Tolerance
  int maxIter = 10000000; // Maximum number of iterations before exiting the program
  std::vector<int> Delta = {1,10,50,100,500,1000,5000,10000,50000,100000}; // Storing different values of delta
  std::vector<int> Lambda = {1,10,50,100,500,1000,5000,10000,50000,100000}; // Storing different values of lambda

  // For loops executing the GaussSeidel algorithm for the required tests
  for (double delta : Delta)
  {
    GaussSeidelTest(N, TOL, maxIter, delta, 0, "GaussSeidel_delta_"+std::to_string(delta)+".dat");
  }
  for (double lambda : Lambda)
  {
    GaussSeidelTest(N, TOL, maxIter, 1, lambda, "GaussSeidel_lambda_"+std::to_string(lambda)+".dat");
  }

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count() << " s\n";
  return 0;
}
