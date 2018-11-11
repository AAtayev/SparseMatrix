#ifndef CLASS_SPARSEMATRIX
#define CLASS_SPARSEMATRIX

#include<iostream>
#include<vector>
#include <string>
#include <algorithm>

class SparseMatrix
{
public:
  SparseMatrix(); // Default Constructor
  SparseMatrix(unsigned int rowSize, unsigned int columnSize); // Default Constructor
  SparseMatrix( const SparseMatrix& matrix); // Copy Constructor
  ~SparseMatrix(); // Destructor

  // Gauss-Seidel
  void GaussSeidel(std::vector<double>& x_0, const double TOL, const int maxIter, const std::vector<double>& b, std::string myName);
  // Boolean operators
  bool operator==(const SparseMatrix& matrix);

  // Operators
  std::vector<double> operator*(const std::vector<double>& vector) const; // Multiplication of SparseMatrix by STL vector

  unsigned int getRowSize() const; // Returns rowSize_ of a SparseMatrix
  unsigned int getColumnSize() const; // Return columnSize_ of a SparseMatrix
  double getEntry(unsigned int rowNum, unsigned int colNum) const; // Returns value of matrix at position (rowNum, colNum)
  void addEntry(unsigned int rowNum, unsigned int colNum, double input); // Inputs values into a SparseMatrix at position (rowNum, colNum)

  void printMatrix(); // Prints SparseMatrix into terminal

private:
  int unsigned rowSize_, columnSize_; // rowSize_ =# of rows, and columnSize_ = # of columns
  std::vector<std::vector<double>* >* rowList_; // The vector containing vectors of rows of a SparseMatrix
  std::vector<std::vector<unsigned int>* >* colIndex_; // The vector containing vectors of column indices for each row of a SparseMatrix
};

std::ostream& operator<<(std::ostream& stream, const SparseMatrix& matrix); // Allows user to type std::cout << A for a SparseMatrix A to print to terminal
double inftyNorm(std::vector<double> vector); // Return the L^{\infty} norm of an STL vector
std::vector<double> v_minus_w(std::vector<double> v, std::vector<double> w); // Performs vector subtraction of two input STL vectors v and w

void GaussSeidelTest(unsigned int N, double TOL, int maxIter, double delta, double lambda, std::string myName); // Constructs matrices and vectors required for Assigment 1 and performs the
                                                                                                                // Gauss-Seidel algorithm on constructed matrix against constructed b vector

#endif
