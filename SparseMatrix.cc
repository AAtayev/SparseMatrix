#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include "SparseMatrix.hh"

// Default constructor
SparseMatrix::SparseMatrix()
{
  // As a default, returns a 1x1 sparse matrix with no assigned entries
  rowSize_ = 1;
  columnSize_ = 1;
  rowList_ = new std::vector<std::vector<double>* > (1);
  colIndex_ = new std::vector<std::vector<unsigned int>* > (1);
}

SparseMatrix::SparseMatrix( unsigned int rowSize, unsigned int columnSize )
{
  rowSize_ = rowSize; // Sets the private variable rowSize_ to rowSize
  columnSize_ = columnSize; // Sets the private variable columnSize_ to columnSize
  if ( rowSize < 1 || columnSize < 1 ) // Prints error message if asked to produce an invalid matrix
  {
    throw std::invalid_argument("Cannot create matrix of size 0 or negative size");
  }
  rowList_ = new std::vector<std::vector<double>* > (rowSize); // Creates a vector of points of a vector of pointers of size rowSize
  colIndex_= new std::vector<std::vector<unsigned int>* > (rowSize); // Creates a vector of points of a vector of pointers of size rowSize
}

// Copy constructor
SparseMatrix::SparseMatrix( const SparseMatrix& matrix )
{
  // Copies and returns a given SparseMatrix
  rowSize_ = matrix.rowSize_;
  columnSize_ = matrix.columnSize_;
  rowList_ = matrix.rowList_;
  colIndex_ = matrix.colIndex_;
}

// Destructor
SparseMatrix::~SparseMatrix()
{
  if (rowList_) // If the rowList_ is non-empty...
  { // ... cycle through the vectors of rowList_ and colIndex_ and delete them
    for (std::vector<double>* v : *rowList_)
    {
      delete v;
    }
    for (std::vector<unsigned int>* v : *colIndex_)
    {
      delete v;
    }
  }
  // Once each of the entries have been deleted, delete rowList_ and colIndex_
  delete(rowList_);
  delete(colIndex_);
}

// getRowSize
unsigned int SparseMatrix::getRowSize() const
{
  return rowSize_; // Returns the number of rows a SparseMatrix variable has been assigned
}

// getColumnSize
unsigned int SparseMatrix::getColumnSize() const
{
  return columnSize_; // Returns the number of columns a SparseMatrix variable has been assigned
}

// addEntry
void SparseMatrix::addEntry(unsigned int rowNum, unsigned int colNum, double input)
{
  if (rowNum > rowSize_) // Performs checks in whether it is possible to add in specified position
  {
    throw std::invalid_argument("Row input out of range");
  }
  else if (colNum > columnSize_) // Performs checks in whether it is possible to add in specified position
  {
    throw std::invalid_argument("Column input out of range");
  }
  else // Otherwise...
  {
    if(input != 0) // ...if the input is 0, do nothing since we do not want to store zero entries. Else...
    {
      std::vector<double>* accessedRowValues = rowList_->at(rowNum); // ... copy the row of values that we want to work with
      std::vector<unsigned int>* accessedColumnIndices = colIndex_->at(rowNum); // ... copy the row of column indices that we want to work with
      if(accessedRowValues == 0) // If accessedRowValues are empty, create it and input the given values into them
      {
        rowList_->at(rowNum) = new std::vector<double>(1,input);
        colIndex_->at(rowNum) = new std::vector<unsigned int>(1,colNum);
      }
      else // Otherwise, push back into each respective vector (notice that this will mean that each row value and its corresponding index will be...
        // stored in the same position within the row. We utilise this in getEntry)
        {
          accessedRowValues->push_back(input);
          accessedColumnIndices->push_back(colNum);
        }
      }
    }
  }

// getEntry
double SparseMatrix::getEntry(unsigned int rowNum, unsigned int colNum) const
{
  if (rowNum > rowSize_) // Error checks whether row/column inputs are possible
  {
    throw std::invalid_argument("Row input out of range");
  }
  else if (colNum > columnSize_) // Error checks whether row/column inputs are possible
  {
    throw std::invalid_argument("Column input out of range");
  }
  if ((colIndex_)->at(rowNum) == 0) // If an entry has not been assigned a value (appears as 0), then return 0
  {
    return 0.0;
  }
  else // otherwise, access the row we want the entry from and then search for the required column
  {
    std::vector<unsigned int>* accessedColumnIndices = colIndex_->at(rowNum); // Copy the row of column indices that we want to work with
    int counter = 0; // Initialise counter to zero
    while (colNum != accessedColumnIndices->at(counter)) // Search for the position of the column index we're looking for in this row of indices
    {
      counter++;
    }
    return (rowList_->at(rowNum))->at(counter); // Return the value in the row of variables at found position
  }
}

// printMatrix: this function prints the given matrix into the terminal in a pleasant format
void SparseMatrix::printMatrix()
{
  for (unsigned int i = 0; i < rowSize_; ++i)
  {
    for (unsigned int j = 0; j < columnSize_; ++j)
    {
      std::cout.width(5);
      std::cout << getEntry(i,j);
    }
    std::cout << std::endl;
  }
}

// Operators
// Multiplication of matrix by a vector on the right
std::vector<double> SparseMatrix::operator*(const std::vector<double>& vector) const
{
  std::vector<double> tempVector(rowSize_);
  for (unsigned int i = 0; i < rowSize_; ++i)
  {
    double sum = 0;
    for (unsigned int j: *colIndex_->at(i)) // Only cycle through the non-zero entries to avoid over calculating
    {
      sum += getEntry(i,j)*vector[j];
    }
    tempVector.at(i) = sum;
  }
  return tempVector;
}

// ostream operator, outputs matrices, without having to type .printMatrix(), into the terminal. Works very similar to printMatrix
std::ostream& operator<<(std::ostream& stream, const SparseMatrix& matrix)
{
  for(unsigned int i = 0; i < matrix.getRowSize(); ++i)
  {
    for (unsigned int j = 0; j < matrix.getColumnSize(); ++j)
    {
      stream.width(5);
      stream << matrix.getEntry(i,j);
    }
    stream << std::endl;
  }
  return stream;
}

// Gauss-Seidel: takes as arguments an initial guess x^{(0)}, a tolerance TOL, maximum number of iterations maxIter, vector to invert A against, b, and a name of output file File
void SparseMatrix::GaussSeidel(std::vector<double>& x_0, const double TOL, const int maxIter, const std::vector<double>& b, std::string myName)
{
  std::ofstream myFile;
  myFile.open(myName, std::ios::out); // Creates and opens a file in the name of myName
  if (!myFile.good())
  {
    throw std::invalid_argument("Failed to open file");
  }
  myFile << "This is a Gauss Seidel data file." << std::endl;
  myFile.width(20);
  myFile << std::left << "Iteration k" << "norm(residual_k)" << std::endl;
  int k = 0;
  double sum; // Define sum variable to be used in for loop
  std::vector<double> residual_k = v_minus_w(b, (*this)*x_0);
  myFile.width(20);
  myFile << k << inftyNorm(residual_k) << std::endl;
  while(inftyNorm(residual_k) > TOL && k <= maxIter)
  {
    for (unsigned int i = 0; i < rowSize_; ++i)
    {
      sum = 0; // Initialise the sum to zero for each row
      for (unsigned int j : *colIndex_->at(i))
      {
        if (j!=i) // produces \sum_{i=1, i!=j}^{N}a_{ij}x_j^{(k)}
        {
          sum += getEntry(i,j)*x_0[j]; // We note that here we have x_0 because we update the value of x_0 by reference
        }
      }
      x_0[i] = (b[i] - sum)/(getEntry(i,i));
    }
    k++;
    residual_k = v_minus_w(b, (*this)*x_0); // Compute the residual at the kth iteration
    myFile.width(20);
    myFile << k << inftyNorm(residual_k) << std::endl; // Print the iteration and its corresponding residual at this iteration
  }
  std::cout << "Number of iterations: " << k << std::endl; // Prints to terminal the number of iterations taken to perform algorithm
  myFile.close(); // Closes file myFile
}

// Non-member class functions
// inftyNorm: computes the L^{\infty} norm of a given STL vector
double inftyNorm(std::vector<double> vector)
{
  double max = fabs(vector[0]); // Initialise the maximum we want to find with the first term of the vector
  for (unsigned int i = 1; i < vector.size(); ++i) // Run through the vector, excluding the first term since this is not needed
  {
    if (max < fabs(vector[i])) // Compare the max variable with the [i]th term in the vector, if term[i] is bigger, then update max with term[i]
    {
      max = fabs(vector[i]);
    }
  }
  return max;
  // Note that we use fabs instead of abs since fabs returns a double but abs returns an int
}

// v_minus_w: performs vector subtraction
std::vector<double> v_minus_w(std::vector<double> v, std::vector<double> w)
{
  std::vector<double> temp(v.size());
  for (unsigned int i = 0; i < v.size(); ++i)
  {
    temp.at(i) = v.at(i) - w.at(i);
  }
  return temp;
}

// GaussSeidelTest: Constructs matrix A and vector b required for Assigment 1, and performs the Gauss-Seidel algorithm on A against b
//                : Takes inputs size of the matrix N, tolerance TOL, maximum number of iterations maxIter, delta and lambda (as in the assignment) and a file name used in GaussSeidel
void GaussSeidelTest(unsigned int N, double TOL, int maxIter, double delta, double lambda, std::string myName)
{
  std::cout << "delta = " << delta << ",\t lambda = " << lambda << std::endl;
  std::vector<double> x_0(N,0);
  double a = 4*(1-delta);
  std::vector<double> w(N);
  std::vector<double> D(N+1);
  std::vector<double> b(N);
  // Initialising w and b
  for (unsigned int i = 0; i < N; ++i)
  {
    w[i] = (i + 1)/(double)(N + 1);
    b[i] = -2*a*(w[i] - 0.5)*w[0]*w[0];
  }
  b[N-1] += 1;
  // Initialising D
  for (unsigned int i = 1; i < N + 1; ++i)
  {
    D[i] = a*(w[i-1] - 0.5)*(w[i-1] - 0.5) + delta;
  }
  D[0] = D[1];
  SparseMatrix A = SparseMatrix(N,N);
  // Initialisng A
  for (unsigned int i = 0; i < N; ++i)
  {
    for (unsigned int j = 0; j < N; ++j)
    {
      if (j == i - 1)
      {
        A.addEntry(i, j, -D[i]);
      }
      else if (j == i)
      {
        A.addEntry(i, j, D[i+1] + D[i] + lambda);
      }
      else if (j == i + 1)
      {
        A.addEntry(i, j, -D[i+1]);
      }
      else
      {continue;}
    }
  }
  // A.printMatrix();
  A.GaussSeidel(x_0, TOL, maxIter, b, myName);
  if (delta == 1 && lambda == 0)
  {
    std::cout << "Error: " << inftyNorm(v_minus_w(w, x_0)) << std::endl;
  }
}
