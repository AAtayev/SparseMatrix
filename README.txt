#######################################################
              C1: Assignment 1 Submission
              Student Number: u1858921
#######################################################

This folder should contain the following following to produce the same outputs
as in the report, report.pdf, along with the report, simply type "make" in the
terminal without the quotation marks "". To clean the data files along with
other files produced when having run the code, in the terminal type
"make clean", again without the use of "".

To produce results for different values of \delta and \lambda used in the
report, access the main.cc file provided and simply type, within the curly
brackets of int main {} the code
GaussSeidelTest(N = Size of Matrix, TOL = tolerance, maxIter = max iterations,
                delta, lambda, File),
where File is the name of data file that will be produced which stores
in two columns the iteration and the norm of the residual at said iteration.
To print the data file for N > 5000, uncomment the if statement stopping this
from occurring within GaussSeidel implementation in the SparseMatrix.cc file.
Once all the changes have been make, simply type "make" within the terminal
once more.
