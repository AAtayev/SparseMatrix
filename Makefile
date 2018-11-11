# input files
TEX=report.tex
PLOTSCRIPTS=plotscript1.gpl plotscript2.gpl

# generated files
REPORT=report.pdf
RESULTS1=GaussSeidel_delta_1.000000.dat GaussSeidel_delta_10.000000.dat GaussSeidel_delta_50.000000.dat GaussSeidel_delta_100.000000.dat GaussSeidel_delta_500.000000.dat GaussSeidel_delta_1000.000000.dat GaussSeidel_delta_5000.000000.dat GaussSeidel_delta_10000.000000.dat GaussSeidel_delta_50000.000000.dat GaussSeidel_delta_100000.000000.dat
RESULTS2=GaussSeidel_lambda_1.000000.dat GaussSeidel_lambda_10.000000.dat GaussSeidel_lambda_50.000000.dat GaussSeidel_lambda_100.000000.dat GaussSeidel_lambda_500.000000.dat GaussSeidel_lambda_1000.000000.dat GaussSeidel_lambda_5000.000000.dat GaussSeidel_lambda_10000.000000.dat GaussSeidel_lambda_50000.000000.dat GaussSeidel_lambda_100000.000000.dat
RESULTS3=$(wildcard *.dat)
PROGRAM=SparseMatrix
OBJS=main.o SparseMatrix.o
PLOTS=GaussSeidel_delta.pdf GaussSeidel_lambda.pdf

# additional variables
CPPFLAGS=-std=c++11

all: $(REPORT)

$(REPORT): $(PLOTS) $(TEX)
	pdflatex -interaction=batchmode report.tex
	pdflatex -interaction=batchmode report.tex  # do latex twice

$(PLOTS): $(RESULTS1) $(RESULTS2) $(PLOTSCRIPTS)
	gnuplot plotscript1.gpl
	gnuplot plotscript2.gpl

$(RESULTS1) $(RESULTS2) $(RESULTS3): $(PROGRAM)
	./$(PROGRAM)

$(PROGRAM): $(OBJS)
	g++ $(CPPFLAGS) $(OBJS) -o $(PROGRAM)

$(OBJS): %.o: %.cc
	g++ -Ofast -Wall -Wfatal-errors $(CPPFLAGS) -c $^ -o  $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(RESULTS1) $(RESULTS2) $(RESULTS3) $(PLOTS)
