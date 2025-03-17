#!/bin/bash

make clean
make

# matrix vertor multiplication
./mat-vec > n1.txt
mpirun -np 2 ./mat-vec-p > parallel_result1_2.txt
mpirun -np 3 ./mat-vec-p > parallel_result1_3.txt
mpirun -np 4 ./mat-vec-p > parallel_result1_4.txt
mpirun -np 5 ./mat-vec-p > parallel_result1_5.txt
mpirun -np 6 ./mat-vec-p > parallel_result1_6.txt
mpirun -np 7 ./mat-vec-p > parallel_result1_7.txt
mpirun -np 8 ./mat-vec-p > parallel_result1_8.txt

# error
./accuracy 2 > err_n1.txt
./accuracy 3 >> err_n1.txt
./accuracy 4 >> err_n1.txt
./accuracy 5 >> err_n1.txt
./accuracy 6 >> err_n1.txt
./accuracy 7 >> err_n1.txt
./accuracy 8 >> err_n1.txt

# plot error
gnuplot err_plot1.plt

# plot performance
gnuplot mat_perf_plot1.plt

# eigenvalue problem
./eigenvalue > eig_n1_1.txt
mpirun -np 2 ./eigenvalue-p > eig_n1_2.txt
mpirun -np 3 ./eigenvalue-p > eig_n1_3.txt
mpirun -np 4 ./eigenvalue-p > eig_n1_4.txt
mpirun -np 5 ./eigenvalue-p > eig_n1_5.txt
mpirun -np 6 ./eigenvalue-p > eig_n1_6.txt
mpirun -np 7 ./eigenvalue-p > eig_n1_7.txt
mpirun -np 8 ./eigenvalue-p > eig_n1_8.txt

# Gershgorin circle
./gershgorin > gersh_n1.txt

# plot circles
gnuplot gersh_plot1.plt

# plot performance
gnuplot eig_perf_plot1.plt

# plot iterates
gnuplot eig_iter_plot1.plt