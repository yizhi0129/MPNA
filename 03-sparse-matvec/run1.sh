#!/bin/bash

make clean
make

# matrix vertor multiplication
./mat-vec > mat_vec_time_n1.txt
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

# eigenvalue problem
./eigenvalue > eig_time_n1.txt
mpirun -np 2 ./eigenvalue-p >> eig_time_n1.txt
mpirun -np 3 ./eigenvalue-p >> eig_time_n1.txt
mpirun -np 4 ./eigenvalue-p >> eig_time_n1.txt
mpirun -np 5 ./eigenvalue-p >> eig_time_n1.txt
mpirun -np 6 ./eigenvalue-p >> eig_time_n1.txt
mpirun -np 7 ./eigenvalue-p >> eig_time_n1.txt
mpirun -np 8 ./eigenvalue-p >> eig_time_n1.txt

# Gershgorin circle
./gershgorin > gersh_n1.txt

# plot circles
gnuplot gersh_plot1.plt
