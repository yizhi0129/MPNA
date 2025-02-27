#!/bin/bash

make clean
make

# matrix vertor multiplication
./mat-vec > n1.txt
mpirun -np 2 ./mat-vec-p >> n1.txt
mpirun -np 4 ./mat-vec-p >> n1.txt
mpirun -np 8 ./mat-vec-p >> n1.txt
mpirun -np 16 ./mat-vec-p >> n1.txt
mpirun -np 32 ./mat-vec-p >> n1.txt
mpirun -np 64 ./mat-vec-p >> n1.txt
mpirun -np 128 ./mat-vec-p >> n1.txt

# error
./accuracy 2 > err_n1.txt
./accuracy 4 >> err_n1.txt
./accuracy 8 >> err_n1.txt
./accuracy 16 >> err_n1.txt
./accuracy 32 >> err_n1.txt
./accuracy 64 >> err_n1.txt
./accuracy 128 >> err_n1.txt

# plot error
gnuplot err_plot1.plt

# plot performance
gnuplot mat_perf_plot1.plt

# eigenvalue problem
./eigenvalue > eig_n1_1.txt
mpirun -np 2 ./eigenvalue-p > eig_n1_2.txt
mpirun -np 4 ./eigenvalue-p > eig_n1_4.txt
mpirun -np 8 ./eigenvalue-p > eig_n1_8.txt
mpirun -np 16 ./eigenvalue-p > eig_n1_16.txt
mpirun -np 32 ./eigenvalue-p > eig_n1_32.txt
mpirun -np 64 ./eigenvalue-p > eig_n1_64.txt
mpirun -np 128 ./eigenvalue-p > eig_n1_128.txt

# Gershgorin circle
./gershgorin > gersh_n1.txt

# plot circles
gnuplot gersh_plot1.plt

# plot performance
gnuplot eig_perf_plot1.plt

# plot iterates
gnuplot eig_iter_plot1.plt