#!/bin/bash
# change file name in .c files

make clean
make

# matrix vertor multiplication
./mat-vec > n2.txt
mpirun -np 2 ./mat-vec-p >> n2.txt
mpirun -np 4 ./mat-vec-p >> n2.txt
mpirun -np 8 ./mat-vec-p >> n2.txt
mpirun -np 16 ./mat-vec-p >> n2.txt
mpirun -np 32 ./mat-vec-p >> n2.txt
mpirun -np 64 ./mat-vec-p >> n2.txt
mpirun -np 128 ./mat-vec-p >> n2.txt

# eigenvalue problem
./eigenvalue > eig_n2_1.txt
mpirun -np 2 ./eigenvalue-p > eig_n2_2.txt
mpirun -np 4 ./eigenvalue-p > eig_n2_4.txt
mpirun -np 8 ./eigenvalue-p > eig_n2_8.txt
mpirun -np 16 ./eigenvalue-p > eig_n2_16.txt
mpirun -np 32 ./eigenvalue-p > eig_n2_32.txt
mpirun -np 64 ./eigenvalue-p > eig_n2_64.txt
mpirun -np 128 ./eigenvalue-p > eig_n2_128.txt

# Gershgorin circle
./gershgorin > gersh_n2.txt

# plot circles
gnuplot gersh_plot2.plt