#!/bin/bash
# change file name in .c files

make clean
make

# matrix vertor multiplication
./mat-vec > n2.txt
mpirun -np 2 ./mat-vec-p >> n2.txt
mpirun -np 3 ./mat-vec-p >> n2.txt
mpirun -np 4 ./mat-vec-p >> n2.txt
mpirun -np 5 ./mat-vec-p >> n2.txt
mpirun -np 6 ./mat-vec-p >> n2.txt
mpirun -np 7 ./mat-vec-p >> n2.txt
mpirun -np 8 ./mat-vec-p >> n2.txt

# eigenvalue problem
./eigenvalue > eig_n2_1.txt
mpirun -np 2 ./eigenvalue-p > eig_n2_2.txt
mpirun -np 3 ./eigenvalue-p > eig_n2_3.txt
mpirun -np 4 ./eigenvalue-p > eig_n2_4.txt
mpirun -np 5 ./eigenvalue-p > eig_n2_5.txt
mpirun -np 6 ./eigenvalue-p > eig_n2_6.txt
mpirun -np 7 ./eigenvalue-p > eig_n2_7.txt
mpirun -np 8 ./eigenvalue-p > eig_n2_8.txt

# Gershgorin circle
./gershgorin > gersh_n2.txt

# plot circles
gnuplot gersh_plot2.plt