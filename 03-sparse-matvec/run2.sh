#!/bin/bash
# change file name in .c files

make clean
make

# matrix vertor multiplication
./mat-vec > mat_vec_time_n2.txt
mpirun -np 2 ./mat-vec-p > parallel_result2_2.txt
mpirun -np 3 ./mat-vec-p > parallel_result2_3.txt
mpirun -np 4 ./mat-vec-p > parallel_result2_4.txt
mpirun -np 5 ./mat-vec-p > parallel_result2_5.txt
mpirun -np 6 ./mat-vec-p > parallel_result2_6.txt
mpirun -np 7 ./mat-vec-p > parallel_result2_7.txt
mpirun -np 8 ./mat-vec-p > parallel_result2_8.txt

# eigenvalue problem
./eigenvalue > eig_time_n2.txt
mpirun -np 2 ./eigenvalue-p >> eig_time_n2.txt
mpirun -np 3 ./eigenvalue-p >> eig_time_n2.txt
mpirun -np 4 ./eigenvalue-p >> eig_time_n2.txt
mpirun -np 5 ./eigenvalue-p >> eig_time_n2.txt
mpirun -np 6 ./eigenvalue-p >> eig_time_n2.txt
mpirun -np 7 ./eigenvalue-p >> eig_time_n2.txt
mpirun -np 8 ./eigenvalue-p >> eig_time_n2.txt

# Gershgorin circle
./gershgorin > gersh_n2.txt

# plot circles
gnuplot gersh_plot2.plt