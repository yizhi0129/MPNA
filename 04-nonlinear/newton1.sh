#!/bin/bash

make clean
make

# sigma = 0.1, beta = 1, kappa = kappa0 * sqrt(u) 

mpirun -np 2 ./bin/main --N 50 --newton > case1_n2_N50_newton.txt
sleep 1
mpirun -np 2 ./bin/main --N 100 --newton > case1_n2_N100_newton.txt
sleep 1
mpirun -np 2 ./bin/main --N 500 --newton > case1_n2_N500_newton.txt
sleep 1
mpirun -np 2 ./bin/main --N 1000 --newton > case1_n2_N1000_newton.txt
sleep 1
mpirun -np 2 ./bin/main --N 5000 --newton > case1_n2_N5000_newton.txt
sleep 1
mpirun -np 2 ./bin/main --N 10000 --newton > case1_n2_N10000_newton.txt
sleep 1

mpirun -np 4 ./bin/main --N 50 --newton > case1_n4_N50_newton.txt
sleep 1
mpirun -np 4 ./bin/main --N 100 --newton > case1_n4_N100_newton.txt
sleep 1
mpirun -np 4 ./bin/main --N 500 --newton > case1_n4_N500_newton.txt
sleep 1
mpirun -np 4 ./bin/main --N 1000 --newton > case1_n4_N1000_newton.txt
sleep 1
mpirun -np 4 ./bin/main --N 5000 --newton > case1_n4_N5000_newton.txt
sleep 1
mpirun -np 4 ./bin/main --N 10000 --newton > case1_n4_N10000_newton.txt
sleep 1

mpirun -np 8 ./bin/main --N 50 --newton > case1_n8_N50_newton.txt
sleep 1
mpirun -np 8 ./bin/main --N 100 --newton > case1_n8_N100_newton.txt
sleep 1
mpirun -np 8 ./bin/main --N 500 --newton > case1_n8_N500_newton.txt
sleep 1
mpirun -np 8 ./bin/main --N 1000 --newton > case1_n8_N1000_newton.txt
sleep 1
mpirun -np 8 ./bin/main --N 5000 --newton > case1_n8_N5000_newton.txt
sleep 1
mpirun -np 8 ./bin/main --N 10000 --newton > case1_n8_N10000_newton.txt
sleep 1

