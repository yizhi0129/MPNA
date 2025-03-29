#!/bin/bash

make clean
make

# sigma = 1, beta = 300, kappa = kappa0 * u * u
# gamma = 10

mpirun -np 2 ./bin/main --N 50 > case2_n2_g10_N50.txt
sleep 1
mpirun -np 2 ./bin/main --N 100 > case2_n2_g10_N100.txt
sleep 1
mpirun -np 2 ./bin/main --N 500 > case2_n2_g10_N500.txt
sleep 1
mpirun -np 2 ./bin/main --N 1000 > case2_n2_g10_N1000.txt
sleep 1
mpirun -np 2 ./bin/main --N 5000 > case2_n2_g10_N5000.txt
sleep 1
mpirun -np 2 ./bin/main --N 10000 > case2_n2_g10_N10000.txt
sleep 1

mpirun -np 4 ./bin/main --N 50 > case2_n4_g10_N50.txt
sleep 1
mpirun -np 4 ./bin/main --N 100 > case2_n4_g10_N100.txt
sleep 1
mpirun -np 4 ./bin/main --N 500 > case2_n4_g10_N500.txt
sleep 1
mpirun -np 4 ./bin/main --N 1000 > case2_n4_g10_N1000.txt
sleep 1
mpirun -np 4 ./bin/main --N 5000 > case2_n4_g10_N5000.txt
sleep 1
mpirun -np 4 ./bin/main --N 10000 > case2_n4_g10_N10000.txt
sleep 1

mpirun -np 8 ./bin/main --N 50 > case2_n8_g10_N50.txt
sleep 1
mpirun -np 8 ./bin/main --N 100 > case2_n8_g10_N100.txt
sleep 1
mpirun -np 8 ./bin/main --N 500 > case2_n8_g10_N500.txt
sleep 1
mpirun -np 8 ./bin/main --N 1000 > case2_n8_g10_N1000.txt
sleep 1
mpirun -np 8 ./bin/main --N 5000 > case2_n8_g10_N5000.txt
sleep 1
mpirun -np 8 ./bin/main --N 10000 > case2_n8_g10_N10000.txt    
sleep 1

