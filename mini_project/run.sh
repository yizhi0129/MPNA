#!/bin/bash

rm -f *.txt
rm -f *.png
rm -f gs

/opt/homebrew/Cellar/gcc/14.2.0_1/bin/gcc-14 grid_solver.c -o gs -I/opt/homebrew/Cellar/open-mpi/5.0.6/include -I/opt/homebrew/Cellar/libomp/19.1.7/include -L/opt/homebrew/Cellar/open-mpi/5.0.6/lib -L/opt/homebrew/Cellar/libomp/19.1.7/lib -lmpi -lomp -fopenmp

for i in {3,6,10,15,20,50,100}
do
    ./gs $i
done

gnuplot plot.gp