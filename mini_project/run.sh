#!/bin/bash

rm -f *.txt
rm -f *.png
rm -f gs

gcc grid_solver.c -o gs -lm

for i in {3,6,10,15,20,50,100}
do
    ./gs $i
done

gnuplot plot.gp