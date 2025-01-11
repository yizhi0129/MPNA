#!/bin/bash

rm -f *.txt
rm -f *.png
rm -f gs

gcc grid_solver.c -o gs

for i in {3,6,10,15,20}
do
    ./gs $i
done

gnuplot plot.gp