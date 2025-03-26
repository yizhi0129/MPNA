#!/bin/bash

make clean
make

mpirun -np 4 ./bin/main

mpirun -np 4 ./bin/main --newton