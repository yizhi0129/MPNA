# Sparse Distributed Matrix-Vector Multiplication

## Description

The goal of this exercise is to implement a function that multiplies a sparse distributed matrix by a vector.

Given a sparse matrix `A` of dimensions `n x n` and a vector `x` of dimension `n`, the product of `A` and `x` is a vector `y` of dimension `n` such that `y[i]` is the dot product of the `i`-th row of `A` and the vector `x`.

## Serial CSR Implementation

In the following, you will implement in C (or C++) such a function. The matrix will be represented in the Compressed Sparse Row (CSR) format.
In order to test your implementation you can read a matrix from a file in Matrix Market format, and a vector from a file in plain text format.
Tests matrices and vectors can be found in the data directory and more can be found on https://sparse.tamu.edu/.

## Distributed CSR Implementation

Once you have a working serial implementation, you can extend it to a distributed memory setting. You will use MPI to distribute the matrix among the processes.

You can distribute the matrix in a block-cyclic fashion, where each process will have a block of rows of the matrix. The block size can be fixed or computed based on the number of processes and the matrix size.

You will comment on the performance of your implementation and compare it to the serial version. In particular, you will measure the speedup and efficiency of your parallel implementation.

## Application: Eigenvalue Problem

You will use the matrix-vector multiplication to solve an eigenvalue problem.

### Power Iteration
You will use the power iteration method to compute the largest magnitude eigenvalue of `A`.

Evaluate the convergence of the method as a function of the matrix size and the number of processes.
