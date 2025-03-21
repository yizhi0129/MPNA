set terminal png
set logscale y
#set logscale x

# Set common settings for the plot
set xlabel "iterations"
set ylabel "relative residual"

# Loop over different grid sizes (3, 6, 10, 15, 20, 50, 100)
do for [n in "3 6 10 15 20 50 100"] {
    set output sprintf('%s.png', n)  # Set the output file name based on the grid size
    plot sprintf("jacobi_%s.txt", n) using 1:2 with lines title "Jacobi", \
         sprintf("gauss_seidel_%s.txt", n) using 1:2 with lines title "Gauss-Seidel", \
         sprintf("conjugate_gradient_%s.txt", n) using 1:2 with lines title "Conjugate Gradient", \
         sprintf("gmres_%s.txt", n) using 1:2 with lines title "GMRES"
}
