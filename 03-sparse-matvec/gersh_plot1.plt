set terminal png

set output "gersh_n1.png"

set format y ""

plot "gersh_n1.txt" using ($2):(0):($3) with circles lc rgb "blue" fill empty title "Gershgorin", \
     "serial_eig1.txt" using 1:(0) with points title "proc 1", \
     "parallel_eig1_2.txt" using 1:(0) with points title "proc 2", \
     "parallel_eig1_3.txt" using 1:(0) with points title "proc 3", \
     "parallel_eig1_4.txt" using 1:(0) with points title "proc 4", \
     "parallel_eig1_5.txt" using 1:(0) with points title "proc 5", \
     "parallel_eig1_6.txt" using 1:(0) with points title "proc 6", \
     "parallel_eig1_7.txt" using 1:(0) with points title "proc 7", \
     "parallel_eig1_8.txt" using 1:(0) with points title "proc 8"