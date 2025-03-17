set terminal png

set output "gersh_n1.png"

unset key

plot "gersh_n1.txt" using ($2):(0):($3) with circles lc rgb "blue" fill empty, \
     "eig_n1_1.txt" using 1:(0) with points title "proc 1", \
     "eig_n1_2.txt" using 1:(0) with points title "proc 2", \
     "eig_n1_3.txt" using 1:(0) with points title "proc 3", \
     "eig_n1_4.txt" using 1:(0) with points title "proc 4", \
     "eig_n1_5.txt" using 1:(0) with points title "proc 5", \
     "eig_n1_6.txt" using 1:(0) with points title "proc 6", \
     "eig_n1_7.txt" using 1:(0) with points title "proc 7", \
     "eig_n1_8.txt" using 1:(0) with points title "proc 8"