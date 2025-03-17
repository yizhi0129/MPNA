set terminal png

set output "gersh_n2.png"

unset key

plot "gersh_n2.txt" using ($2):(0):($3) with circles lc rgb "blue" fill empty, \
     "eig_n2_1.txt" using 1:(0) with points title "proc 1", \
     "eig_n2_2.txt" using 1:(0) with points title "proc 2", \
     "eig_n2_3.txt" using 1:(0) with points title "proc 3", \
     "eig_n2_4.txt" using 1:(0) with points title "proc 4", \
     "eig_n2_5.txt" using 1:(0) with points title "proc 5", \
     "eig_n2_6.txt" using 1:(0) with points title "proc 6", \
     "eig_n2_7.txt" using 1:(0) with points title "proc 7", \
     "eig_n2_8.txt" using 1:(0) with points title "proc 8"