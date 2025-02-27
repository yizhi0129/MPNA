set terminal png

set output "gersh_n1.png"

unset key

plot "gersh_n1.txt" using ($2):(0):($3) with circles lc rgb "blue" fill empty