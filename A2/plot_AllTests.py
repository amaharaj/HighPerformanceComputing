set key right top box

set title "Wall Times for Various System Sizes"
set xlabel "d"
set ylabel "Wall Time" 

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5
set style line 2 lc rgb '#2F4F4F' lt 1 lw 2 pt 7 ps 1.5
set style line 3 lc rgb '#8B0000' lt 1 lw 2 pt 7 ps 1.5
set style line 4 lc rgb '#006400' lt 1 lw 2 pt 7 ps 1.5
set style line 5 lc rgb '#FF1493' lt 1 lw 2 pt 7 ps 1.5
set style line 6 lc rgb '#FF4500' lt 1 lw 2 pt 7 ps 1.5
set style line 7 lc rgb '#CD5C5C' lt 1 lw 2 pt 7 ps 1.5
set style line 8 lc rgb '#696969' lt 1 lw 2 pt 7 ps 1.5

set xtics 1

plot 'n500.txt' u 2:3 with linespoints ls 1 title " n=500", 'n1000.txt' u 2:3 with linespoints ls 6 title ' n=1000', 'n5000.txt' u 2:3 with linespoints ls 4 title ' n=5000'
