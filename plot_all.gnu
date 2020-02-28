set terminal pngcairo size 1800,1000 font ",16"
set output "out.png"


set key o
set multiplot layout 2,2

set xlabel "step"
set ylabel "arbitrary"



p "output.txt" u 2:3 w l lw 2 title "s_x",\
"" u 2:4 w l lw 3 title "s_y",\
"" u 2:5 w l lw 4 title "s_z",\
"" u 2:6 w l lw 5 title "|s|",\

p "output.txt" u 2:7 w l lw 2 title "t_x",\
"" u 2:8 w l lw 3 title "t_y",\
"" u 2:9 w l lw 4 title "t_z",\
"" u 2:10 w l lw 5 title "|t|",\

p "output.txt" u 2:11 w l lw 2 title "grad_sx",\
"" u 2:12 w l lw 3 title "grad_sy",\
"" u 2:13 w l lw 4 title "grad_sz",\
"" u 2:14 w l lw 5 title "grad_lx",\
"" u 2:15 w l lw 6 title "grad_ly",\
"" u 2:16 w l lw 7 title "grad_lz",\
"" u 2:17 w l lw 8 title "grad_l1",\
"" u 2:18 w l lw 9 title "|grad_s|",\
"" u 2:19 w l lw 10 title "|grad_l|"


p "output.txt" u 2:20 w l lw 2 title "l_x",\
"" u 2:21 w l lw 3 title "l_y",\
"" u 2:22 w l lw 4 title "l_z",\
"" u 2:23 w l lw 5 title "l1",\


unset multiplot
