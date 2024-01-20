# 残基颜色 amino shapely
RESINFO="ALA A #C8C8C8 #8CFF8C ARG R #145AFF #00007C ASN N #00DCDC #FF7C70 ASP D #E60A0A #A00042 ASH D #FF69B4 #FF00FF CYS C #E6E600 #EEEE70 CYM C #B4B400 #E6E62E CYX C #FAFA00 #F1F186 GLN Q #00DCDC #FF4C4C GLH E #FF69B4 #FF00FF GLU E #E60A0A #660000 GLY G #EBEBEB #EEEEEE HIS H #8282D2 #7070FF HID H #4D4DBF #8F8FFF HIE H #3A3AA1 #2222FF HIP H #8282D2 #7070FF ILE I #0F820F #004C00 LEU L #0F820F #455E45 LYS K #145AFF #4747B8 LYN K #477EFF #5F6FC7 MET M #E6E600 #B8A042 PHE F #3232AA #534C52 PRO P #DC9682 #525252 SER S #FA9600 #FF7042 THR T #FA9600 #B84C00 TRP W #B45AB4 #4F4600 TYR Y #3232AA #8C704C VAL V #0F820F #FF8CFF UNK U #BEA06E #FF00FF "
resLett(name) = ( _i_=strstrt(RESINFO, name." "), \
	strlen(name)==3 ? RESINFO[_i_+4:_i_+5] : RESINFO[_i_-4:_i_-1])
resColor(name, n) = ( _i_=strstrt(RESINFO, name." ")+8*(n-1), \
	strlen(name)==3 ? RESINFO[_i_+6:_i_+13] : RESINFO[_i_+2:_i_+9])

#HBond #resDonor #resAcceptor Occupancy%
$HBdat <<EOD
   1    2-ILE  15-LYS     98.501
   2    5-ASP   5-ASP      0.200
   3    5-ASP   5-ASP      0.300
   4    6-TYR  33-CYS     19.381
   5    7-GLY   4-GLU     75.624
   6    8-LYS   9-CYS      0.100
   7    8-LYS  29-THR      0.100
   8    8-LYS  30-ASN      2.098
   9    8-LYS  30-ASN      3.297
  10    8-LYS  30-ASN      1.499
  11    8-LYS  30-ASN      0.100
  12    8-LYS  30-ASN      0.400
  13    8-LYS  30-ASN      0.300
  14    9-CYS  31-CYS     98.302
  15   10-THR  14-THR     88.811
  16   10-THR  14-THR      0.400
  17   10-THR   9-CYS      0.400
  18   10-THR  14-THR      0.500
  19   11-TRP  29-THR     10.390
  20   11-TRP  22-CYS     21.279
  21   11-TRP  24-CYS      0.200
  22   11-TRP  28-GLY      1.898
  23   13-GLY  10-THR     57.942
  24   13-GLY  11-TRP      8.392
  25   14-THR  10-THR     16.983
  26   14-THR  10-THR     49.151
  27   14-THR  10-THR      1.399
  28   14-THR  10-THR      2.697
  29   14-THR  13-GLY      0.200
  30   14-THR  14-THR      0.100
  31   15-LYS  12-GLY      0.200
  32   15-LYS  12-GLY      0.500
  33   15-LYS  13-GLY      1.099
  34   15-LYS  13-GLY      2.498
  35   15-LYS  13-GLY      2.198
  36   17-CYS   2-ILE     95.005
  37   18-ARG   3-ALA      0.100
  38   18-ARG   5-ASP      0.999
  39   18-ARG   5-ASP      1.299
  40   18-ARG   3-ALA     13.287
  41   18-ARG   5-ASP     28.072
  42   18-ARG   5-ASP     27.373
  43   19-GLY  16-CYS      0.599
  44   20-ARG  16-CYS      1.199
  45   20-ARG  17-CYS     48.851
  46   20-ARG  18-ARG      0.100
  47   20-ARG  35-PRO      2.498
  48   20-ARG   5-ASP     16.184
  49   20-ARG   5-ASP     29.271
  50   20-ARG  17-CYS      3.097
  51   20-ARG   5-ASP      0.200
  52   20-ARG  18-ARG      0.100
  53   20-ARG  35-PRO      1.099
  54   23-ARG  21-PRO     14.785
  55   23-ARG  32-GLU     44.256
  56   23-ARG  34-THR      2.498
  57   23-ARG  25-SER      3.497
  58   23-ARG  32-GLU      0.999
  59   23-ARG  25-SER     36.164
  60   23-ARG  25-SER      0.400
  61   23-ARG  32-GLU      0.200
  62   23-ARG  32-GLU     11.788
  63   25-SER  23-ARG     32.967
  64   25-SER  30-ASN     18.581
  65   25-SER  24-CYS      0.500
  66   25-SER  32-GLU     14.386
  67   25-SER  32-GLU      1.998
  68   26-MET  24-CYS      1.399
  69   26-MET  25-SER      0.100
  70   27-ILE  24-CYS      2.298
  71   27-ILE  25-SER     34.565
  72   28-GLY  24-CYS      8.292
  73   28-GLY  25-SER     50.250
  74   28-GLY  26-MET      2.098
  75   29-THR  27-ILE     21.778
  76   29-THR  30-ASN      3.197
  77   30-ASN  27-ILE     11.189
  78   30-ASN  28-GLY      2.897
  79   30-ASN  29-THR      0.100
  80   30-ASN   9-CYS     13.187
  81   30-ASN  29-THR      2.897
  82   30-ASN  29-THR      2.498
  83   31-CYS   9-CYS     20.480
  84   32-GLU  23-ARG     99.401
  85   33-CYS   7-GLY     99.800
  86   34-THR  20-ARG      0.100
  87   34-THR  21-PRO     98.701
  88   34-THR  21-PRO      2.498
  89   34-THR  32-GLU     83.916
  90   34-THR  33-CYS      0.200
EOD

resD(i) = int(substr(word($HBdat[i],2), 1, strstrt(word($HBdat[i],2), "-")-1))
resA(i) = int(substr(word($HBdat[i],3), 1, strstrt(word($HBdat[i],3), "-")-1))
occ(i)  = word($HBdat[i],4)

resMin=1E9; resMax=0
do for[i=1:|$HBdat|] {
	resMin=min(resMin, resD(i)); resMax=max(resMax, resD(i))
	resMin=min(resMin, resA(i)); resMax=max(resMax, resA(i))
}
array resName[resMax]
do for[i=1:|$HBdat|] {
	str=word($HBdat[i],2); k=strstrt(str,"-"); idx=int(str[1:k-1]); resName[idx]=str[k+1:k+3]
	str=word($HBdat[i],3); k=strstrt(str,"-"); idx=int(str[1:k-1]); resName[idx]=str[k+1:k+3]
}

Radi=10; dR=1.5; Rext=Radi+dR; mode=3; color=1
eval set_pal('cm_plasma')
eval set_pal('cm_kindlmann 1 0')

set size square
set angle degrees
#set key out noautotitle below
#set grid xtics ls -1 lc rgb "gray"
#set grid ytics ls -1 lc rgb "gray"
#set xrange [-Rext:Rext]; set xtics 1; set format x "";
#set yrange [-Rext:Rext]; set ytics 1; set format y "";
unset border; unset xtic; unset ytics;

ang(i) = 360.*(i-resMin+dt)/(resMax-resMin+2)+90
td(i)  = ang(resD(i)+wd)
ta(i)  = ang(resA(i)+wa)

do for [i=resMin:resMax] {
	if(mode==1) { dt=1; wd=wa=0; 
		set obj i circ at Radi*cos(ang(i)), Radi*sin(ang(i)) size .8 \
			fs solid fc rgb resColor(resName[i], color)
	}
	if(mode==2) { dt=.5; wd=0; wa=0.5; 
		set obj i circ at 0,0 size Radi arc [ang(i):ang(i+1)] lw 10 \
		nowedge fs solid border lc rgb resColor(resName[i], color)
	}
	if(mode==3) { dt=1; w=0.5; wd=-w/2; wa=w/2
		set obj i polygon from \
		(Radi-dR/2)*cos(ang(i-w/2)), (Radi-dR/2)*sin(ang(i-w/2)) to \
		(Radi+dR/2)*cos(ang(i-w/2)), (Radi+dR/2)*sin(ang(i-w/2)) to \
		(Radi+dR/2)*cos(ang(i+w/2)), (Radi+dR/2)*sin(ang(i+w/2)) to \
		(Radi-dR/2)*cos(ang(i+w/2)), (Radi-dR/2)*sin(ang(i+w/2)) to \
		(Radi-dR/2)*cos(ang(i-w/2)), (Radi-dR/2)*sin(ang(i-w/2))    \
		fs solid fc rgb resColor(resName[i], color)
	}
}

# 立方贝塞尔曲线
Bp(t, p0, p1, p2, p3) = p0*(1-t)**3 + 3*p1*t*(1-t)**2 + 3*p2*t*t*(1-t) + p3*t*t*t

Bx(t, x0, x3, r0, a0, r3, a3) = Bp(t, x0, x0+r0*cos(a0), x3+r3*cos(a3), x3)
By(t, y0, y3, r0, a0, r3, a3) = Bp(t, y0, y0+r0*sin(a0), y3+r3*sin(a3), y3)

Bezier(t, xy, p0, p3, r0, a0, r3, a3) = ( \
	Bp(t, p0, \
	(xy eq "x" ? p0+r0*cos(a0) : p0+r0*sin(a0)), \
	(xy eq "x" ? p3+r3*cos(a3) : p3+r3*sin(a3)), p3) )
#plot [0:1] \
	Bx(t, 0, 5, r0, 180, r3, 90),      \
	By(t, 0, 5, r0, 180, r3, 90) w l,  \
	Bezier(t, "x", 0, 4, r0, 180, r3, 90 ), \
	Bezier(t, "y", 0, 4, r0, 180, r3, 90 ) w l

set parametric
set angle degrees
plot \
	for [i=1:|$HBdat|] [0:1:.1] '+'  \
		u (Bx(t, Radi*cos(td(i)), Radi*cos(ta(i)), Radi/2, td(i)+180, Radi/2, ta(i)+180)) \
		: (By(t, Radi*sin(td(i)), Radi*sin(ta(i)), Radi/2, td(i)+180, Radi/2, ta(i)+180)) \
		w l lw 1 lc rgb rgba(occ(i),0,70,.5 ) , \
	for [i=resMin:resMax] [0:0] '+' \
		u ((Radi+dR)*cos(ang(i))):((Radi+dR)*sin(ang(i))):(resName[i])\
		: (ang(i)<270?ang(i)-180:ang(i)) \
		w labels rot var t""

exit

# 弯曲贝塞尔箭头

#    P0x  P0y   A0   R0   P3x  P3y  A3    R3   mode  color
$Arrows <<EOD
 1    1    1     0   0.5   3    3     0   0.5  0   0xff0000
 2    3    1     0   0.5   5    3     0   0.5  1   0x00c000
 3    5    1     0   0.5   7    3     0   0.5  2   0x0000ff
 4    7    1     0   0.5   9    3     0   0.5  3   0xff00ff
 5    1    4     0   0.5   3    6    90   0.5  0   0xff0000
 6    3    4     0   0.5   5    6    90   0.5  1   0x00c000
 7    5    4     0   0.5   7    6    90   0.5  2   0x0000ff
 8    7    4     0   0.5   9    6    90   0.5  3   0xff00ff
 9    1    7    90   0.5   3    9     0   0.5  0   0xff0000
10    3    7    90   0.5   5    9     0   0.5  1   0x00c000
11    5    7    90   0.5   7    9     0   0.5  2   0x0000ff
12    7    7    90   0.5   9    9     0   0.5  3   0xff00ff
13   11    1    45   0.5  13    3   -45   0.5  0   0xff0000
14   13    1    45   0.5  15    3   -45   0.5  1   0x00c000
15   15    1    45   0.5  17    3   -45   0.5  2   0x0000ff
16   17    1    45   0.5  19    3   -45   0.5  3   0xff00ff
17   11    4   -45   0.5  13    6   -45   0.5  0   0xff0000
18   13    4   -45   0.5  15    6   -45   0.5  1   0x00c000
19   15    4   -45   0.5  17    6   -45   0.5  2   0x0000ff
20   17    4   -45   0.5  19    6   -45   0.5  3   0xff00ff
21   11    7     0   0.5  15    9    90   0.5  1   0x00c000
22   15    7     0   0.5  19    9     0   0.5  1   0x00c000
EOD

# 参数 mode 决定为哪些端点指定相对或绝对角度, 二进制位选项
# 0: 两端点皆为相对角度  0 0
# 1: 始点相对, 终点绝对  0 1
# 2: 始点绝对, 终点相对  1 0
# 3: 两端点皆为绝对角度  1 1

AngleMode(i) = int(word($Arrows[i],10))

# 始点 P0 与终点 P3 之间的角度(-90° <= angle < 270°, 若两点重合, 则为NaN)
AngleP0P3(n)  = (dy=p3y-p0y, dx=p3x-p0x)==0 ? (dy==0 ? NaN : sgn(dy)*90) : \
                (dx<0 ? 180 : 0) + atan(dy/dx)

# 根据 mode 计算相应的角度
Angle(i, p) = word($Arrows[i], p) \
			+ ((p==4 && AngleMode(i)&2) || (p==8 && AngleMode(i)&1) ? 0 : AngleP0P3(0))

# 半径为 P0-P3 长度的倍数
Radiius(i, p) = word($Arrows[i], p) * sqrt((p3x-p0x)**2 + (p3y-p0y)**2)

# 初始化箭头 #i
ArrowInit(i) = (p0x=word($Arrows[i],2), p0y=word($Arrows[i],3), \
				p3x=word($Arrows[i],6), p3y=word($Arrows[i],7), \
				p1x=p0x+Radiius(i,5)*cos(Angle(i,4)), \
				p1y=p0y+Radiius(i,5)*sin(Angle(i,4)), \
				p2x=p3x-Radiius(i,9)*cos(Angle(i,8)), \
				p2y=p3y-Radiius(i,9)*sin(Angle(i,8))  )

Color(i)  = word($Arrows[i],11) # 颜色

# 立方贝塞尔曲线, 参数 t[0:1], 始 P0, 终 P3, 控制 P1, P2
Bx(t) = p0x*(1-t)**3 + 3*p1x*t*(1-t)**2 + 3*p2x*t*t*(1-t) + p3x*t*t*t
By(t) = p0y*(1-t)**3 + 3*p1y*t*(1-t)**2 + 3*p2y*t*t*(1-t) + p3y*t*t*t

N=|$Arrows|
# 设置线型和箭头类型
do for [i=1:N] {
	set style line i lw 1 lc rgb Color(i)
	set style arrow i head size .2,15,45 fixed filled ls i
}

set xr [0:20]
set yr [0:10]

set parametric # 参数方程
plot \
	for [i=1:N] [0:1:.05] '+' u (ArrowInit(i), Bx($1)):(By($1)) w l ls i t"", \
	for [i=1:N] [0:0]     '+' u (ArrowInit(i),Bx(.9)):(By(.9)):(Bx(1)-Bx(.9)):(By(1)-By(.9)) w vec as i t"", \
	$Arrows u 2:3:1 w labels offset -.5,-1 t"", \
	keyentry w l ls 1 t "relative", \
	keyentry w l ls 2 t "relative -> absolute           ", \
	keyentry w l ls 3 t "absolute -> relative           ", \
	keyentry w l ls 4 t "absolute"

exit
