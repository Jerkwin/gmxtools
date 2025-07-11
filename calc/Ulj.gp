print 5+lambertw(-5/exp(5))
print 5+lambertw_1(-5/exp(5))

p lambertw(x), lambertw_1(x), lambertw_11(x)

exit
set xl"r(σ)"
set yl"U_{LJ}(ε{/symbolpi e})"
set label "({/symbolpi s}, 0)" at 1.05,0.1
set label "(R_{min}, -{/symbolpi e})" at 1.2,-1
$dat1<<EOD
1 0
EOD
$dat2<<EOD
1.122462048309373 -1
EOD

p [0:3][-1.2:2] 4*(x**-12-x**(-6)) w l lc rgb "#D62728" lw 2 t"U = 4ε {/symbolpi e} [ (σ{/symbolpi s}/r)^{12} - ({/symbolpi s}/r)^6 ]\n   = {/symbolpi e}[(R_{min}/r)^{12}-2(R_{min}/r)^6]", \
0 ls 1 dt 2 t"", $dat1 w p pt 7  lc rgb"#FF7400" t"{/symbolpi s} = R_{min}/2^{1/6} = 0.890899R_{min}", $dat2 w p pt 5  lc rgb"#00A13B" t"R_{min} = {/symbolpi s}*2^{1/6} = 1.122462{/symbolpi s}"
