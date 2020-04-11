set sel1 [ atomselect [molinfo top] "serial 1 to 6" ]
set sel2 [ atomselect [molinfo top] "serial 13 to 18" ]

source pistack_linalg.tcl

set Nfrm [molinfo top get numframes]

set log [ open "pistack.xvg" w ]
puts -nonewline $log "#Frame  Xcnt(1)  Ycnt(1)  Zcnt(1)  Xcnt(2)  Ycnt(2)  Zcnt(2)"
puts $log "   a(1)     b(1)     c(1)     a(2)     b(2)     c(2)   A(n1,n2)  D(c1,p2) D(c2,p1)"

for {set i 1} {$i<$Nfrm} {incr i} {
	$sel1 frame $i; $sel1 update
	$sel2 frame $i; $sel2 update

	set cnt1 [measure center $sel1]
	set cnt2 [measure center $sel2]

	set A {}; set B {}
	foreach xyz [$sel1 get {x y z}] {
		lappend A $xyz
		lappend B 1
	}
	set d1 -1
	set n1 [ ::math::linearalgebra::leastSquaresSVD $A $B ]
	if {[lindex $n1 2] < 0} {
		set d1 1
		set n1 [vecscale -1 $n1]
	}

	set A {}; set B {}
	foreach xyz [$sel2 get {x y z}] {
		lappend A $xyz
		lappend B 1
	}
	set d2 -1
	set n2 [ ::math::linearalgebra::leastSquaresSVD $A $B ]
	if {[lindex $n2 2] < 0} {
		set d2 1
		set n2 [vecscale -1 $n2]
	}

	set ang [expr [veclength $n1]*[veclength $n2] ]
	set ang [expr [vecdot $n1 $n2]/$ang ]
	set ang [expr acos($ang)*180/$M_PI ]

	set r12 [expr abs([vecdot $cnt1 $n2]+$d2)/[veclength $n2] ]
	set r21 [expr abs([vecdot $cnt2 $n1]+$d1)/[veclength $n1] ]

	puts -nonewline $log [format "%6d" $i]
	foreach x $cnt1 { puts -nonewline $log [format "%9.4f" $x] }
	foreach x $cnt2 { puts -nonewline $log [format "%9.4f" $x] }
	foreach x $n1   { puts -nonewline $log [format "%9.4f" $x] }
	foreach x $n2   { puts -nonewline $log [format "%9.4f" $x] }
	puts $log [format "%9.3f%9.4f%9.4f" $ang $r12 $r21]
}
flush $log
close $log

global vmd_frame
trace variable vmd_frame([molinfo top]) w frame_draw

draw delete all
proc vmd_draw_arrow {mol start end} {
	set middle [vecadd $start [vecscale 0.75 [vecsub $end $start]]]
	graphics $mol cylinder $start $middle radius 0.05
	graphics $mol cone $middle $end radius 0.2
}

proc frame_draw {name element op} {
	global M_PI; global sel1; global sel2

	set i [molinfo top get frame]

	$sel1 frame $i; $sel1 update
	$sel2 frame $i; $sel2 update

	set cnt1 [measure center $sel1]
	set cnt2 [measure center $sel2]

	set A {}; set B {}
	foreach xyz [$sel1 get {x y z}] {
		lappend A $xyz
		lappend B 1
	}
	set d1 -1
	set n1 [ ::math::linearalgebra::leastSquaresSVD $A $B ]
	if {[lindex $n1 2] < 0} {
		set d1 1
		set n1 [vecscale -1 $n1]
	}

	set A {}; set B {}
	foreach xyz [$sel2 get {x y z}] {
		lappend A $xyz
		lappend B 1
	}
	set d2 -1
	set n2 [ ::math::linearalgebra::leastSquaresSVD $A $B ]
	if {[lindex $n2 2] < 0} {
		set d2 1
		set n2 [vecscale -1 $n2]
	}

	set ang [expr [veclength $n1]*[veclength $n2] ]
	set ang [expr [vecdot $n1 $n2]/$ang ]

	draw delete all
	draw color yellow
	draw material Transparent

	draw arrow $cnt1 [vecadd $cnt1 [vecscale 5 [vecnorm $n1]] ]
	draw arrow $cnt2 [vecadd $cnt2 [vecscale 5 [vecnorm $n2]] ]
	draw cylinder $cnt1 $cnt2 radius 0.05

	set d [expr -1*([vecdot $cnt1 $n2]+$d2)/[veclength $n2] ]
	set p [vecadd $cnt1 [vecscale $d [vecnorm $n2]] ]
	set txt [format "D(cnt1-plane2)=%6.2f" [expr abs($d)]]
	draw text $cnt1 "cnt1" size 1 thickness 2
	draw text $p    "$txt" size 1 thickness 2
	draw triangle $cnt1 $cnt2 $p

	set d [expr -1*([vecdot $cnt2 $n1]+$d1)/[veclength $n1] ]
	set p [vecadd $cnt2 [vecscale $d [vecnorm $n1]] ]
	set txt [format "D(cnt2-plane1)=%6.2f" [expr abs($d)]]
	draw text $cnt2 "cnt2" size 1 thickness 2
	draw text $p    "$txt" size 1 thickness 2
	draw triangle $cnt1 $cnt2 $p

	set p [vecscale 0.5 [vecadd $cnt1 $cnt2] ]
	set txt [format "%6.1f deg" [expr acos($ang)*180/$M_PI]]
	draw text $p "$txt" size 1 thickness 2
	puts $xxx
}
