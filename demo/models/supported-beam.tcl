# Reference: https://opensees.berkeley.edu/wiki/index.php/Simply_supported_beam_modeled_with_two_dimensional_solid_elements

model BasicBuilder -ndm 2 -ndf 2

# Young's Modulus is enjected
nDMaterial ElasticIsotropic 1 {{{ :E }}} 0.25 3.0

set Quad quad

#set Quad bbarQuad
#set Quad enhancedQuad

if {$Quad == "enhancedQuad" } {
	set eleArgs "PlaneStress2D 1"
}
if {$Quad == "quad" } {
	set eleArgs "1 PlaneStress2D 1"
}
if {$Quad == "bbarQuad" } {
	set eleArgs "1"
}

set nx 8
set ny  2

block2D $nx $ny 1 1 $Quad $eleArgs {
	1  0   0
	2  40  0
	3  40  10
	4  0   10
}

set bn [expr $nx + 1]
set l1 [expr $nx/2 + 1]
set l2 [expr $l1 + $ny*($nx+1)]

fix   1 1 1
fix $bn 0 1

recorder Node -file displacement.out -time -node $l1 -dof 2 disp

pattern Plain 1 Linear {
    load $l1 0.0 -1.0
    load $l2 0.0 -1.0
}

integrator LoadControl 1.0 1 1.0 10.0
test EnergyIncr 1.0e-12  10  0
algorithm Newton
numberer RCM
constraints Plain
system ProfileSPD
analysis Static
analyze 10

wipeAnalysis
setTime 0.0
remove loadPattern 1

rayleigh 0. 0. 0. [expr 2*0.02/sqrt([eigen 1])];

test EnergyIncr 1.0e-12 10 0
algorithm Newton
numberer RCM
constraints Plain
integrator Newmark 0.5 0.25
system BandGeneral
analysis Transient

analyze 1500 0.5