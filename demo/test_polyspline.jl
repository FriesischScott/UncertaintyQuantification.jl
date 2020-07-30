using UncertaintyQuantification, Plots

#Test polyspline for function f(x, y) = x*y

A = [1 2; 1 3; 3 6; 4 3; 5 7; 5 2; 6 8; 6 9; 3 7]
f = [2; 3; 18; 12; 35; 10; 48; 54; 21]
k = 1
test = createpolyspline(A, f, k, :fkt)

test_ps_1 = calcpolyspline(test, [4; 3])

#Test polyspline with plot
#compare wikipedia: https://en.wikipedia.org/wiki/File:Polyharmonic-splines-example1.png

B = [0; 0.4; 1.35; 2]
fb = [0; -0.3; 0.4; -0.6]
kb = 3

test2 = createpolyspline(B, fb, kb, :fkt)

test_ps_2 = calcpolyspline(test2, [2])

x = collect(0:.1:2)
y =  Vector{Float64}(undef, 21)
for i in 1:21
    y[i] = calcpolyspline(test2, [x[i]])
end
plot(x, y)
