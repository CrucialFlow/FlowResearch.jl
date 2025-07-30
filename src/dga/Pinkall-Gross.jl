# github.com/chakravala

# 1.4 (i) quarter circle

-1/sqrt(2):0.01:1/sqrt(2) → (t->Chain.(t,sqrt(1-t^2)))

# 1.4 (ii) quarter circle

π/4:0.01:3π/4 → (t->Chain.(cos(t),sin(t)))

# 1.4 (iii) full circle

circle(t) = Chain(cos(t),sin(t))
0:0.01:2π → (t->Chain.(cos(t),sin(t)))

# 1.4 (iv) helix

helix(t) = Chain(cos(t),sin(t),t)
0:0.01:2π → (t->Chain.(cos(t),sin(t),t))

# 1.4 (v) cartesian leaf

cartesianleaf(t) = Chain(x^3-4x,x^2-4)
-π:0.01:π → (t->Chain.(t^3-4t,t^2-4))

# cartesian leaf derivative

tangent(-π:0.01:π → (t->Chain.(t^3-4t,t^2-4)))
-π:0.01:π → (t->Chain.(3t^2-4,2t))

# 1.4 (vi) Neil's parabola

neilparabola(t) = Chain(x^3,x^2)

# 5.3 da Rios

start0(x) = Chain(cos(x),sin(x),sin(3x))
start1(x) = Chain(cos(x)+0.1sin(5x),sin(x)+0.1cos(5x),0.3sin(3x))

darios(t,dt=tangent(fiber(t))) = ⋆(dt∧tangent(dt))

sol1 = odesolve(darios,x1,1.0,11)


