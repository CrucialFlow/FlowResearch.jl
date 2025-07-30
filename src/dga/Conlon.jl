# github.com/chakravala
# DGA: Differentiable Manifolds, Conlon

# 1.5.6

rev = revolve(unitcircle(TensorField(LinRange(-π/2,π/2,60))))
f(x) = Chain(x[2]*x[3],x[1]*x[3],x[1]*x[2])
wireframe(f.(rev))

# 2.6 Construction of Smooth Functions

t = TensorField(0:0.1:5)
h(x) = x[1]>0 ? exp(-1/x[1]^2) : zero(x[1])
k(x,a,b) = h.(b-x)*h.(x-a)
lines(h.(t))
lines(k(t,1,4))

g(x) = k(x[1])*k(x[2]) # *k(x[i])*...

integral(k(t,1,4))/integrate(k(t,1,4))

# 4.5.16

f(x) = (1-x[1]^2-x[2]^2)*exp(x[3])
set = fol.(OpenParameter(-1.5:0.1:1.5,-1.5:0.1:1.5,0:0.1:2))
contour(set)

