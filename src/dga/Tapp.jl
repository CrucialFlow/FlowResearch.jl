
# Example 1.2

lines(unitcircle())

# Example 1.3

t = TensorField(range(0,5pi,100))
lines(unithelix(t))

# Example 1.4

t = TensorField(-2:0.01:2)
lines(Chain.(t,t^2))

# Example 1.5

t = TensorField(-1:0.01:3)
lines(Chain(2,4)+t.*Chain(3,-1))

# Example 1.10

t = TensorField(-2:0.01:2)
lines(Chain.(t^3,t^2))

# Exercise 1.7

t = TensorField(pi/2,5pi/2,100)
lines(unitcircle(t))

# Exercise 1.8

t = TensorField(0:0.001:1)
arclength(Chain.(2t,3t^2))

# Exercise 1.12

t = TensorField(range(0,2pi,100))

lemniscate(t) = Chain.(cos(t)/(1+sin(t)^2),sin(t)*cos(t)/(1+sin(t)^2))
lines(lemniscate(t))

deltoid(t,n) = Chain.(2n*cos(t)*(1+cos(t))-n,2n*sin(t)*(1-cos(t)))
lines(deltoid(t,1))
lines!(deltoid(t,2))
lines!(deltoid(t,3))

astroid(t) = Chain.(cos(t)^3,sin(t)^3)
lines(astroid(t))

epitrochoid(t,n,c) = unitcircle(t)-c*unitcircle(n*t)
lines(epitrochoid(t,2,1))

# Exercise 1.13

toroidal(t,n) = Chain.((4+sin(n*t))*cos(t),(4+sin(n*t))*sin(t),cos(n*t))
lines(toroidal(t,1))

trefoil(t,n=1.5) = Chain.((2+cos(n*t))*cos(t),(2+cos(n*t))*sin(t),sin(n*t))
lines(trefoil(2t))

t = TensorField(-1:0.01:1)
lines(Chain.(t,t^2,t^3))

