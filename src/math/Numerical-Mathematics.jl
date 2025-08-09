# github.com/chakravala
# Math: Numerical Mathematics, Grasselli-Pelinovsky

# Figure 6.3

t = TensorField(LinRange(0,2π,50))
lines(tangent(sin(t)))
lines(tangent(tangent(sin(t))))

lines(tangent(sin(t))-cos(t))
lines(tangent(tangent(sin(t)))+sin(t))

# Figure 6.5

t = OpenParameter(100)
lines(sqrt(1-t^2))
integrate(sqrt(1-t^2))

# Figure 6.7

t = OpenParameter(100)
inv(sqrt(1-t^2))

t = TensorField(LinRange(0,π/2,100))
lines(cos(t)/sin(t))
integral(cos(t)/sin(t))

# Exercise 6.2

t = TensorField(0:0.01:2)
lines(log(1+t^2))
lines!(tangent(1+t^2))

# Exercise 6.3

t = TensorField(0:0.01:2)
lines(exp(-t^2))
lines!(tangent(exp(-t^2)))
lines!(tangent(tangent(exp(-t^2))))

# Exercise 6.4

t = OpenParameter(100)
lines(sin(t^3))
lines!(tangent(sin(t^3)))

# Exercise 6.5

t = TensorField(-1:0.0003:1)
lines(t*sin(1/t))
lines!(tangent(t*sin(1/t)))

# Exercise 6.6

t = TensorField(0:0.01:2)
lines(t*(1-t^2+2t^4-t^6))
lines!(tangent(t*(1-t^2+2t^4-t^6)))

# Exercise 6.7

t = TensorField(0:0.01:2)
lines(sech(t)^2)
lines!(tangent(sech(t)^2))
lines!(tangent(tangent(sech(t)^2)))

# Exercise 6.8

t = TensorField(LinRange(0,π/2,100))
lines(integral(cos(t)))

t = TensorField(1:0.01:2)
lines(integral(log(t)))

# Exercise 6.9

t = OpenParameter(100)
lines(integral(t*1-t^2))
t = OpenParameter(4)
lines!(integral(t*1-t^2))

# Exercise 6.10

t = OpenParameter(100)
lines(integral(inv(1+t^2)))

# Exercise 6.11

t = OpenParameter(100)
lines(integral(sqrt(t)^3))
lines!(integral(sqrt(t)^5))

# Exercise 6.12

t = OpenParameter(100)
lines(integral(cos(10π*t)))

# Exercise 6.13

t = OpenParameter(100)
lines(integral(inv(1+t)))

# Exercise 6.16

t = TensorField(LineRange(0,π/2,100))
lines(integral(sqrt(1+tangent(sin(t))^2)))
lines(arclength(graph(sin(t))))

# Exercise 6.18

t = TensorField(LinRange(0,π,100))
lines(integral(sin(t)^3))

# Exercise 6.19

t = TensorField(LinRange(0,π,100))
J0(x) = integrate(cos(x[1]*sin(t)))/π

x = TensorField(0:0.01:10)
J0x = TensorField(x,float(fiber(J0.(x))))
lines(J0x)

# Exercise 6.21

t = OpenParameter(100)
lines(integral(t^20))
(integrate(t^20),1/21)

# Exercise 6.22

t = TensorField(1:0.01:3)
mylog(x,n=100) = integrate(inv(TensorField(LinRange(1,x[1],n))))
mylogx = TensorField(t,float.(fiber(mylog.(t))))
lines(mylogx)
lines!(log(t))

# 7.1 function graph

t = TensorField(-1:0.01:1)
lines(graph(cbrt(t^2)-t^2))

# 7.2 surface function

XY = TensorField(ProductSpace(-4:0.5:4,-4:0.5:4))
surface(abs2(XY))
contour(abs2(XY))
contour3d(abs2(XY))

# 7.3 Cobb-Douglas production function

CobbDouglas(a=1,α=0.4,β=0.6) = x->a*x[1]^α*x[2]^β
XY = TensorField(ProductSpace(0:10:200,0:10:200))
surface(CobbDouglas().(XY))
contour(CobbDouglas().(XY))

# 7.4 volume

XYZ = TensorField(ProductSpace(-2:0.2:2,-2:0.2:2,-2:0.2:2))
fun(x) = x[1]^2+x[2]^2-x[3]^2
contour(fun.(XYZ),alpha=0.05)

# volumeslice plots todo

# ezsurf

XY = TensorField(ProductSpace(-2π:0.03:2π,-2π:0.03:2π))
fun(x) = x[1]^2/abs2(x)
surface(fun.(XY))

# 7.2 Partial Derivatives

XY = TensorField(ProductSpace(-2:0.5:2,-2:0.5:2))
fun(x) = x[1]^2-x[1]*x[2]-x[2]^2
surface(fun.(XY))

XY2 = TensorField(ProductSpace(-2:0.5:2,-5:5))
yplane(x) = Chain(x[1],-1,x[2])
xplane(x) = Chain(1,x[1],x[2])
mesh!(yplane.(XY2))
mesh!(xplane.(XY2))

XY = TensorField(ProductSpace(-1:0.2:1,-1:0.2:1))
fun(x) = sin(π*x[1])*x[2]^2
surface(fun.(XY))
surface(gradient(fun.(XY),1))
surface(gradient(fun.(XY),2))
surface(gradient(gradient(fun.(XY),1),1))
surface(gradient(gradient(fun.(XY),2),2))
surface(gradient(gradient(fun.(XY),1),2))
surface(gradient(gradient(fun.(XY),2),1))

XY = TensorField(ProductSpace(-0.5:0.01:0.5,-0.5:0.01:0.5))
fun(x) = iszero(x[1]*x[2]) ? 0.0 : x[1]*x[2]/sqrt(x[1]^2+x[2]^2)
surface(fun.(XY))
surface(gradient(fun.(XY),2))

# 7.3 The Gradient Vector

XY = TensorField(ProductSpace(-4:0.5:4,-4:0.5:4))
fun(x) = x[1]^2-x[2]^2
arrows(0.05gradient(fun.(XY)))
contour!(fun.(XY))

XYZ = TensorField(ProductSpace(-2:0.2:2,-2:0.2:2,-2:0.2:2))
fun(x) = x[1]^2+x[2]^2+x[3]^2
contour(fun.(XYZ),alpha=0.05)
arrows!(0.05gradient(fun.(XYZ)))

# 7.4 Paths

circ = unitcircle()
lines(circ)
arrows!(circ,tangent(circ);gridsize=13)

t = TensorField(0:0.1:4π)
hlx = unithelix(t)
lines(hlx)
arrows!(hlx,tangent(hlx);gridsize=26)

# 7.5 Vector Fields

XY = TensorField(ProductSpace(-4:0.5:4,-4:0.5:4))
fun(x) = Chain(-x[2],x[1])
arrows(fun.(XY))
streamplot!(fun.(XY))

# Figure 7.16

fun(x) = Chain(x[1]*x[2],x[1]^2)
surface(∂(fun.(XY)))
arrows!(0.05fun.(XY))

fun(x) = Chain(-x[1],-x[2])
mesh(graph(∂(fun.(XY))))
arrows!(0.05fun.(XY))

# Figure 7.17

XYZ = TensorField(ProductSpace(-4:4,-4:4,-4:4))
fun(x) = Chain(-x[2],x[1],0)
streamplot(fun.(XYZ),gridsize=(11,11,11))
streamplot(curl(fun.(XYZ)),gridsize=(11,11,11))

# Figure 7.18

XY = TensorField(ProductSpace(-2:0.2:2,-2:0.2:2))
arrows(-0.05*!XY/abs2(XY))

# 7.6 Line Integrals

fun(x) = 4 + 2x[2]*sin(x[1])

t = TensorField(0:0.05:3π)
Chain.(t,cos(t),4+2sin
graph(cos(t))

t = TensorField(0:0.05:2π)

cardiod(t) = Chain.((sin(t)-1)*cos(t),(sin(t)-1)*sin(t))
lines(cardiod(t))

limacon(t) = Chain.((3sin(t)+2)*cos(t),(3sin(t)+2)*sin(t))
lines(limacon(t))

trefoil(t) = Chain.(4cos(2t)+2cos(t),4sin(2t)-2sin(t),sin(3t))
lines(trefoil(t))

# 7.7 Surface Integrals

# 7.36

t = TorusParameter(60)
mycirc = Chain.(2cos(t)+5,2sin(t))
wireframe(revolve(mycirc))

# 7.40

t = TorusParameter(60)
c(u) = 2Chain(cos(u[1]),sin(u[1]),0)
r(u) = Chain(cos(u[1]/2)*cos(u[1]),cos(u[1]/2)*sin(u[1]),sin(u[1]/2))
mobi = MobiusTopology(ruledsurface(c.(t),r.(t)))
mesh(mobi,normalnorm)
wireframe!(mobi)

# asphere

t = TensorField(LinRange(0,1π,60))
circ = 2Chain.(sin(t),cos(t))
rev = revolve(circ)
surfacearea(revolve(circ))

# msphere

mass(x) = 4*(cos(x[1]))^2
integrate(rev,mass.(base(rev)))

# flux integral

field(q=10,ε0=1) = x-> (q/ε0/4π)*fiber(x)/abs(fiber(x))^3
fluxintegrate(rev,field().(rev))

# 7.8

t = TensorField(LinRange(0,2π,100))
f1 = Chain.(1t+0,sin(t),0t+0)
f2 = Chain.(1t+0,-sin(t),0t+0)
wireframe(scrollsurface(f1,f2))

t = TensorField(-1:0.1:1)
g1 = Chain.(0t+0,1t+0,0t+0)
g2 = Chain.(t^3-t+1,1t+0,0t+0)
wireframe(scrollsurface(g1,g2))

# Green

fun(x) = Chain(x[1]^2*x[2]^5,x[1]^3*x[2]^2)
integrate(unitcircle(),fun)

XY = TensorField(ProductSpace(-1:0.01:1,-1:0.01:1))
disk = (x->float(abs(x)<1)).(XY) # compact support
gf1 = gradient(fun.(XY),1)
gf2 = gradient(fun.(XY),2)
integrate(disk*(getindex.(gf1,2)-getindex.(gf2,1)))

# Stokes

t = TensorField(LinRange(0,π,20))
rev = revolve(2Chain.(sin(t),cos(t)),unitcircle(TensorField(LinRange(-π/2,π/2,100))))
fun(x) = Chain(exp(x[1]*x[2]),x[2]*x[3]^2,exp(x[1]^2)*x[2])
t2 = TensorField(LinRange(0,2π,40))
circ = Chain.(0t2+0,2sin(t2),2cos(t2))

XYZ = TensorField(ProductSpace(-1:0.01:3,-3:0.03:3,-3:0.03:3))
streamplot(fun.(XYZ),gridsize=(11,17,17))
wireframe!(rev)

integrate(circ,fun.(XYZ))
fluxintegrate(rev,curl(fun.(XYZ))

# Gauss

XYZ = TensorField(ProductSpace(0:0.01:1,0:0.02:2,-1:0.02:1))
fun(x) = Chain(x[1]^4*exp(x[2]),exp(x[1])*cos(x[2]),sin(x[1])+x[2]^2*x[3])

integrate(∂(fun.(XYZ)))

boundarycomponents(fun.(XYZ))
boundarycomponents(XYZ)

#fluxintegrate(X1,F1)

