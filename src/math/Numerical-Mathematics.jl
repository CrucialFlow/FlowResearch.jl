# github.com/chakravala
# Math: Numerical Mathematics, Grasselli-Pelinovsky

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
wireframe(ruledsurface(c.(t),r.(t)))

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

# Green

fun(x) = Chain(x[1]^2*x[2]^5,x[1]^3*x[2]^2)
integrate(unitcircle(),fun)

XY = TensorField(ProductSpace(-1:0.1:1,-1:0.1:1))
disk = (x->float(abs(x)<1)).(XY) # compact support
gf1 = gradient(fun.(XY),1)
gf2 = gradient(fun.(XY),2)
integrate(disk*(getindex.(gf1,2)-getindex.(gf2,1)))

# Stokes

t = TensorField(LinRange(0,π,20))
rev = revolve(2Chain.(sin(t),cos(t)),unitcircle(TensorField(LinRange(-π/2,π/2,20))))
fun(x) = Chain(exp(x[1]*x[2]),x[2]*x[3]^2,exp(x[1]^2)*x[2])
t2 = TensorField(LinRange(0,2π,40))
circ = Chain.(0t2+0,2sin(t2),2cos(t2))

XYZ = TensorField(ProductSpace(-1:0.2:3,-3:0.2:3,-3:0.2:3))
streamplot(fun.(XYZ),gridsize=(11,17,17))
wireframe!(rev)

integrate(circ,fun.(XYZ))




