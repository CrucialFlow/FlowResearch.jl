# github.com/chakravala
# AEM: O'Neil
# 12 Vector Integral Calculus

# 12.1 Line Integrals

# Example 12.1

t = TorusParameter(100)
circ = Chain.(2cos(t),2sin(t),0t+4)

t = TensorField(ProductSpace(0:0.01:3π))
circ = Chain.(2cos(t),2sin(t),0t+4)

# Example 12.2

t = TensorField(1:0.01:2)
curv = Chain.(t^3,-t+0,t^2)
fun(x) = Chain(x[1],-x[2]*x[3],exp(x[3]))
(integrate(curv,fun),111/4+exp(4)-exp(1))

# Example 12.3

t = OpenParameter(100)
curv = Chain.(1-3t,0t+1,1+2t)
fun(x) = Chain(x[1]*x[2]*x[3],-cos(x[2]*x[3]),x[1]*x[3])
(integrate(curv,fun),3/2)

# Example 12.4

t = TensorField(-1:0.01:4)
curv = Chain.(t^2,t+0)
fun(x) = Chain(x[1]*x[2],-x[2]*sin(x[1]))
(integrate(curv,fun),410+cos(16)/2-cos(1)/2)

# Example 12.5

fun(x) = Chain(x[1]^2*x[2],x[2]^2)

t = TensorField(LinRange(0,π/2,100))
curv1 = Chain.(cos(t),sin(t))
(integrate(curv1,fun),1/3-π/16)

p = TensorField(0:0.01:2)
curv2 = Chain.(p+0,0p+1)
(integrate(curv2,fun),8/3)

(integrate(curv1,fun)+integrate(curv2,fun),3-π/16)

# Example 12.6

fun(x) = Chain(1.0,-x[2],x[1]*x[2]*x[3])
t = OpenParameter(100)
curv = Chain.(t+0,-t^2,t+0)
(integrate(curv,fun),3/10)

# Example 12.7

fun(x) = Chain(x[1]^2,-x[3]*x[2],x[1]*cos(x[3]))
t = TensorField(0:0.01:3)
curv = Chain.(t^2,t+0,π*t+0)
(-integrate(curv,fun),6/π+9π-243)

# Example 12.8

t = TensorField(LinRange(0,π/2,100))
curv = Chain.(4cos(t),4sin(t),0t-3)
fun(x) = x[1]*x[2]
D = TensorField(ProductSpace(-5:0.1:5,-5:0.1:5,-4:0.1:-2))
(integrate(curv,fun.(D)),32)

# Example 12.9

t = TensorField(LinRange(0,π/2,100))
curv = Chain.(2cos(t),2sin(t),0t+3)
fun(x) = x[1]*x[2]^2
D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3,2:0.1:4))
(integrate(curv,fun.(D)),16/3)

fun1(x) = x[1]*fun(x)
(integrate(curv,fun1.(D))/integrate(curv,fun.(D)),3π/8)

fun2(x) = x[2]*fun(x)
(integrate(curv,fun2.(D))/integrate(curv,fun.(D)),3/2)

fun3(x) = x[3]*fun(x)
(integrate(curv,fun3.(D))/integrate(curv,fun.(D)),3)

# 12.2 Green's Theorem

# Example 12.10

D = TensorField(ProductSpace(0:0.01:1,1:0.01:3))
V = Cartan.affmanifold(2)
fun(x) = Chain{V}(x[2]-x[1]^2*exp(x[1]),cos(2x[2]^2)-x[1])
fD = TensorField(D,fun.(fiber(D)))
gf1 = gradient(fD,1)
gf2 = gradient(fD,2)
(integrate(getindex.(gf1,2)-getindex.(gf2,1)),-4)

# Example 12.11

V = Cartan.affmanifold(2)
fun(x) = Chain{V}(2x[1]*cos(2x[2]),-2x[1]^2*sin(2x[2]))
D = OpenParameter(100,100)
fD = TensorField(D,fun.(fiber(D)))
gf1 = gradient(fD,1)
gf2 = gradient(fD,2)
(integrate(getindex.(gf1,2)-getindex.(gf2,1)),0)

# 12.3 Independence of Path and Potential Theory in the Plane

# Example 12.13

pot(x) = x[1]^2*cos(2x[2]) # + g(y)
fun(x) = Chain(2x[1]*cos(2x[2]),-2x[1]^2*sin(2x[2]))
t = OpenParameter(100)
curv = Chain.(t+0,t*(π/8)+0)
(integrate(curv,fun),sqrt(2)/2)

# Example 12.14

pot(x) = x[1]^2*cos(2x[2])
fun(x) = Chain(2x[1]*cos(2x[2]),-(2x[1]^2*sin(2x[2])+4x[2]^2))

# Example 12.15

pot(x) = x[1]*x[2]*exp(x[1]*x[2]) - x[2]^2 + x[1]^2
fun(x) = Chain(x[2]*exp(x[1]*x[2])+x[1]*x[2]^2*exp(x[1]*x[2])+2x[1],x[1]*exp(x[1]*x[2])+x[1]^2*x[2]*exp(x[1]*x[2])-2x[2])

# Example 12.16

fun(x) = Chain(2x[1]*x[2]^2+x[2],2x[1]^2*x[2]+exp(x[1])*x[2])

# Example 12.17

pot(x) = x[1]*x[2]*exp(x[1]*x[2]) - x[2]^2 + x[1]^2
fun(x) = Chain(x[2]*exp(x[1]*x[2])+x[1]*x[2]^2*exp(x[1]*x[2])+2x[1],x[1]*exp(x[1]*x[2])+x[1]^2*x[2]*exp(x[1]*x[2])-2x[2])

# Example 12.18

fun(x) = Chain(2x[1]*x[2]^2+x[2],2x[1]^2*x[2]+exp(x[1])*x[2])

# 12.3.1

fun(x) = Chain(-x[2],x[1])/(x[1]^2+x[2]^2)
t = TensorField(LinRange(0,π,100))
(integrate(unitcircle(t),fun),1π)
(integrate(unitcircle(-t),fun),-π)

# 12.4 Surfaces in 3-Space and Surface Integrals

# Example 12.19

t = TensorField(LinRange(0,2π,100))
curv = Chain.(2t+0,3t+0,t+0);
wireframe(revolve(curv))

# Example 12.20

x = TensorField(LinRange(0,2,100))
curv = Chain.(x+0,x+0,x^2/2)
t = TensorField(LinRange(0,2π,100))
wireframe(revolve(curv,Chain.(cos(t),sin(t),sin(2t))))
mesh(revolve(curv,Chain.(cos(t),sin(t),sin(2t))),normalnorm)

# Example 12.21

D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3,-1:0.1:3))
fun(x) = (a = 4-x[1]^2-x[2]^2; a<0 ? 0.0 : sqrt(a)-x[3])
contour(fun.(D),levels=[0.2])

D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3))
fun(x) = (a = 4-x[1]^2-x[2]^2; a<0 ? 0.0 : sqrt(a))
surface(fun.(D))

# Example 12.22

D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3,-1:0.1:3))
fun(x) = sqrt(x[1]^2+x[2]^2)-x[3]
contour(fun.(D))

D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3))
fun(x) = sqrt(x[1]^2+x[2]^2)
surface(fun.(D))

# Example 12.23

D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3,-1:0.1:3))
fun(x) = x[1]^2+x[2]^2-x[3]
contour(fun.(D))

D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3))
fun(x) = x[1]^2+x[2]^2
surface(fun.(D))

# Figures

D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3))

fun(x) = cos(x[1]^2-x[2]^2)
surface(fun.(D))

fun(x) = 6sin(x[1]-x[2])/sqrt(1+x[1]^2+x[2]^2)
surface(fun.(D))

fun(x) = 4cos(x[1]^2+x[2]^2)/(1+x[1]^2+x[2]^2)
surface(fun.(D))

fun(x) = x[1]^2*cos(x[1]^2-x[2]^2)
surface(fun.(D))

fun(x) = cos(x[1]*x[2])*log(4+x[2])
surface(fun.(D))

# Example 12.24

t = TensorField(LinRange(0,2π,100))
curv = Chain.(2t+0,3t+0,t+0);
normal(revolve(curv))(1/2,π/6)
Chain(-3sqrt(3)/4,-2/4,2*3/2)

# Example 12.25

D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3))
fun(x) = sqrt(x[1]^2+x[2]^2)
normal(fun.(D))(3,1)
-Chain(3/sqrt(10),1/sqrt(10),-1)

# Example 12.26

t = TensorField(LinRange(0,2π,100))
curv = Chain.(2t+0,3t+0,t+0);
jacobian(revolve(curv))(1/2,π/6)

# Example 12.27

t = TensorField(0:0.001:3)
wireframe(revolve(sqrt(9-t^2)))
(surfacearea(revolve(sqrt(9-t^2))),18π)

# Example 12.28

D = TensorField(ProductSpace(0:0.01:2,0:0.01:1))
fun(x) = 4-x[1]-x[2]
(surfacearea(graph(fun.(D))),2sqrt(3))

# Example 12.29

t = TensorField(2:0.01:3)
θ = TensorField(LinRange(0,π/2,100))
surf = revolve(t^2,unitcircle(θ))
wireframe(surf)
fun(x) = x[1]*x[2]/x[3]
(integrate(surf,fun),(sqrt(37)^3-sqrt(17)^3)/24)

# Example 12.30

x = TensorField(LinRange(1,2,100))
curv = Chain.(x+0,x+0,x^2/2)
t = TensorField(LinRange(0,π,100))
surf = revolve(curv,Chain.(cos(t),sin(t),sin(2t)))
wireframe(surf)
mesh(surf,normalnorm)
fun(x) = x[1]*x[2]*x[3]
(integrate(surf,fun),π/4*(100sqrt(5)/21-11sqrt(2)/105))

# 12.5 Applications of Surface Integrals

# Example 12.31

t = TensorField(0:0.01:2)
surf = revolve(t)
wireframe(surf)

fun(x) = x[1]^2+x[2]^2
(integrate(surf,fun),8sqrt(2)*π)

funz(x) = x[3]*(x[1]^2+x[2]^2)
(integrate(surf,funz)/integrate(surf,fun),8/5)

# Example 12.32

fun(x) = Chain(x[1],x[2],x[3])
t = TensorField(LinRange(0,sqrt(3),100))
surf = revolve(sqrt(4-t^2))
wireframe(surf)
(fluxintegrate(surf,fun),8π)

# 12.6 Preparation for the Integral Theorems of Gauss and Stokes

# 12.7 The Divergence Theorem of Gauss

# Example 12.33

fun(x) = Chain(x[1],x[2],x[3])
t = TensorField(0:0.01:1)
disk = revolve(Chain.(t+0,0t+1))
surf = revolve(t)
wireframe(disk)
wireframe!(surf)
(fluxintegrate(disk,fun)+fluxintegrate(surf,fun),1π)

# Example 12.34

D = OpenParameter(60,60,60);
fun(x) = Chain(x[1]^2,x[2]^2,x[3]^2)
(integrate(∂(fun.(D))),3)

# 12.8 The Integral Theorem of Stokes

# Example 12.35

fun(x) = Chain(-x[2],x[1],-x[1]*x[2]*x[3])
t = TorusParameter(100)
circ = Chain.(3cos(t),3sin(t),0t+3)
(integrate(circ,fun),18π)

t = TensorField(0:0.01:3)
surf = revolve(t)
streamplot(fun.(D))
wireframe!(surf)

D = TensorField(ProductSpace(-3:0.05:3,-3:0.05:3,0:0.05:3))
(fluxintegrate(surf,curl(fun.(D))),18π)

# Example 12.36

pot(x) = exp(x[1]*x[2]*x[3]) + x[2]*x[3] - 2x[1]^2 + sin(x[2])
fun(x) = Chain(x[2]*x[3]*exp(x[1]*x[2]*x[3])-4x[1],x[1]*x[3]*exp(x[1]*x[2]*x[3])+x[3]+cos(x[2]),x[1]*x[2]*exp(x[1]*x[2]*x[3])+x[2])

