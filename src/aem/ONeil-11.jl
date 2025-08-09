# github.com/chakravala
# AEM: O'Neil
# 11 Vector Differential Calculus

# 11.1 Vector Functions of One Variable

# Example 11.1

t = TensorField(LinRange(0,3π,100))
lines(tangent(Chain.(t^2,sin(t),-t^2)))

# Examples 11.2

t = TensorField(-4π:0.07:4π)
hlx = Chain.(cos(t),sin(t),t/3+0)
lines(hlx)
speed(hlx)
sqrt(10)/3
(totalarclength(hlx),8π/3*sqrt(10))

# Example 11.3

arclength(hlx)
speed(arcparametrize(hlx))
arcparemetrize(hlx)
speed(arcparametrize(hlx))

# 11.2 Velocity, Acceleration, Curvature, and Torsion

# Example 11.4

t = TensorField(LinRange(0,3π,100))
lines(Chain.(sin(t),2exp(-t),t^2))
lines(speed(Chain.(sin(t),2exp(-t),t^2)))
lines(sqrt(cos(t)^2+4exp(-2t)+4t^2))

# Example 11.6

t = TorusParameter(100)
circ = Chain.(4cos(t),0t+3,4sin(t))
curvature(circ)

# Example 11.7

t = TensorField(-4π:0.07:4π)
curv = Chain.(cos(t)+t*sin(t),sin(t)-t*cos(t),t^2)
lines(curv)

curvature(curv) # abs(1/(5t))

# Example 11.8

lines(arclength(curv))
lines!(sqrt(5)/2*t*abs(t))

# Example 11.9

lines(tangent(curv))
lines!(tangent(tangent(curv)))
lines!(unittangent(curv))
lines!(unitnormal(curv))

# Example 11.10

t = TensorField(0:0.01:3)
curv = Chain.(t^2,t^3,1t+0)
lines(curvature(curv))
lines(sqrt(36t^2+4+36t^4)/sqrt(4t^2+9t^4+1)^3)

# 11.3 Vector Fields and Streamlines

D = TensorField(ProductSpace(-4:0.1:4,-4:0.1:4))
G(x) = Chain(x[1]*x[2],x[1]-x[2])
H(x) = Chain(x[2]*cos(x[1]),x[1]^2-x[2]^2)
streamplot(G.(D))
streamplot(H.(D))

D = TensorField(ProductSpace(-4:0.1:4,-4:0.1:4,-4:0.1:4))
F(x) = Chain(cos(x[1]+x[2]),-x[1],x[1]-x[3])
Q(x) = Chain(-x[2],x[3],x[1]+x[2]+x[3])
M(x) = Chain(cos(x[1]),exp(-x[1])*sin(x[2]),x[3]-x[2])
streamplot(F.(D),gridsize=(11,11,11))
streamplot(Q.(D),gridsize=(11,11,11))
streamplot(M.(D),gridsize=(11,11,11))

# Example 11.11

D = TensorField(ProductSpace(-4:0.1:4,-4:0.1:4,-4:0.1:4))
fun(x) = Chain(x[1]^2,2x[2],-1.0)
streamplot(fun.(D),gridsize=(11,11,11))

# Example 11.12

D = TensorField(ProductSpace(-4:0.1:4,-4:0.1:4))
fun(x) = Chain(-x[1],x[2])
streamplot(fun.(D))

# 11.4 The Gradient Field and Directional Derivatives

D = TensorField(ProductSpace(-2π:0.3:2π,-2π:0.3:2π,-2π:0.3:2π))
fun(x) = x[1]^2*x[2]*cos(x[2]*x[3])
Makie.volume(fun.(D))
contour(fun.(D),levels=10,alpha=0.1)
streamplot(gradient(fun.(D)),gridsize=(11,11,11))

# Example 11.13

D = TensorField(ProductSpace(-2π:0.3:2π,-2π:0.3:2π,-2π:0.3:2π))
fun(x) = x[1]^2*x[2]-x[1]*exp(x[3])
Makie.volume(fun.(D))
streamplot(gradient(fun.(D)),gridsize=(11,11,11))

gradient(fun.(D))(2,-1,π)⋅Chain(1,-2,1)/sqrt(6)
-3/sqrt(6)*(4+exp(π))

# Example 11.14

D = TensorField(ProductSpace(-2π:0.3:2π,-2π:0.3:2π,-2π:0.3:2π))
fun(x) = 2x[1]*x[3]+exp(x[2])*x[3]^2
contour(fun.(D))
streamplot(gradient(fun.(D)),gridsize=(11,11,11))

abs(gradient(fun.(D))(2,1,1))
sqrt(4+ℯ^2+(4+2ℯ)^2)

# Example 11.15

D = TensorField(ProductSpace(-2:0.1:2,-2:0.1:2,-2:0.1:2))
fun(x) = x[3]-sqrt(x[1]^2+x[2]^2)
Makie.volume(fun.(D))
contour(fun.(D))
streamplot!(gradient(fun.(D)),gridsize=(11,11,11))

# Example 11.16

D = TensorField(ProductSpace(-2:0.1:2,-2:0.1:2,-2:0.1:2))
fun(x) = sin(x[1]*x[2])-x[3]
Makie.volume(fun.(D))
contour(fun.(D))
streamplot!(gradient(fun.(D)),gridsize=(11,11,11))

# Example 11.17

D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3,-3:0.1:3))
fun(x) = sin(x[1]^2+x[2]^2)-x[3]
Makie.volume(fun.(D))
contour(fun.(D),levels=2)
streamplot!(gradient(fun.(D)),gridsize=(11,11,11))

# Example 11.19

D = TensorField(ProductSpace(-2:0.1:2,-2:0.1:2,6:0.1:10))
fun(x) = x[1]^2+x[2]^2-x[3]
Makie.volume(fun.(D))
contour(fun.(D))
streamplot!(gradient(fun.(D)),gridsize=(11,11,11))

gradient(fun.(D))(2,-2,8)
abs(gradient(fun.(D))(2,-2,8))
sqrt(33)

# 11.5 Divergence and Curl

D = TensorField(ProductSpace(-2:0.1:2,-2:0.1:2,-2:0.1:2))
fun(x) = Chain(2x[1]*x[2],x[1]*x[2]*x[3]^2-sin(x[2]*x[3]),x[3]*exp(x[1]+x[2]))
contour(∂(fun.(D)))
streamplot!(fun.(D),gridsize=(11,11,11))

fun(x) = Chain(x[2],2x[1]*x[3],x[3]*exp(x[1]))
streamplot(fun.(D),gridsize=(11,11,11))
streamplot(curl(fun.(D)),gridsize=(11,11,11))

