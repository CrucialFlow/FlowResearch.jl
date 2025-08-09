# https://gitbub.com/chakravala

# Example 1.3

χ1(x) = Chain(x[1],x[2],0.0)/(1-x[3])
χ2(x) = Chain(x[1],x[2],0.0)/(1+x[3])

XYZ = TensorField(ProductSpace(-1.1:0.05:1.1,-1.1:0.05:1.1,-1.1:0.05:1.1))
fun(x) = x[1]^2+x[2]^2+x[3]^2
contour(fun.(XYZ))
contour(fun.(XYZ),levels=[1])

# Example 1.4

XY = TensorField(ProductSpace(-1.1:0.05:1.1,-1.1:0.05:1.1))
fun(x) = x[1]^2+x[2]^2
contour(fun.(XY))
contour(fun.(XY),levels=[1])

# Example 1.5

t = TorusParameter(100)
unitcircle(t)

# Example 1.6

tor = revolve(unitcircle()+Chain(sqrt(2),0))

fun(x) = 2sqrt(2(x[1]^2+x[2]^2))-(x[1]^2+x[2]^2+x[3]^2)
XYZ = TensorField(ProductSpace(-3:0.1:1.1,-3:0.1:3,-3:0.1:3))
Makie.volume(fun.(XYZ),algorithm=:iso,isovalue=1)

# Definition 1.9

t = TensorField(-1:0.01:1)
lines(Chain.(t^2,t^3))
lines(tangent(Chain.(t^2,t^3)))


