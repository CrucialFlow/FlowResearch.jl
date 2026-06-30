# https://github.com/chakravala
# math: Colley

# Chapter 3

# Chapter 3.1

# Example 3.1.1

t = TensorField(range(0,2pi,100))
lines(Chain(1,1,1)+t.*Chain(1,1,1))

# Example 3.1.2

t = TensorField(range(0,2pi,100))
lines(3unitcircle(t))

# Example 3.1.3

t = TensorField(range(0,2pi,100))
lines(3unithelix(t))

# Example 3.1.4

t = TensorField(range(0,2pi,100))
minimum(speed(unithelix(t)))

# Example 3.1.5

t = TensorField(range(0,2pi,100))
lines(Chain.(3t+2,t^2-7,t-t^2))

# Example 3.1.6

t = TensorField(range(0,2pi,100))
integrate((-9.8*v2)*t+50Chain(cos(π/4),sin(π/4)))

# Chapter 3.2

# Example 3.2.1

t = TensorField(range(0,2pi,100))
(totalarclength(unitcircle(t)),2pi)

# Example 3.2.2

t = TensorField(range(0,2pi,100))
(totalarclength(unithelix(t)),2pi*sqrt(2))


# Chapter 5

# Chapter 5.1

# Example 5.1.1

XY = TensorField(ProductSpace{2}(0:0.01:2,0:0.01:3))
(integrate(TensorField(XY,1.0)),2*3*1)
mesh(graph(TensorField(XY,1.0)))
mesh!(graph(TensorField(XY,0.0)))

# Example 5.1.2

XY = TensorField(ProductSpace{2}(-1:0.01:1,-1:0.01:1))
fun(x) = 4-x[1]^2-x[2]^2
(integrate(fun.(XY)),40/3)
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.1

XY = TensorField(ProductSpace{2}(0:0.01:2,1:0.01:3))
fun(x) = x[1]^2+x[2]
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.2

XY = TensorField(ProductSpace{2}(0:0.01:pi,1:0.01:2))
fun(x) = x[2]*sin(x[1])
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.3

XY = TensorField(ProductSpace{2}(-2:0.01:4,0:0.01:1))
fun(x) = x[1]*exp(x[2])
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.4

XY = TensorField(ProductSpace{2}(0:0.01:pi/2,0:0.01:1))
fun(x) = exp(x[1])*cos(x[2])
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.5

XY = TensorField(ProductSpace{2}(1:0.01:2,0:0.01:1))
fun(x) = exp(x[1]+x[2])+x[1]^2+log(x[2])
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.6

XY = TensorField(ProductSpace{2}(1:0.01:9,1:0.01:exp(1)))
fun(x) = log(sqrt(x[1]))/(x[1]*x[2])
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.7

XY = TensorField(ProductSpace{2}(-1:0.01:2,0:0.01:2))
fun(x) = x[1]^2+x[2]^2+2
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.8

XY = TensorField(ProductSpace{2}(0:0.01:3,1:0.01:2))
fun(x) = x[1]+3x[2]+1
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.9

XY = TensorField(ProductSpace{2}(0:0.01:1,-1:0.01:2))
fun(x) = 2x[1]^2+x[2]^4*sin(pi*x[1])
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.10

XY = TensorField(ProductSpace{2}(0:0.01:2,1:0.01:3))
integrate(TensorField(XY,2))
mesh(graph(TensorField(XY,2)))
mesh!(XY)

# Exercise 5.1.11

XY = TensorField(ProductSpace{2}(1:0.01:3,-2:0.01:2))
fun(x) = 16-x[1]^2-x[2]^2
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.12

XY = TensorField(ProductSpace{2}(-pi/2:0.01:pi/2,0:0.01:pi))
fun(x) = sin(x[1])*cos(x[2])
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.13

XY = TensorField(ProductSpace{2}(0:0.01:5,-2:0.01:2))
fun(x) = 4-x[1]^2
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.14

XY = TensorField(ProductSpace{2}(-2:0.01:3,0:0.01:1))
fun(x) = abs(x[1])*sin(pi*x[2])
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Exercise 5.1.15

XY = TensorField(ProductSpace{2}(-5:0.01:5,-1:0.01:2))
fun(x) = 5-abs(x[2])
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Chapter 5.2

# Example 5.2.1

XY = TensorField(ProductSpace{2}(-2:0.01:2,-1:0.01:3))
fun(x) = x[1]
integrate(fun.(XY))
surface(fun.(XY))
mesh!(XY)

# Chapter 5.4

# Example 5.4.1

XYZ = TensorField(ProductSpace{3}(-2:0.05:3,0:0.05:1,0:0.05:5))
fun(x) = x[1]^2*exp(x[2])+x[1]*x[2]*x[3]
(integrate(fun.(XYZ)),175*(exp(1)-1)/3+125/8)
contour(fun.(XYZ))




