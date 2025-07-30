# github.com/chakravala
# DGA: An Introductory Course on Differentiable Manifolds, Shahshahani

# 2.12 Smooth Bumb Functions

t = TensorField(0:0.1:5)
α(x) = x[1]>0 ? exp(-1/x[1]) : zero(x[1])
β(x,a,b) = α.(b-x)*α.(x-a)
γ = integral(β(t,1,4))/integrate(β(t,1,4))

lines(α.(t))
lines(β(t,1,4))
lines(1-γ.(abs(t-1)))

# Exercise 2.14

D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3))
f(x) = sin(x[1])*α(x)
fun(x) = Chain(x[1]*f(norm(x)),x[2]*f(norm(x)))

# Exercise 2.15

fun(x) = Chain(exp(x[1]*x[2]^2),0)

# Exercise 2.16

fun(x) = Chain((x[2]^2+1)*sin(x[1]),(x[1]^2+1)*sin(x[2]))

# Exercise 3.20

fun(x) = Chain(1+x[1]^2+x[2]^2,2x[1]*x[2])

# 5.8b

α(t) = Chain(t[1]^2-1,t[1]^3-t[1])

# Exercise 5.17

f(x) = Chain(g(x)^2,g(x)*h(x[2]))
g(x) = (a+b*cos(x[2]))*exp(x[1]*im)
h(x) = sin(x[1])+im*sin(2x[1])

# 6.26b

fun(x) = Chain(cos(x[1])^2,sin(x[1]))

