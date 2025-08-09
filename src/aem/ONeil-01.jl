# github.com/chakravala
# AEM: O'Neil
# 01 First-Order Differential Equations

# 1.1 Preliminary Concepts

# Example 1.1

D = TensorField(ProductSpace(-1.5:0.1:6,-30:0.1:30))
fun(x) = Chain(1.0,2-x[2])
streamplot(fun.(D))

x = TensorField(-1.5:0.01:6)
gensol(k) = 2+k*exp(-x)

lines(gensol(0))
lines!(gensol(3))
lines!(gensol(6))
lines!(gensol(-3))
lines!(gensol(-6))

# Example 1.2

D = TensorField(ProductSpace(0.5:0.1:3,-25:0.1:40))
fun(x) = Chain(1.0,exp(x[1])-x[2]/x[1])

x = TensorField(0.5:0.01:3)
gensol(c) = ((x-1)*exp(x)+c)/x

lines(gensol(0))
lines!(gensol(5))
lines!(gensol(20))
lines!(gensol(-6))
lines!(gensol(-10))

# Example 1.3

D = TensorField(ProductSpace(-4:0.1:4,-15:0.1:15))
fun(x) = Chain(1.0,2-x[1]*x[2])
streamplot(fun.(D))

x = TensorField(0:0.01:4)
x2 = TensorField(-4:0.01:0)
gensol(k,x=x) = exp(-x^2/2)*(integral(2exp(x^2/2))+k)
gensol2(k) = gensol(-fiber(gensol(-k,x2))[end],x2)
bothlines(k) = (display(lines(gensol(k))); lines!(gensol2(k)))
bothlines!(k) = (display(lines!(gensol(k))); lines!(gensol2(k)))

bothlines(0)
bothlines!(4)
bothlines!(13)
bothlines!(-7)
bothlines!(-11)
bothlines!(-15)

# Example 1.4

fun(x) = 2-x[1]
ic = InitialCondition(fun,1.0↦-5.0,4.0)
sol = odesolve(ic,MultistepIntegrator{4}(2^-7))
lines(sol)

# Example 1.5

D = TensorField(ProductSpace(-4:0.1:4,-4:0.1:4))
fun(x) = Chain(1.0,x[2]^2)
streamplot(fun.(D))

# Example 1.6

D = TensorField(ProductSpace(-4:0.1:4,-4:0.1:4))
fun(x) = Chain(1.0,sin(x[1]*x[2]))
streamplot(fun.(D))





