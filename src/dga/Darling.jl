# https://github.com/chakravala

# Figure 2.1

D = TensorField(ProductSpace(0:0.1:π,0:0.1:π))
fun(x) = Chain(sin(x[1]*x[2]),cos(x[1]*x[2]))
streamplot(fun.(D))

# Example 2.1.3

D = TensorField(ProductSpace(0:0.1:π,0:0.1:π,0:0.1:π))
fun(x) = Chain(x[1],x[2],x[3])/sqrt(x[1]^2+x[2]^2+x[3]^2)^3
streamplot(Lie[fun.(D),TensorField(D,Chain(0.0,1.0,0.0))],gridsize=(11,11,11))


