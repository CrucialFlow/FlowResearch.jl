# github.com/chakravala
# AEM: Malek-Madani

# 1.2

roots(3,-4,1)

roots(-1,1,2,3,-1,1)

t = TensorField(0:0.01:1)
integrate(1/(1+t^4))
lines(integral(1/(1+t^4)))

x = TensorFiel(0:0.01:3)
lines(integral(exp(-x^2)))
lines(integral(sin(x^2)))

t = TensorField(0:0.01:1)
h(x,n::Int=100) = h(x,TensorField(LinRange(0,1,n)))
h(x,t::TensorField) = integrate(1/sqrt(1-x[1]*sin(t)^2))
h(x::TensorField) = h.(x,Ref(x))
h(0.1)
lines(h(t))
tangent(h(t))(0.1)

