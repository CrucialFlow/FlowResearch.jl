# https://github.com/chakravala

# 1 Differential Equations

D = TensorField(ProductSpace(LinRange(0,2π,100),LinRange(0,2π,100)))
fun(x) = sin(x[1]^2+x[2])
variation(fun.(D),0.05,lines,lines!)


t = TensorField(LinRange(-1,1,1000))
SimpleEvolver(profile,tsize=0.1) = profile + tsize*profile*(profile - tangent(profile))

# 2.2

fun(k,x) = sin(k*x[1]+k*x[2]) + sin(k*x[1]-k*x[2])
fun(x) = fun(1,x) + 0.4fun(2,x) - 0.2fun(5,x)

# Burgers

t = TensorField(-2:0.01:8)
u = Chain.(t,1+0.5exp(-t^2))
fun(x) = Chain.(getindex.(localfiber(x),2),0)
variation(odesolve(IC(fun,u,10),ExplicitIntegrator{3}(0.1)),0.1,lines,lines!)

#

halfPeriods(28,-24)

WeierstrassP(z,k1,k2) = real(wp(t;omega=halfPeriods(k1,k2)))
WeierstrassPPrime(z,k1,k2) = imag(wp(t;omega=halfPeriods(k1,k2)))
WeierstrassInvariants(γ1,γ2) = ellipticInvariants(-γ2,γ1)

t = TensorField(-10:0.1:10)
lines(bound(weierstrassp(,70)))

