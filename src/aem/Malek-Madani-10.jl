# github.com/chakravala
# AEM: Malek-Madani
# 10 Integral Theorems of Vector Calculus

## 10.2 Double Integrals

# 10.2.1

D = OpenParameter(-1:0.1:1,2:0.1:3)
f(x) = x[1]^2-2x[2]
(integrate(f.(D)),-28/3)

# 10.2.2

D = OpenParameter(0:0.01:1,0:0.01:1)
g(x) = 0 < x[1] < 1-x[2] && 0 < x[2] < 1
surface(g.(D))
g(x) = 0 < x[2] < 1-x[1] && 0 < x[1] < 1
surface(g.(D))

D = OpenParameter(-1:0.01:1,-1:0.01:1)
g(x,a=1) = -sqrt(a^2-x[2]^2) < x[1] < sqrt(a^2-x[2]^2) && -a < x[2] < a
surface(g.(D))
g(x,a=1) = -sqrt(a^2-x[1]^2) < x[2] < sqrt(a^2-x[1]^2) && -a < x[1] < a
surface(g.(D))

# 10.2.3

D = OpenParameter(0:0.01:1,0:0.01:1)
g(x) = x[1]+x[2] < 1 && x[1]>0 && x[2]>0
f(x) = x[1]^2-x[2]^2
(integrate(g.(D)*f.(D)),0.0)

# 10.2.4

D = OpenParameter(0:0.01:1,-1:0.01:1)
g(x) = -sqrt(1-x[1]^2) < y < sqrt(1-x[1]^2) && 0 < x[1] < 1
f(x) = x[1]
(integrate(g.(D)*f.(D)),2/3)

# 10.2.5

D1 = OpenParameter(-1:0.01:1,-1:0.01:1)
g(x) = x[1]^2+x[2]^2 ≤ 1
surface(g.(D))

D2 = OpenParameter(LinRange(0,1,30),LinRange(0,2π,120))
h(x) = x[1]*Chain(cos(x[2]),sin(x[2]))
mesh(f.(D))

f1(x) = exp(x[1]^2+x[2]^2)
f2(x) = exp(x[1]^2)*x[1]
(integrate(g.(D1)*f1.(D1)),integrate(f2.(D2)),integrate(h.(D2),f1),π*(ℯ-1))

# 10.2.6

D1 = OpenParameter(-1.5:0.01:1.5,-1.5:0.01:1.5)
g(x) = 2x[1]^2+3x[2]^2 ≤ 4

D2 = OpenParameter(LinRange(0,1,30),LinRange(0,2π,120))
h(x) = Chain(sqrt(2)*x[1]*cos(x[2]),(2/sqrt(3))*x[1]*sin(x[2]))

f1(x) = x[1]*x[2]
f2(x) = (9/3)*x[1]^3*cos(x[2])*sin(x[2])

(integrate(g.(D1)*f1.(D1)),integrate(f2.(D2)),integrate(h.(D2),f1.(D1)))

centerofmass(ρ) = integrate(ρ*TensorField(base(ρ)))/integrate(ρ)
centerofmass(S,ρ) = integrate(S,ρ*TensorField(base(ρ)))/integrate(S,ρ)
centerofmass(g.(D1))

_momentofinertia(ρ,x) = ρ(x)*Chain(x[2]^2,x[1]^2)
momentofinertia(ρ) = integrate(_momentofinertia.(Ref(ρ),TensorField(base(ρ))))
momentofinertia(g.(D1))

# 10.3 Surface Integrals

# 10.3.1

X(x) = Chain(x[1]-x[2],x[2]-x[3],x[3]-x[1])
r(x) = Chain(x[1],x[2],0)
D = OpenParameter(-1:0.01:1,-1:0.01:1)
fluxintegrate(r.(D),X)

DD = OpenParameter(-1:0.01:1,-1:0.01:1,-1:0.01:1)
streamplot(X.(DD),gridsize=(11,11,11))
mesh!(r.(D))

# 10.3.2

D = OpenParameter(0:0.01:π/2,0:0.01:π/2)
r(x) = Chain(sin(x[1])*cos(x[2]),sin(x[1])*sin(x[2]),cos(x[1]))
X(x) = Chain(0,0,x[3])
(fluxintegrate(r.(D),X),π/6)

DD = OpenParameter(0:0.01:1,0:0.01:1,0:0.01:1)
streamplot(X.(DD),gridsize=(11,11,11))
mesh!(r.(D))

# 10.3.3

ψ(x) = x[2]-x[2]/norm(fiber(x))
r(x,x0) = Chain(x0,x[1]*cos(x[2]),x[1]*sin(x[2]))
S(a) = OpenParameter(0:0.01:a,0:0.01:2π)

DD = OpenParameter(-2:0.02:2,-2:0.02:2,-2:0.02:2)
streamplot(gradient(ψ.(DD)),gridsize=(11,11,11))
mesh!(r.(S(1),-1.5))
fluxintegrate(r.(S(1),-1.5),gradient(ψ.(DD)))

# 10.3.4

F(T,λ=1) = (-λ)*gradient(T)
T(x) = 300 - x[3]/50
r(x) = 500Chain(cos(x[1])*sin(x[2]),sin(x[1])*sin(x[2]),cos(x[2]))
D = OpenParameter(0:0.02:2π,0:0.02:π/2)
DD = OpenParameter(-500:10:500,-500:10:500,0:5:500)

fluxintegrate(r.(D),F(T.(DD)))

# 10.3.5

r(x,a=1) = a*Chain(cos(x[1])*sin(x[2]),sin(x[1])*sin(x[2]),cos(x[2]))
D = OpenParameter(0:0.01:2π,0:0.01:π)
(surfacearea(r.(D)),4π)

# 10.4

# 10.4.2

D = OpenParameter(0:0.01:1,0:0.01:1,0:0.01:1)
g(x) = x[1]+x[2]≤x[3]≤1 && 0≤x[2]≤1-x[1] && 0≤x[1]≤1
(integrate(g.(D)),1/6)

Makie.volume(g.(D),algorithm = :iso, isovalue = 1)

# 10.4.3

D = OpenParameter(-1:0.02:1,-1:0.02:1,0:0.02:1)
g(x) = x[3] ≥ x[1]^2 + x[2]^2 && x[3] ≤ 1
f(x) = x[3]
(integrate(g.(D)*f.(D)),π/3)


h(x) = Chain(x[1]*cos(x[2]),x[1]*sin(x[2]),x[3])
r(x) = x[3] ≥ x[1]^2
DD = OpenParameter(0:0.01:1,0:0.02:2π,0:0.01:1)
integrate(h.(DD),r.(DD)*f.(DD))

# 10.4.4

D = OpenParameter(-1:0.02:1,-1:0.02:1,0:0.01:1)
g(x) = x[1]^2+x[2]^2+(x[3]-1)^2≤1
f(x) = norm(fiber(x))
(integrate(g.(D)*f.(D)),8π/5)

h(x) = Chain(x[1]*cos(x[2])*sin(x[3]),x[1]*sin(x[2])*sin(x[3]),1+x[1]*cos(x[3]))
DD = OpenParameter(0:0.01:1,0:0.03:π,0:0.06:2π)
