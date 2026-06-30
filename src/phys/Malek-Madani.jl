# https://github.com/chakravala
#

# 1.1 Introduction to MATLAB

x = TensorField(-π:π/10:π)
lines(tan(sin(x)) - sin(tan(x)))

x = TensorField(-π:π/100:π)
lines(tan(sin(x)) - sin(tan(x)))

x = TensorField(-10π:π/100:10π)
lines(tan(sin(x)) - sin(tan(x)))

# 1.2

A = [1,2,3,-1]
B = [0,3,-2,7]
A+B
A.*B

A = [1 2; 3 -1]
B = [0 3; -2 7]
A+B

A = [1,3,4,6]
B = [4,6,-7,12]
A./B

A = [1,2,3,-1]
B = [2,3,3,2]
A.^B

A = [1,2,3,7]
A.^2-3A
sin.(A)
cos.(A)
exp.(A)
log.(A)

# 1.3

x = TensorField(LinRange(-2π,2π,1000))
f(a,b,c,x) = a*x^2+b*x+c
lines(f.(1,-2,3,x))

x = TensorField(LinRange(-6π,6π,1000))
f(a,b,c,x) = a*x^2+b*x+c
lines(f.(1,-2,3,x))

x = TensorField(0:0.1:4π)
lines(x*sin(2x))
lines(x*sin(x))

x = TensorField(0:0.1:1)
lines(sin(x))
scatter!(sinh(x))
scatter!(x^3)

x = TensorField(0:0.01:4π)
lines(x*exp(sin(x)))
lines(tangent(x*exp(sin(x))))
lines(tangent(tangent(x*exp(sin(x)))))

t = TensorField(0:0.01:2π)
lines(Chain.(sin(t)^2,cos(t)))

# 1.4

t = TensorField(0:0.01:10)
lines(Chain.(t,cos(t),sin(2t)))

XY = TensorField(ProductSpace(-5:0.5:5,-3:0.1:3))
wireframe(graph(besselj(3,norm(XY))))
contour(besselj(3,norm(XY)))

# 1.6

t = TensorField(1:100)
lines(cumsum(1/t^2))
sum(1/t^2)

t = TensorField(1:1000)
(sum(1/t^2),π^2/6)

t = TensorField(1:1000000)
(sum(1/t^2),π^2/6)

x = TensorField(0:0.01:10)
lines(besselj(1,x))
for i ∈ 2:10
    lines!(besselj(i,x))
end

xk = TensorField(ProductSpace(0:0.01:2π,1:16))
fun(x) = sin(x[2]*x[1])
variation(fun.(xk),0.1,lines,lines!)

# 1.7

x = TensorField(-2:0.1:2)
absvalue(x) = x[1]<0 ? -x : x
lines(absvalue.(x))

x = TensorField(-3:0.1:3)
vecabsval(x) = -x*(1-heaviside(x))+x*heaviside(x)
lines(vecabsval.(x))

# 1.8

function FourierSine(N,f)
    x = TensorField(base(f))
    L = fiber(x[end]-x[1]) # T/2
    ω = Cartan.FourierSpace((π/L)*(1:N),points(x))
    TensorField(ω,[(2/L)*integrate(f*sin(nω*x)) for nω ∈ ω])
end

x = TensorField(0:0.01:5)
b = FourierSine(8,1+x^2)

lines(1+x^2)
lines!(sum(b.*sin.(((π/5)*(1:8)).*Ref(x))))

lines(1+x^2)
for i ∈ 1:3
    lines!(sum(b.*sin.(((π/5)*(1:2^(2+i))).*Ref(x))))
end

b = FourierSine(2^(2+3),1+x^2)
lines(sum(b.*sin.(((π/5)*(1:2^(2+3))).*Ref(x))))

function FourierSine(N,M,f)
    xy = TensorField(base(f))
    x,y = getindex.(xy,1),getindex.(xy,2)
    L1 = points(x).v[1][end]-points(x).v[1][1]
    L2 = points(x).v[2][end]-points(x).v[2][1]
    ω1 = Cartan.FourierSpace((π/L1)*(1:N),points(x).v[1])
    ω2 = Cartan.FourierSpace((π/L2)*(1:N),points(x).v[1])
    TensorField(ProductSpace(ω1,ω2),
        [(4/(L1*L2))*integrate(f*sin(n*x)*sin(m*y)) for n ∈ ω1, m ∈ ω2])
end

xy = TensorField(ProductSpace(0:0.01:2,0:0.01:3))
x,y = getindex.(xy,1),getindex.(xy,2)
f = exp(-(1-x)^2*(1-y)^2)*sin(π*x)*sin(2π*y)+exp(-(1.5-x)^2*(2.3-y)^2)*sin(2π*x)*sin(π*y)
b = FourierSine(10,10,f)
S = sum([fiber(b)[n,m]*(sin((n*π/2)*x)*sin((m*π/3)*y)) for n ∈ 1:10, m ∈ 1:10])
contour(S)
maximum(abs(f-S))
contour(norm(f-S))

# 1.9

pendulum(α,ω,A) = x -> Chain(x[2],-α*x[2]-sin(x[1])-A*cos(ω*point(x)))
ic = IC(pendulum(0.1,3.2,14),Chain(0,1),30)
sol = odesolve(ic,ExplicitIntegrator{4}(2^-7))

# 2.3

A = [1 -2 3; 0 1 -1]
B = [0 1; -1 1; 1 3]
A*B

zeros(10,10)

A = spdiagm(-2=>[2,-7,3,-5])
B = spdiagm(1=>[2,2])

A = 2I + spdiagm(-1=>-ones(9),1=>-ones(9))

A = [1 2; 3 4]
B = A'

# 2.6

A = [1 1 0; 0 1 1; 1 0 1]
inv(A)

# 2.9

A = [1 1 1 1 6; 2 -3 1 1 -1; -3 0 1 1 0 ; 2 -1 0 -1 0]
rref(A)

rref([1 1 1; 1 1 0])
rref([1 1 1; 2 2 2])

# 2.10

t = TensorField(0:2π/100:2π)
lines(Chain.(cos(t),sin(t)))
lines!(Chain.(cos(t)+2sin(t),2cos(t)+sin(t)))

D,V = eigen([1 2; 2 1])

hilb(N) = [inv(i+j-1) for i ∈ 1:N, j ∈ 1:N]
A = hilb(6)
D,V = eigen(hilb(6))

# 3.2

x = TensorField(0:0.01:2π)
lines(sin(x))

z = 0*x
for i ∈ 1:4
    z += (-1)^(i-1)*x^(2i-1)/factorial(2i-1)
end
lines!(z)

# 3.3

XY = TensorField(ProductSpace(0:0.2:3,0:0.15:π))
x,y = getindex.(XY,1),getindex.(XY,2)
z = exp(-x)*sin(y)
streamplot(gradient(z))
contour!(z)

# 3.6

t = TensorField(LinRange(0,π/2,100))
integrate(-1/sqrt(cos(t)^2+2sin(t)^2))

uv = TensorField(ProductSpace(LinRange(0,2π,100),LinRange(0,π/2,100)))
fun(x) = cos(x[1])^2*sin(x[2])^2*cos(x[2])
integrate(fun.(uv))

# 4.2

roots(3,0.1,1)

x = TensorField(0:0.1:2π)
lines(exp(-0.05x)*(0.016625cos(1.7313x)+0.57664sin(1.7313x))-0.332502sin(3x)-0.0166251cos(3x))

# 4.4

fun(x) = -x[1]^2+point(x)
ic = IC(fun,0.1↦-0.3,4)
lines(odesolve(ic,MultistepIntegrator{4}(2^-7)))

rhs(x) = Chain(x[2],-0.1x[2]-sin(x[1])+cos(point(x)))
ic = IC(rhs,Chain(1,2),3)
sol = odesolve(ic,MultistepIntegrator{4}(2^-7))
lines(sol)
linegraph(sol)

ic = IC(rhs,Chain(1,2),30)
lines(odesolve(ic,MultistepIntegrator{4}(2^-7)))

# 4.5

odedef(α=0.3,β=1) = x -> Chain(x[2],-β*x[2]-α*sin(x[1]))
lines(odesolve(IC(odedef(),Chain(-5,0),30),MultistepIntegrator{4}(2^-7)))
for i ∈ -3:2:5
    lines!(odesolve(IC(odedef(),Chain(i,0),30),MultistepIntegrator{4}(2^-7)))
end

lines(odesolve(IC(odedef(0.3,0.1),Chain(-5,0),30),MultistepIntegrator{4}(2^-7)))
for i ∈ -3:2:5
    lines!(odesolve(IC(odedef(0.3,0.1),Chain(i,0),30),MultistepIntegrator{4}(2^-7)))
end

# 4.6

XY = TensorField(ProductSpace(-3:0.01:4,0:0.01:2.5))
function fpc(x)
    term = 1/(x[1]^2+x[2]^2)^2
    Chain(1-(x[1]^2-x[2]^2)*term,-2x[1]*x[2]*term)
end
streamplot(fpc.(XY))
lines(unitcircle(0:0.01:π))
lines!(Flow(fpc,1.0),[Chain(-2,0.2+0.3i)+0.1unitcircle(40) for i ∈ 0:6],5)

# 5.1



# 6.1

s = TensorField(0:2π/100:2π)
chi(x) = chi(point(x),x)
chi(t,x) = Chain(x[1]+t^2*x[2],x[2]-t*x[1])
lines(unitcircle(s))
for i ∈ 1:5
    lines!(chi.(0.3i,unitcircle(s)))
end

# 6.4

ψ(x) = (x[1]^2+x[2]^2)/2
XY = TensorField(ProductSpace(-3:0.01:3,-3:0.01:3))
contour(ψ.(XY))

ic = IC(-!gradient(ψ.(XY)),Chain(1.2,2.3),3)
lines!(odesolve(ic,MultistepIntegrator{4}(2^-7)))


