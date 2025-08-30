# github.com/chakravala
# AEM: Malek-Madani
# 09 Laplace's Equation--Steady-State Flows and Normal Modes

# 9.2

bay(h=1,w=2) = TensorField(ProductSpace(LinRange(-w,0,100),LinRange(0,h,100)))
ψ(x,n,cn,h) = cn*sinh((π*n/h)*x[1])*sin((π*n/h)*x[2])
ψ1(x) = ψ(x,1,0.01,1)
contourf(-ψ1.(bay()),levels=LinRange(0,1,40))

initdata = [Chain(-2,0.001,0.002j) for j ∈ 0:5]
lines(!gradient(ψ1.(bay())),initdata)


ψ2(x) = ψ(x,2,0.001,π)
contourf(-ψ2.(bay(π)),levels=40)

# Exercise 5

!gradient(ψ.(bay(),1,-0.1,π))(-1.2,0.1)
!gradient(ψ.(bay(),1,0.1,π))(-1.2,0.1)
!gradient(ψ.(bay(),1,-0.1,π)+ψ.(bay(),2,0.2,π))(-1.2,0.1)
!gradient(ψ.(bay(),1,-0.1,π)+ψ.(bay(),2,-0.2,π))(-1.2,0.1)

# Exercise 20

V = Submanifold(2)
D = TensorField(ProductSpace{V}(-1:0.1:4,-1:0.1:7))
XY = TensorField(ProductSpace{V}(0.1:0.1:1,1:0.1:5))
ψ(x) = x[1]*x[2]

lines(Flow(-!gradient(ψ.(D)),1),boundarycomponents(XY),1)

# 9.3

# Exercise 12

xt(c=3,L=2,n=1) = TensorField(ProductSpace(LinRange(0,L,100),LinRange(0,2L/(n*c),100)))
fun(x,A=1.3,c=3) = A*sin((c*π/2)*x[2])*sin((π/2)*x[1])
variation(fun.(xt(),1.3,3),0.01,lines,lines!)
surface(fun.(xt(),1.3,3))
surface!(fun.(xt(),1.3,4.3))

# 9.4

# Example 9.4.2

fun(x,n,m) = sin(n*x[1])*sin(m*x[2]/sqrt(2))
D = TensorField(ProductSpace(-0.1:0.1:π+0.1,-0.1:0.1:sqrt(2)*π+0.1))
contour(fun.(D,1,2))
contour(fun.(D,1,12)+0.2fun.(D,12,1)-0.4fun.(D,8,9)+fun.(D,9.8),levels=[0])

# 9.5

# 9.5.2

f1(x) = ψ(Chain(sqrt(x[1]^2+x[2]^2),atan(x[2]/x[1])))

# 9.5.4

D = TensorField(ProductSpace(-3:0.1:3,-2:0.1:2))
ψ(x) = x[2]-x[2]/(x[1]^2+x[2]^2)
contour(ψ.(D),levels=200)

D = TensorField(ProductSpace(-4:0.1:4,-2:0.1:2))
initdata = [Chain(-4,sign(j)*sqrt(abs(0.1+0.2j))) for j ∈ -10:10]
lines(-!gradient(ψ.(D)),initdata,6)

D = TensorField(ProductSpace(-6:0.1:6,-2:0.1:2))
ic = InitialCondition(-!gradient(ψ.(D)),Chain(-6,1),10)
sol = odesolve(ic,MultistepIntegrator{4}(2^-7))
lines(sol)
lines!(Chain.(getindex.(sol,1)+0,speed(sol)))

fun(x) = -(1-2x[1]^2+2x[2]^2)/(2*(x[1]^2+x[2]^2))^2
lines(Chain.(getindex.(sol,1)+0,fun.(D)(sol)))

# 9.6

t = TensorField(0:0.01:15)
lines(besselj(2,t))

findroot(besselj(2,TensorField(11:0.0001:12)))

plate(x) = Chain(x[1]*cos(x[2]),x[1]*sin(x[2]),(0.1sin(2x[2])-0.4cos(2x[2]))*besselj(2,11.6198*x[1]))

mesh(plate.(base(unitdisk)))
wireframe!(plate.(base(unitdisk())),color=:black)

D = TensorField(ProductSpace(-1:0.01:1,-1:0.01:1))
contour(getindex.(plate.(vectorize(polarize(complexify(D)))),3))

# Project A

function wind_driven(;λ=1e9,b=π*2e8,d=2e4,F=1,R=0.02,coriolis=1e-13)
    D = TensorField(ProductSpace(LinRange(0,λ,100),LinRange(0,b,100)))
    γ = F*π/(R*b)
    μ = π/b
    α = d/R*gradient(coriolis*getindex.(D,2),2)
    k1 = -α/2 + sqrt(α^2/4+μ^2)
    k2 = -α/2 - sqrt(α^2/4+μ^2)
    p = (1-exp(k2*λ))/(exp(k1*λ) - exp(k2*λ))
    q = 1-p
    ψ(x) = (γ*(b/π)^2)*sin((π/b)*x[2])*(p(fiber(x))*exp(k1(fiber(x))*x[1])+q(fiber(x))*exp(k2(fiber(x))*x[1])-1)
    ψ.(D)
end

contourf(wind_driven();levels=30)

# Project B

D = TensorField(ProductSpace(0:0.01:3,0:0.01:3))
ψ(x) = 2sin(3x[1])*cos(4x[2])-0.5sin(4x[2])*cos(3x[1])-1.3cos(3x[1])*cos(4x[2])
contourf(ψ.(D))
streamplot(-!gradient(ψ.(D)))

ψ2(x) = ψ(x)-sin(5x[1])-cos(5x[2])
contourf(ψ2.(D))
streamplot(-!gradient(ψ2.(D)))

