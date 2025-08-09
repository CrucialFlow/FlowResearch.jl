# github.com/chakravala
# AEM: O'Neil
# 10 Qualitative Methods and Systems of Nonlinear Differential Equations

# Nonlinear Systems and Existence of Solutions

# Example 10.1

D = TensorField(ProductSpace(-15:0.1:15,-7:0.1:9))
fun(x,ω=sqrt(10),γ=0.3) = Chain(x[2],-ω^2*sin(x[1])-γ*x[2])
streamplot(fun.(D))

# Example 10.2

D = TensorField(ProductSpace(-9:0.1:11,-40:0.1:30))
fun(x,α=0.2,k=4,c=2,m=1) = Chain(x[2],(-k*x[1]+α*x[1]^3-c*x[2])/m)
streamplot(fun.(D))

# The Phase Plane, Phase Portraits, and Direction Fields

# Example 10.3

D = TensorField(ProductSpace(-3:0.1:2,-6:0.1:6))
fun(x) = Chain(x[2],x[1]^2*x[2]^2)
streamplot(fun.(D))

# Example 10.4

D = TensorField(ProductSpace(-2:0.1:2,-2:0.1:2))
fun(x) = Chain(-2x[2]-x[1]*sin(x[1]*x[2]),2x[1]+x[2]*sin(x[1]*x[2]))
streamplot(fun.(D))

# Figure

D = TensorField(ProductSpace(0:0.1:2.3,-3:0.1:2))
fun(x) = Chain(x[1]*cos(x[2]),x[1]^2-x[2]^3+sin(x[1]+x[2]))
streamplot(fun.(D))

# Example 10.5

fun(x) = Chain(x[1]/point(x),x[1]-x[2]/point(x))

# 10.3 Phase Portraits of Linear Systems

V = Cartan.affmanifold(2)

# Example 10.6

D = TensorField(ProductSpace(-3:0.1:2,-4.5:0.1:4))
C = Endomorphism{V}([-6 -2; 5 1])
eigvals(C)
eigvecs(C)
streamplot(C⋅D)

# Example 10.7

D = TensorField(ProductSpace(-3:0.1:2,-4.5:0.1:4))
C = Endomorphism{V}([3 3; 1 5])
eigvals(C)
eigvecs(C)
streamplot(C⋅D)

# Example 10.8

D = TensorField(ProductSpace(-5:0.1:5,-5:0.1:5))
C = Endomorphism{V}([-1 3; 2 -2])
eigvals(C)
eigvecs(C)
streamplot(C⋅D)

D = TensorField(ProductSpace(-160:1:50,-30:0.1:200))
streamplot(C⋅D)

# Example 10.9

D = TensorField(ProductSpace(-100:1:110,-15:1:60))
C = Endomorphism{V}([-10 6; -6 2])
eigvals(C)
eigvecs(C)
streamplot(C⋅D)

# Example 10.10

D = TensorField(ProductSpace(-40:1:20,-30:1:50))
C = Endomorphism{V}([-1 -2; 4 3])
eigvals(C)
eigvecs(C)
streamplot(C⋅D)

# Example 10.11

D = TensorField(ProductSpace(-50:1:50,-11:0.1:11))
C = Endomorphism{V}([3 18; -1 -3])
eigvals(C)
eigvecs(C)
streamplot(C⋅D)

# 10.4 Critical Points and Stability

# Example 10.12

D = TensorField(ProductSpace(-15:0.1:15,-7:0.1:9))
fun(x,ω=sqrt(10),γ=0.3) = Chain(x[2],-ω^2*sin(x[1])-γ*x[2])
streamplot(fun.(D))

# Example 10.13

D = TensorField(ProductSpace(-9:0.1:11,-40:0.1:30))
fun(x,α=0.2,k=4,c=2,m=1) = Chain(x[2],(-k*x[1]+α*x[1]^3-c*x[2])/m)
streamplot(fun.(D))

# 10.5 Almost Linear Systems

# Example 10.14

D = TensorField(ProductSpace(-1:0.1:5,-3:0.1:6))
fun(x) = Chain(4x[1]-2x[2]-4x[1]*x[2],x[1]+6x[2]-8x[1]^2*x[2])
streamplot(fun.(D))
jacobian(fun.(D))(0,0)

D = TensorField(ProductSpace(-6:0.1:8,-33:1:30))
fun(x) = Chain(4x[1]-2x[2],x[1]+6x[2])
streamplot(fun.(D))

# Example 10.15

D = TensorField(ProductSpace(-0.6:0.01:0.5,-2:0.01:7))
fun(x) = Chain(-x[1]-x[2]+x[1]^2*x[2]^2,-x[1]-3x[2]+x[1]^2-x[2]^2)
streamplot(fun.(D))
jacobian(fun.(D))(0,0)
eigvals(jacobian(fun.(D))(0,0))

D = TensorField(ProductSpace(-4:0.1:5,-22:1:30))
fun(x) = Chain(-x[1]-x[2],-x[1]-3x[2])
streamplot(fun.(D))

# Example 10.16

D = TensorField(ProductSpace(-1.1:0.01:1,-5:0.1:5))
fun(x) = Chain(3x[1]-4x[2]+x[1]^2*cos(x[2]),6x[1]+2x[2]-x[2]^3)
streamplot(fun.(D))
jacobian(fun.(D))(0,0)
eigvals(jacobian(fun.(D))(0,0))

D = TensorField(ProductSpace(-2.1:0.01:2.1,-5:0.1:5))
fun(x) = Chain(3x[1]-4x[2],6x[1]+2x[2])
streamplot(fun.(D))

# Example 10.17

D = TensorField(ProductSpace(-75:1:75,-33:1:33))
fun(x) = Chain(-x[1]+2x[2]+x[1]*sin(x[2]),2x[1]+3x[2]+8sin(x[1]))
streamplot(fun.(D))
jacobian(fun.(D))(0,0)
eigvals(jacobian(fun.(D))(0,0))

D = TensorField(ProductSpace(-80:1:65,-30:1:35))
fun(x) = Chain(-x[1]+2x[2],2x[1]+3x[2])
streamplot(fun.(D))

# Example 10.18

D = TensorField(ProductSpace(-33:1:10,-35:1:10))
fun(x) = Chain(4x[1]+11x[2]+x[1]*sin(x[2]),-2x[1]-4x[2]+sin(x[2]))
streamplot(fun.(D))
jacobian(fun.(D))(0,0)
eigvals(jacobian(fun.(D))(0,0))

D = TensorField(ProductSpace(-20:1:20,-6:0.1:6))
fun(x) = Chain(4x[1]+11x[2],-2x[1]-4x[2])
streamplot(fun.(D))

# Example 10.19

D = TensorField(ProductSpace(-0.8:0.01:0.8,-8:0.1:6))
fun(x,α=1/3,h=0.4,k=0.7) = Chain(α*x[1]-x[2]+h*x[1]*(x[1]^2+x[2]^2),x[1]-α*x[2]+k*x[2]*(x[1]^2+x[2]^2))
streamplot(fun.(D,1/3))
jacobian(fun.(D))(0,0)
eigvals(jacobian(fun.(D))(0,0))

D = TensorField(ProductSpace(-0.5:0.01:0.5,-8:0.1:6))
streamplot(fun.(D,2))
jacobian(fun.(D,2))(0,0)
eigvals(jacobian(fun.(D,2))(0,0))

# Example 10.20

D = TensorField(ProductSpace(-1.5:0.01:1.5,-5:0.1:9))
fun(x,ϵ=-0.2) = Chain(x[2]+ϵ*x[1]*(x[1]^2+x[2]^2),-x[1]+ϵ*x[2]*(x[1]^2+x[2]^2))
streamplot(fun.(D,-0.2))
streamplot(fun.(D,0.2))
jacobian(fun.(D))(0,0)
eigvals(jacobian(fun.(D))(0,0))

D = TensorField(ProductSpace(-9:0.01:9,-9:0.1:9))
streamplot(fun.(D,0))

# Example 10.21

D = TensorField(ProductSpace(-3:0.01:1,-1.5:0.01:3))
fun(x) = Chain(sin(π*x[1])-x[1]^2+x[2]^2,cos((x[1]+x[2]+1)*π/2))
streamplot(fun.(D))

D = TensorField(ProductSpace(-2:0.1:8,-2.2:0.01:2.2))
streamplot(fun.(D))

D = TensorField(ProductSpace(1.4:0.1:4.6,1:0.1:5))
streamplot(fun.(D))

# Example 10.22

D = TensorField(ProductSpace(-3.3:0.1:3.3,-4.2:0.1:4.2))
fun(x,ω=sqrt(10),γ=0.3) = Chain(x[2],-ω^2*sin(x[1])-γ*x[2])

streamplot(fun.(D,0.44,sqrt(0.8)))
jacobian(fun.(D,0.44,sqrt(0.8)))(0,0)
eigvals(jacobian(fun.(D,0.44,sqrt(0.8)))(0,0))

streamplot(fun.(D,sqrt(0.2),sqrt(0.8)))
jacobian(fun.(D,sqrt(0.2),sqrt(0.8)))(0,0)
eigvals(jacobian(fun.(D,sqrt(0.2),sqrt(0.8)))(0,0))

D = TensorField(ProductSpace(-1.5:0.01:1.5,-3:0.1:3))
streamplot(fun.(D,sqrt(0.3),sqrt(0.6)))
jacobian(fun.(D,sqrt(0.3),sqrt(0.6)))(0,0)
eigvals(jacobian(fun.(D,sqrt(0.3),sqrt(0.6)))(0,0))
jacobian(fun.(D,sqrt(0.3),sqrt(0.6)))(π,0)
eigvals(jacobian(fun.(D,sqrt(0.3),sqrt(0.6)))(π,0))

# Example 10.23

D = TensorField(ProductSpace(-2.2:0.1:2,-2:0.1:3))
fun(x,α=0.2,k=4,c=2,m=1) = Chain(x[2],(-k*x[1]+α*x[1]^3-c*x[2])/m)
streamplot(fun.(D,1,5,2,3))
jac = jacobian(fun.(D,1,5,2,3))
jac(0,0)
eigvals(jac(0,0))
jac(sqrt(5/1),0)
eigvals(jac(sqrt(5/1),0))

# Example 10.24

fun(x,k,K,L,m) = Chain(x[2],(-k*x[1]+K*(inv(L-x[1])^2-inv(L^2)))/m)

# 10.6 Predator/Prey Population Models

# Example 10.25

D = TensorField(ProductSpace(0:1:160,0:1:70))
fun(x,a,b,k,c) = Chain(a*x[1]-b*x[1]*x[2],c*x[1]*x[2]-k*x[2])
streamplot(fun.(D,0.2,0.02,0.02,1.2))
jacobian(fun.(D,0.2,0.02,0.02,1.2))(0,0)
eigvals(jacobian(fun.(D,0.2,0.02,0.02,1.2))(0,0))
jacobian(fun.(D,0.2,0.02,0.02,1.2))(60,10)
eigvals(jacobian(fun.(D,0.2,0.02,0.02,1.2))(60,10))

# Example 10.26

D = TensorField(ProductSpace(0:0.1:6,0:0.1:6))
fun(x) = Chain(x[1]-(x[1]^2)/4-x[1]*x[2]/4,x[1]*x[2]-x[2]-(x[2]^2)/2)
streamplot(fun.(D))
jacobian(fun.(D))(0,0)
eigvals(jacobian(fun.(D))(0,0))
jacobian(fun.(D))(4,0)
eigvals(jacobian(fun.(D))(4,0))
jacobian(fun.(D))(2,2)
eigvals(jacobian(fun.(D))(2,2))

# 10.7 Competing Species Models

# Example 10.27

D = TensorField(ProductSpace(0:0.1:9,0:0.1:7))
fun(x,a,b,k,c) = Chain(a*x[1]-b*x[1]*x[2],k*x[2]-c*x[1]*x[2])
streamplot(fun.(D,2,0.3,4,0.7))
jac = jacobian(fun.(D,2,0.3,4,0.7))
jac(0,0)
eigvals(jac(0,0))
jac(40/7,20/3)
eigvals(jac(40/7,20/3))

D = TensorField(ProductSpace(0:1:33,0:1:30))
streamplot(fun.(D,2,0.3,4,0.7))

# Example 10.28

fun(x,G,D,C,g,d,c) = Chain(x[1]*(G-D*x[1]-C*x[2]),x[2]*(g-d*x[2]-c*x[1]))
D = TensorField(ProductSpace(0:0.1:4,0:0.1:5))
streamplot(fun.(D,3/5,1/6,1/4,1,1/4,1/2))
contour!(getindex.(fun.(D,3/5,1/6,1/4,1,1/4,1/2),1),levels=[0])
contour!(getindex.(fun.(D,3/5,1/6,1/4,1,1/4,1/2),2),levels=[0])
jac = jacobian(fun.(D,3/5,1/6,1/4,1,1/4,1/2))
jac(0,0)
eigvals(jac(0,0))
jac(0,4)
eigvals(jac(0,4))
jac(18/5,0)
eigvals(jac(18/5,0))
jac(6/5,8/5)
eigvals(jac(6/5,8/5))

# Example 10.29

D = TensorField(ProductSpace(0.5:0.1:4,1:0.1:5))
streamplot(fun.(D,3,1,1/4,2,1/2,1/6))
contour!(getindex.(fun.(D,3,1,1/4,2,1/2,1/6),1),levels=[0])
contour!(getindex.(fun.(D,3,1,1/4,2,1/2,1/6),2),levels=[0])
jac = jacobian(fun.(D,3,1,1/4,2,1/2,1/6))
jac(0,0)
eigvals(jac(0,0))
jac(0,4)
eigvals(jac(0,4))
jac(3,0)
eigvals(jac(3,0))
jac(24/11,36/11)
eigvals(jac(24/11,36/11))

# 10.8 Lyapunov's Stability Criteria

# Example 10.30

D = TensorField(ProductSpace(-1:0.01:1,-2:0.01:2))
fun(x) = Chain(-x[1]^3,-4x[1]^2*x[2])
streamplot(fun.(D))

# Example 10.31

fun(x,g,L) = Chain(x[2],(-g/L)*sin(x[1]))

# Example 10.32

D = TensorField(ProductSpace(-1.1:0.01:1,-1:0.01:2))
fun(x) = Chain(x[1]^3-x[1]*x[2]^2,x[2]^3+6x[1]^2*x[2])
streamplot(fun.(D))

# Example 10.33

D = TensorField(ProductSpace(-14:0.1:10,-65:1:45))
fun(x,α=1,β=1/4,γ=1/6) = Chain(x[2],-α*x[2]-x[1]-β*x[1]^2-γ*x[1]^3)
streamplot(fun.(D))

# Example 10.34

function fun(x,fg)
    fgx = fg(point(x))
    f,g = fgx[1],fgx[2]
    Chain(f*x[2]+g*x[1]*(x[1]^2+x[2]^2),-f*x[1]+g*x[2]*(x[1]^2+x[2]^2))
end

# 10.9 Limit Cycles and Periodic Solutions

# Example 10.35

D = TensorField(ProductSpace(-11:0.1:11,-11:0.1:11))
fun(x) = Chain(x[1]+x[2]-x[1]*sqrt(x[1]^2+x[2]^2),-x[1]+x[2]-x[2]*sqrt(x[1]^2+x[2]^2))
streamplot(fun.(D))
jacobian(fun.(D))(0,0)
eigvals(jacobian(fun.(D))(0,0))

D = TensorField(ProductSpace(-3:0.01:3,-4:0.01:4))
streamplot(fun.(D))

# Example 10.36

D = TensorField(ProductSpace(-15:0.1:15,-15:0.1:15))
fun(x) = Chain(x[2]-x[1]*sin(sqrt(x[1]^2+x[2]^2)),-x[1]-x[2]*sin(sqrt(x[1]^2+x[2]^2)))
streamplot(fun.(D))

# Example 10.37

D = TensorField(ProductSpace(-4:0.01:4,-4:0.01:4))
fun(x) = Chain(3x[1]+4x[2]-x[1]^3,5x[1]-2x[2]+x[2]^3)
streamplot(fun.(D))

# Lienard

function fun(x,pq)
    pqx = pq(x)
    Chain(x[2],-pqx[1]*x[2]-pqx[2])
end

# Example 10.38

fun(x,α) = Chain(x[2],-α*(x[1]^2-1)*x[2]-x[1])

D = TensorField(ProductSpace(-6:0.01:6,-6:0.01:6))
streamplot(fun.(D,0.2))
streamplot(fun.(D,0.5))

D = TensorField(ProductSpace(-4:0.01:4,-4:0.01:4))
streamplot(fun.(D,1))

D = TensorField(ProductSpace(-2:0.01:2,-3:0.01:3))
streamplot(fun.(D,3))

