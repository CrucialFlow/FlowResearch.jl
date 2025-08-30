# github.com/chakravala

# Chapter 7

# 7.1 Sinks, Saddles, and Sources

V = Cartan.affmanifold(2)
B = Endomorphism{V}([-1 -5; 5 -1])
D = TensorField(ProductSpace(-5:0.1:5,-5:0.1:5))
streamplot(B⋅D)
streamplot(-B⋅D)

C = Endomorphism{V}([2 1; -1 -3])
streamplot(C⋅D)

# 7.2 Phase Portraits of Sinks

C = Endomorphism{V}([-1 2; -5 0])
streamplot(C⋅D)

C = Endomorphism{V}([-1 1; 0 -2])
streamplot(C⋅D)

C = Endomorphism{V}([-1/2 0; 0 -1/2])
streamplot(C⋅D)

# 7.3 Phase Portraits of Nonhyperbolic Systems

eigmat(z) = Endomorphism{Cartan.affmanifold(2)}([real(z) -imag(z); imag(z) real(z)])
streamplot(eigmat(0+3im)⋅D)
streamplot(eigmat(1+3im)⋅D)
streamplot(eigmat(-1+3im)⋅D)

C = Endomorphism{V}([-2 4; 1 -2])
streamplot(C⋅D)

system(μ) = Endomorphism{Cartan.affmanifold(2)}([μ 0; 0 -1])
streamplot(system(-0.2)⋅D)
streamplot(system(0)⋅D)
streamplot(system(0.2)⋅D)

C = Endomorphism{V}([0 1; 0 0])
streamplot(C⋅D)

# Chapter 11

# 11.1 Introduction

D = TensorField(ProductSpace(-5.5:0.1:5.5,-0.6:0.1:1.6))
fun(x) = Chain(1.0,x[2]^3-3x[2]^2+2x[2])
streamplot(fun.(D))

V = Cartan.affmanifold(2)
C = Endomorphism{V}([0 1; 2.2 2.1])
eigvecs(C)
eigvals(C)

D = TensorField(ProductSpace(-5.5:0.1:5.5,-5.5:0.1:5.5))
streamplot(C⋅D)

fun(x) = Chain(x[2],2.2x[1]+2.1x[2]+x[1]^2+x[1]*x[2])
streamplot(fun.(D))

# Exercise 11.1.1

D = TensorField(ProductSpace(-3:0.1:1,-4:0.1:3))
fun(x) = Chain(1.0,2+3x[2]-4x[2]^2+x[2]^3-2x[2]^4-x[2]^5)
streamplot(fun.(D))

# Exercise 11.1.2

D = TensorField(ProductSpace(-5.5:0.1:5.5,-5.5:0.1:5.5))
fun(x) = Chain(x[2],2.2x[1]+2.1x[2]-x[1]^3+1.6x[1]*x[2])
streamplot(fun.(D))

# 11.2 Equilibria and Linearization

D = TensorField(ProductSpace(-5.5:0.1:5.5,-5.5:0.1:5.5))
fun(x) = Chain(x[2],2.2x[1]+2.1x[2]+x[1]^2+x[1]*x[2])
streamplot(fun.(D))
jacobian(fun.(D))(-2.2,0)

D = TensorField(ProductSpace(-2.5:0.1:10.5,-4.5:0.1:4.5))
fun(x) = Chain(1+x[1]-x[2]^2,-1+6x[2]+x[1]^2-5x[2]^2)
streamplot(fun.(D))
jacobian(fun.(D))(-1,0)
jacobian(fun.(D))(0,1)
jacobian(fun.(D))(3,2)
jacobian(fun.(D))(8,-3)

D = TensorField(ProductSpace(-5.5:0.1:5.5,-5.5:0.1:5.5))
fun(x) = Chain(x[1]+x[2]+x[1]^2-x[2]^2+0.1,x[2]-2x[1]*x[2]+0.5x[1]^2+x[2]^2)
streamplot(fun.(D))

# Exercise 11.2.14

D = TensorField(ProductSpace(-5.5:0.1:5.5,-5.5:0.1:5.5))
fun(x) = Chain(2x[1]-x[2],11x[1]+x[1]^2-3x[2]^2)
streamplot(fun.(D))

# Exercise 11.2.15

D = TensorField(ProductSpace(-2.5:0.1:2.5,-2.5:0.1:2.5))
fun(x) = Chain(-x[1]-x[2]+x[1]^2-x[2]^2,0.25x[1]-3x[2]-2x[1]*x[2]+x[2]^2)
streamplot(fun.(D))

# Exercise 11.2.16

D = TensorField(ProductSpace(-2.5:0.1:2.5,-2.5:0.1:2.5))
fun(x) = Chain(3x[1]-2x[2]+3x[1]^2+x[2]^2,-x[1]+x[2]-3x[1]*x[2])
streamplot(fun.(D))

# 11.3

V = Cartan.affmanifold(2)
C = Endomorphism{V}([0 -3; 3 0])
D = TensorField(ProductSpace(-0.5:0.1:0.5,-0.5:0.1:0.5))
streamplot(C⋅D)

D = TensorField(ProductSpace(-2:0.1:2,-2:0.1:2))
fun(x) = Chain(-2x[2]-(x[1]^2+x[2]^2)*x[1],2x[1]-(x[1]^2+x[2]^2)*x[2])
streamplot(fun.(D))
fun(x) = Chain(x[1]-2x[2]-(x[1]^2+x[2]^2)*x[1],2x[1]+x[2]-(x[1]^2+x[2]^2)*x[2])
streamplot(fun.(D))

# Exercise 11.3.6

D = TensorField(ProductSpace(-2:0.1:2,-2:0.1:2))
fun(x) = Chain(x[1]-2x[2]-(x[1]^2+x[2]^2)*x[1]+0.1x[1]^2,2x[1]+x[2]-(x[1]^2+x[2]^2)*x[2]-0.2(x[2]^2-x[1]^2))
streamplot(fun.(D))

# Exercise 11.3.7

D = TensorField(ProductSpace(-2.5:0.1:2.5,-2.5:0.1:2.5))
fun(x) = Chain(x[1]-2x[2]+(x[1]^2+x[2]^2)^2*x[1]-4x[1]^3,2x[1]+x[2]+(x[1]^2+x[2]^2)^2*x[2]-7(x[2]^2-x[1]^2)*x[2])
streamplot(fun.(D))

# 11.4

D = TensorField(ProductSpace(-2:0.1:4,-4:0.1:2))
fun(x) = Chain(2x[1]-x[2]+3(x[1]^2-x[2]^2)+2x[1]*x[2],x[1]-3x[2]-3(x[1]^2-x[2]^2)+3x[1]*x[2])
streamplot(fun.(D))

# Exercise 11.4.1

D = TensorField(ProductSpace(-10:0.1:10,-10:0.1:10))
fun(x) = Chain(2x[1]-x[2]+3cos(x[1]+x[2]),x[1]-3x[2]-2x[1]*x[2])
streamplot(fun.(D))

# Exercise 11.4.2

D = TensorField(ProductSpace(-2:0.1:4,-3:0.1:6))
fun(x) = Chain(2x[1]-x[2]+3cos(x[1]+x[2]),x[1]-x[2]-2sin(x[1]-x[2]))
streamplot(fun.(D))

# Exercise 11.4.3

D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3))
fun(x) = Chain(2-x[2]^2-((x[1]^2-2)^2+(x[2]^2-2)^2-1)*(x[1]^2-2)-0.2x[1],x[1]^2-2-((x[1]^2-2)^2+(x[2]^2-2)^2-1)*(x[2]^2-2))
streamplot(fun.(D))

# Exercise 11.4.5

D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3))
fun(x) = Chain(0.04x[1]-x[2]-3x[2]^2+2.5x[1]*x[2],x[1]-3x[1]^2+2x[1]^2*x[2])
streamplot(fun.(D))
fun(x) = Chain(0.05x[1]-x[2]-3x[2]^2+2x[1]*x[2],x[1]-3x[1]^2+x[1]^2*x[2])
streamplot(fun.(D))



