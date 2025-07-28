# github.com/chakravala
# DGA: Differential Geometry of Curves and Surfaces, Kobayashi


# Example 2.1.1 Sphere

tt = TensorField(LinRange(-pi/2,pi/2,60))
rev = revolvesphere(unitcircle(tt))

# Example 2.1.2 Ellipsoid

ellipsoidmap(r::Real=1) = ellipsoidmap(r,r,r)
ellipsoidmap(x) = ellipsoidmap(x[1],x[2],x[3])
ellipsoidmap(a,b,c) = x->Chain(a*cos(x[2])*cos(x[1]),b*cos(x[2])*sin(x[1]),c*sin(x[2]))

# Example 2.1.3 Hyperboloid of one sheet

hyperbolamap(x,n::Int=60) = hyperbolamap(x,TensorField(LinRange(-π/2,π/2,60)))
hyperbolamap(x,t) = Chain.(x[1].*cosh.(t),x[2].*sinh.(t))

revolve(hyperbolamap((3,4)))

# Example 2.1.4 Hyperboloid of two sheets

revolve(!hyperbolamap((3,4)))

# Example 2.1.5 Elliptic paraboloid

ellipticparaboloidfun(a,b) = x->(x[1]/a)^2+(x[2]/b)^2
ellipticparaboloidfun(1,1).(OpenParameter())

# Example 2.1.6 Hyperbolic paraboloid

hyperbolicparaboloidfun(a,b) = x->-(x[1]/a)^2+(x[2]/b)^2
hyperbolicparaboloidfun(1,1).(OpenParameter())

# Example 2.1.7 Torus

t = TorusParameter(60)
circ = Chain.(0.1cos(t)+0.5,0.1sin(t))
revolve(circ)

# Problem 2.1.1 (i)

ellipsemap(x,n::Int=60) = ellipsemap(x,TorusParameter(n))
ellipsemap(x,t) = Chain.(x[1].*cos.(t),x[2].*sin.(t))

ellipsoid(a::Real,b,c,n::Int=60) = ellipsoid(a,b,c,TensorField(LinRange(-π/2,π/2,60)),n)
ellipsoid(a::Real,b,c,tt,n::Int=60) = ellipsoid(ellipsemap((1,c),tt),a,b,n)
ellipsoid(f,a,b,n::Int=60) = revolvesphere(f,ellipsemap((a,b),n))
ellipsoid(2,3,4)

ellipticparaboloid(a::Real,b,c,n::Int=60) = ellipticparaboloid(a,b,c,TensorField(LinRange(0,π/2,60)),n)
ellipticparaboloid(a::Real,b,c,tt,n::Int=60) = ellipticparaboloid(c*tt^2,a,b,n)
ellipticparaboloid(f,a,b,n::Int=60) = revolvepolar(f,ellipsemap((a,b),n))
ellipticparaboloid(2,3,4)

# Problem 2.1.1 (ii)

hyperbolicparaboloid(a::Real,b,c,n::Int=60) = hyperbolicparaboloid(a,b,c,TensorField(LinRange(0,π/2,60)),n)
hyperbolicparaboloid(a::Real,b,c,tt,n::Int=60) = hyperbolicparaboloid(c*tt^2,a,b,n)
hyperbolicparaboloid(f,a,b,n::Int=60) = revolvepolar(f,hyperbolamap((a,b),n))
hyperbolicparaboloid(1,1,1)

# Problem 2.1.1 (iii)

fun(a,b,c) = x->Chain(a*(x[1]+x[2]),b*(x[1]-x[2]),c*x[1]*x[2])
fun(1,1,4).(OpenParameter())

# Problem 2.1.1 (iv)

fun(a,b,c) = x->Chain(a*(x[1]-x[2]),b*(x[1]*x[2]+1),c*(x[1]*x[2]-1))/(x[1]+x[2])
fun(1,1,1).(OpenParameter()*10)

# Example 2.3.7 curvature K = -c^2

t = OpenParameter(60)
c = 1
revolve(Chain.(exp(-c*t)/c,integral(sqrt(1-exp((-2c)*t)))))

# Example 2.3.8

t = TensorField(-1:0.01:1)
a = 1
revolve(Chain.(sqrt(t^2+a^2),a*asinh(t/a)))

# Problem 2.3.3

enneper(x) = Chain(3x[1]+3x[1]*x[2]^2-x[1]^3,x[2]^3-3x[2]-3x[1]^1*x[2],3*(x[1]^2-x[2]^2))


