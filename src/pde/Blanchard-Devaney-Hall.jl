# github.com/chakravala
# Math: Differential Equations, Blanchard-Devaney-Hall

# 5 Nonlinear Systems

# Van der Pol

D = TensorField(ProductSpace(-3:0.1:3,-3:0.1:3))
fun(x) = Chain(x[2],-x[1]+(1-x[1]^2)*x[2])
streamplot(fun.(D))
jacobian(fun.(D))(0,0)
eigvals(jacobian(fun.(D))(0,0))

# Competing Species

D = TensorField(ProductSpace(-1:0.1:4,-1:0.1:4))
fun(x) = Chain(2x[1]*(1-x[1]/2)-x[1]*x[2],3x[2]*(1-x[2]/3)-2x[1]*x[2])
streamplot(fun.(D))
jacobian(fun.(D))(0,0)
jacobian(fun.(D))(1,1)
jacobian(fun.(D))(2,0)
jacobian(fun.(D))(0,3)

D = TensorField(ProductSpace(-2:0.1:3,-2:0.1:3))
fun(x) = Chain(-x[1]-x[2]-x[1]^2-x[1]*x[2],-2x[1]-x[2]-2x[1]*x[2]-x[2]^2)
streamplot(fun.(D))
jacobian(fun.(D))(0,0)

# Nonpolynomial

D = TensorField(ProductSpace(-2π:0.1:2π,-2π:0.1:2π))
fun(x) = Chain(x[2],-x[2]-sin(x[1]))
streamplot(fun.(D))
jacobian(fun.(D))(0,0)
eigvals(jacobain(fun.(D))(0,0))

# More Examples

D = TensorField(ProductSpace(-1:0.1:2,-1:0.1:2))
fun(x) = Chain(-2x[1]+2x[1]^2,-3x[1]+x[2]+3x[1]^2)
heatmap(det(gradient(fun.(D))))
streamplot!(fun.(D))

jacobian(fun.(D))(0,0)
jacobian(fun.(D))(1,0)

# When Linearization Fails

D = TensorField(ProductSpace(-0.7:0.01:0.7,-0.7:0.01:0.7))
fun(x) = Chain(x[2]-(x[1]^2+x[2]^2)*x[1],-x[1]-(x[1]^2+x[2]^2)*x[2])
streamplot(fun.(D))
jacobian(fun.(D))(0,0)

fun(x) = Chain(x[2]+(x[1]^2+x[2]^2)*x[1],-x[1]+(x[1]^2+x[2]^2)*x[2])
streamplot(fun.(D))
jacobian(fun.(D))(0,0)

# 5.2 Qualitative Analysis

D = TensorField(ProductSpace(-1:0.1:3,-1:0.1:4))
fun(x) = Chain(2x[1]*(1-x[1]/2)-x[1]*x[2],3x[2]*(1-x[2]/3)-2x[1]*x[2])
streamplot(fun.(D))

heatmap(getindex.(fun.(D),1))
contour!(getindex.(fun.(D),1),levels=[0],color=:black)

heatmap(getindex.(fun.(D),2))
contour!(getindex.(fun.(D),2),levels=[0],color=:black)

contour(getindex.(fun.(D),1),levels=[0],color=:black)
contour!(getindex.(fun.(D),2),levels=[0],color=:black)

# Nullcines That Are Not Lines

D = TensorField(ProductSpace(-1:0.01:3,-1:0.01:3))
fun(x) = Chain(2x[1]*(1-x[1]/2)-x[1]*x[2],x[2]*(9/4-x[2]^2)-x[1]^2*x[2])
contour(getindex.(fun.(D),1),levels=[0],color=:black)
contour!(getindex.(fun.(D),2),levels=[0],color=:black)

fun(x) = Chain(2x[1]*(1-x[1]/2),x[2]*(9/4-x[2]^2))

# Using All Our Tools

D = TensorField(ProductSpace(-2:0.01:2,-2:0.01:2))
fun(x) = Chain(x[1]+x[2]-x[1]^3,x[1]/-2)
contour(getindex.(fun.(D),1),levels=[0],color=:black)
contour!(getindex.(fun.(D),2),levels=[0],color=:black)

streamplot(fun.(D))
jacobian(fun.(D))(0,0)

# 5.3 Hamiltonian Systems




