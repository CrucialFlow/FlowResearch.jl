
# my implementation
using Grassmann, Cartan, Adapode

# plotting software (optional)
using Makie, GLMakie

sphere(x) = Chain(cos(x[2])*sin(x[1]), sin(x[2])*sin(x[1]), cos(x[1]))

# assign sphere with 60 by 60 points
sph = sphere.(SphereParameter(60,60));
sphmet = surfacemetric(sph);
coef = secondkind(sphmet);

# initial conditions
x0 = Chain(1.0,1.0)
v0 = Chain(1.0,1.0)
sol = geosolve(coef,x0,v0,2π,15,1,4);

geosphere(x) = Chain(cos(x[2])*cos(x[1]), cos(x[2])*sin(x[1]), sin(x[2]))

# alternatively, compute coefficients explicitly:
sphcoef(x) = TensorOperator(Chain(Chain(Chain(-sin(x[1])*cos(x[1]),0),Chain(0,cot(x[1]))),Chain(Chain(0,cot(x[1])),Chain(0,0))))
coef = sphcoef.(SphereParameter(60,60))

circle2(x) = Chain(cos(3x[1]),sin(2x[1]))
f(x) = sin(x[1]/2)*sin(x[2])
vf(x) = Chain(cos(x[1])*cos(x[2]),sin(x[2])*sin(x[1]))

cf = circle2.(TorusParameter(100,100));
tf = f.(TorusParameter(100,100));
gtf = gradient(tf);
vft = vf.(TorusParameter(100,100));

streamplot(𝓛[cf,gtf,vft])

halfplane(x) = TensorOperator(Chain(Chain(Chain(0.0,inv(x[2])),Chain(-inv(x[2]),0.0)),Chain(Chain(-inv(x[2]),0.0),Chain(0.0,-inv(x[2])))))

