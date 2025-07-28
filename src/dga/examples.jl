
# parametric
θ = LinRange(0, pi, 30)
φ2 = LinRange(0, 2pi, 60)
fdom = GridFrameBundle(φ2⊕φ2);

sphere(x) = 0.5Chain(cos(x[2])*sin(x[1]), sin(x[2])*sin(x[1]), cos(x[1]))
torus(x) = Chain((2+0.5cos(x[1]))*cos(x[2]), (2+0.5cos(x[1]))*sin(x[2]), 0.5sin(x[1]))
mobius(x) = Chain((2+x[2]*cos(x[1]/2))*cos(x[1]), (2+x[2]*cos(x[1]/2))*sin(x[1]), x[2]*sin(x[1]/2))

SPH = sphere.(GridFrameBundle(θ⊕φ2))
TOR = torus.(fdom)

mesh(TensorField(SPH, rand(30, 60)), colormap=:blues)
mesh(TensorField(TOR, rand(60, 60)), colormap=:blues)

rad(x) = (1+0.3sin(3x[1])+0.1cos(7x[2]))
rad2(x) = (3+0.5cos(x[2]))
wiggle(x) = Chain((rad2(x)+rad(x)*cos(x[1]))*cos(x[2]), (rad2(x)+rad(x)*cos(x[1]))*sin(x[2]), rad(x)*sin(x[1]))
WIG = fdom → wiggle
arr = -⋆det(gradient(WIG));
kurv = Real(abs(det(gradient(arr/abs(arr)))));

mesh(WIG → kurv)

# Klein

klein(x) = klein(x[1],x[2]/2)
function klein(v,u)
    x = cos(u)*(-2/15)*(3cos(v)-30sin(u)+90sin(u)*cos(u)^4-60sin(u)*cos(u)^6+5cos(u)*cos(v)*sin(u))
    y = sin(u)*(-1/15)*(3cos(v)-3cos(v)*cos(u)^2-48cos(v)*cos(u)^4+48cos(v)*cos(u)^6-60sin(u)+5cos(u)*cos(v)*sin(u)-5cos(v)*sin(u)*cos(u)^3-80cos(v)*sin(u)*cos(u)^5+80cos(v)*sin(u)*cos(u)^7)
    z = sin(v)*(2/15)*(3+5cos(u)*sin(u))
    Chain(x,y,z)
end

# Hopf

proj(x) = Chain(x[2],x[3],x[4])/Real(abs(x))
stereo(x) = Chain(x[2],x[3],x[4])/(1-x[1])

hopf(x) = hopf(x[1],x[2],x[3])
function hopf(θ,φ,ψ)
    a = cos(θ)*exp((im/2)*(ψ-φ))
    b = sin(θ)*exp((im/2)*(ψ+φ))
    Chain(real(a),imag(a),real(b),imag(b))
end

hs = stereo.(hopf.(HopfParameter()));
hp = proj.(hopf.(HopfParameter()));

hopfstereo(x) = stereo(hopf(x))
hopfproj(x) = proj(hopf(x))

# φ: 0,2π ; ψ: 0,4π ; θ: 0,π
# θ: 0,π  ; ψ: 0,4π ; φ: 0,2π
# θ: 0,π  ; φ: 0,2π ; ψ: 0,4π

HopfParameter(x=7,y=9,z=60) = HopfParameter(Values(x,y,z))
HopfParameter(n::Values{3,Int}) = TensorField(GridFrameBundle(PointArray((LinRange(7π/4/n[1],7π/4,n[1])⊕LinRange(-π/2,π/2,n[2]))⊕LinRange(0,4π,n[3])),OpenTopology(n)))

HopfParameter(x=7,y=18,z=60) = HopfParameter(Values(x,y,z))
HopfParameter(n::Values{3,Int}) = TensorField(GridFrameBundle(PointArray((LinRange(7π/8/n[1],7π/8,n[1])⊕LinRange(0,2π,n[2]))⊕LinRange(0,4π,n[3])),OpenTopology(n)))

# bracket

circle2(x) = Chain(cos(3x[1]),sin(2x[1]))
f(x) = sin(x[1]/2)*sin(x[2])
vf(x) = Chain(cos(x[1])*cos(x[2]),sin(x[2])*sin(x[1]))

cf = circle2.(TorusParameter(100,100));
tf = f.(TorusParameter(100,100));
gtf = gradient(tf);
vft = vf.(TorusParameter(100,100));

𝓛[cf,gtf]
Lie[gtf,vft,cf]
𝓛[cf,gtf,vft]
streamplot(𝓛[cf,gtf,vft])

# Video

using Grassmann, Cartan, Adapode, FlowGeometry, MATLAB, Makie, GLMakie

# experiment with t
t = TensorField(0:0.01:2π)
exp(v12*t)+exp(v12*t/2)+exp(v12*2t)

# torus
φ = LinRange(0, 2π, 60)
fundom = TensorField(φ⊕φ);
wobble(x) = (1+0.3sin(3x[1])+0.1cos(7x[2]))
wumble(x) = (3+0.5cos(x[2]))
wiggle(x) = Chain((wumble(x)+wobble(x)*cos(x[1]))*cos(x[2]), (wumble(x)+wobble(x)*cos(x[1]))*sin(x[2]), wobble(x)*sin(x[1]))
obj = wiggle.(fundom)
arr = -⋆det(gradient(obj));
gauss = Real(abs(det(gradient(arr/abs(arr)))));

mesh(TensorField(fiber(obj),fiber(gauss)))

met = surfacemetric(obj)
@basis MetricTensor([1 1; 1 1])
circ = (2+cos(t))*v1+(2+sin(t))*v2
circ2 = TensorField(met(circ),fiber(circ))
circ3 = obj(circ)

# integration
square = (-2:0.003:2)⊕(-2:0.003:2)
disk = TensorField(square, x->abs(x)<1)
integrate(disk)

# differential equations
x0 = Chain(10.0,10.0,10.0)
Lorenz(σ,r,b) = x -> Chain(σ*(x[2]-x[1]), x[1]*(r-x[3])-x[2], x[1]*x[2]-b*x[3])
sol = odesolve(Lorenz(10.0,60.0,8/3),Chain(10.0,10.0,10.0),2π,15,3,4)

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

halfplane(x) = TensorOperator(Chain(Chain(Chain(0.0,inv(x[2])),Chain(-inv(x[2]),0.0)),Chain(Chain(-inv(x[2]),0.0),Chain(0.0,-inv(x[2])))))
