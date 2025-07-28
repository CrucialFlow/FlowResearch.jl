
# utility

function myplot(x,lw::Int=10;args...)
    fig = lines(x;colormap=:greys,linewidth=lw,args...)
    lines!(x;color=:black,linestyle=:dash,args...)
    return fig
end
function myplot!(x,lw::Int=10;args...)
    lines!(x;colormap=:greys,linewidth=lw,args...)
    lines!(x;color=:black,linestyle=:dash,args...)
end

# plane curves

t = TensorField(0:0.01:4*pi)
lin = Chain.(cos(t)*t,sin(t)*11+t)
fig = myplot(lin,rasterize=2)
scaledarrows!(lin,unitframe(lin),gridsize=50)
save("lin1.pdf",fig)
save("lin2.pdf",myplot(arclength(lin),rasterize=2))
save("lin3.pdf",myplot(speed(lin),rasterize=2))
save("lin4.pdf",myplot(curvature(lin),rasterize=2))
save("lin5.pdf",myplot(planecurve(t),rasterize=2))
save("lin6.pdf",myplot(planecurve(cos(t)*t),rasterize=2))
save("lin8.pdf",myplot(planecurve(cos(t)-t*sin(t)),rasterize=2))
t = TensorField(0:0.001:4*pi)
save("lin7.pdf",myplot(planecurve(cos(t*t)*t),rasterize=2))

# Lorenz

Lorenz2(s,r,b) = x -> Chain(
    s*(-x[1]-x[2]),
    x[2]*(r-x[3])+x[1],
    -x[2]*x[1]-b*x[3])
Lorenz(s,r,b) = x -> Chain(
    s*(x[2]-x[1]),
    x[1]*(r-x[3])-x[2],
    x[1]*x[2]-b*x[3])
p = TensorField(ProductSpace(-40:0.2:40,-40:0.2:40,10:0.2:90))
rescal2(x) = Chain(-x[2],x[1],x[3])
vf = Lorenz2(10.0,60.0,8/3).(p)
save("lor1.pdf",streamplot(vf,gridsize=(10,10),colormap=:grays,rasterize=3))
x0 = Chain(10.0,10.0,10.0)
sol = odesolve(Lorenz(10.0,60.0,8/3),x0,2*pi,15,3,4)
save("lor2.pdf",myplot(rescal2.(sol),2,rasterize=4))

# Riemann sphere

pts = TensorField(-2*pi:0.0001:2*pi)
@basis S"∞+++"
f(t) = (↓(exp(π*t*((3/7)*v12+v∞3))>>>↑(v1+v2+v3)))
save("explin1.pdf",myplot(V(2,3,4).(f.(pts)),2,rasterize=3))
@basis S"∞∅+++"
f(t) = (↓(exp(π*t*((3/7)*v12+v∞3))>>>↑(v1+v2+v3)))
rescal(x) = Chain(x[1],x[2],x[3]/10)
out = V(3,4,5).(vector.(f.(pts)))
save("explin2.pdf",myplot(rescal.(out),2,rasterize=3))

@basis S"∞+++"
f(t) = ↓(exp(t*v∞*(sin(3t)*3v1+cos(2t)*7v2-sin(5t)*4v3)/2)>>>↑(v1+v2-v3))
save("explin3.pdf",myplot(V(2,3,4).(f.(pts)),2,rasterize=3))


@basis S"∞+++"
f(t) = ↓(exp(t*(v12+0.07v∞*(sin(3t)*3v1+cos(2t)*7v2-sin(5t)*4v3)/2))>>>↑(v1+v2-v3))
save("explin4.pdf",myplot(V(2,3,4).(f.(pts)),2,rasterize=3))

# bivector

basis"2"
vdom = TensorField(ProductSpace{V}(-1.5:0.1:1.5,-1.5:0.1:1.5))
save("exp1.pdf",streamplot(tensorfield(exp(pi*v12/2)).(vdom),colormap=:grays,rasterize=2))
save("exp2.pdf",streamplot(tensorfield(exp((pi/2)*v12/2)).(vdom),colormap=:grays,rasterize=2))
save("exp3.pdf",streamplot(tensorfield(exp((pi/4)*v12/2)).(vdom),colormap=:grays,rasterize=2))
save("exp4.pdf",streamplot(tensorfield(v1*exp((pi/4)*v12/2)).(vdom),colormap=:grays,rasterize=2))
@basis S"+-"
vdom = TensorField(ProductSpace{V}(-1.5:0.1:1.5,-1.5:0.1:1.5))
save("exp5.pdf",streamplot(tensorfield(exp((pi/8)*v12/2)).(vdom),colormap=:grays,rasterize=2))
save("exp6.pdf",streamplot(tensorfield(v1*exp((pi/4)*v12/2)).(vdom),colormap=:grays,rasterize=2))

@basis S"∞+++"
vdom1 = TensorField(ProductSpace{V(1,2,3)}(-1.5:0.1:1.5,-1.5:0.1:1.5,-1.5:0.1:1.5))
vf1 = tensorfield(exp((pi/4)*(v12+v∞3)),V(2,3,4)).(vdom1)
fig = streamplot(vf1,gridsize=(10,10),colormap=:grays,rasterize=3)
ff,ax,pl = fig
rotate_cam!(ax.scene,π/16,0.9π/2,0)
save("exp7.pdf",fig)

vdom2 = TensorField(ProductSpace{V(2,3,4)}(-1.5:0.1:1.5,-1.5:0.1:1.5,-1.5:0.1:1.5))
vf2 = tensorfield(exp((pi/4)*(v12+v∞3)),V(2,3,4)).(vdom2)
fig = streamplot(vf2,gridsize=(10,10),colormap=:grays,rasterize=3)
ff,ax,pl = fig
rotate_cam!(ax.scene,1.2π/16,1.2π/16,0)
save("exp8.pdf",fig)

# Bracket

circle2(x) = Chain(cos(3x[1]),sin(2x[1]))
f(x) = sin(x[1]/2)*sin(x[2])
vf(x) = Chain(cos(x[1])*cos(x[2]),sin(x[2])*sin(x[1]))

cf = circle2.(TorusParameter(100,100));
tf = f.(TorusParameter(100,100));
gtf = gradient(tf);
vft = vf.(TorusParameter(100,100));

ltf1 = Lie[cf,gtf]
ltf2 = Lie[cf,gtf,vft]

lie1,lie2 = ltf1,ltf2

save("lie1.pdf",streamplot(lie1,colormap=:grays,rasterize=2))
save("lie2.pdf",streamplot(lie2,colormap=:grays,rasterize=2))

# integrals

t = TensorField(0:0.01:2pi)
circ = Chain.(cos(t),sin(t))
sphere(x) = Chain(
    cos(x[2])*sin(x[1]), sin(x[2])*sin(x[1]), cos(x[1]))
sph = sphere.(SphereParameter(15,30))

save("circ.pdf",lines(circ,color=:black,linestyle=:dash))
save("sphere.pdf",wireframe(sph,color=:black))

t = TensorField(0:0.01:2pi);
gamma(t) = Chain(cos(t[1]),sin(t[1]),0)
sigma(t) = Chain(0,1+cos(t[1]),sin(t[1]))
(linknumber(gamma.(t),sigma.(t)), 1.0)
fig = lines(gamma.(t),color=:black,linestyle=:dash)
lines!(sigma.(t),color=:black,linestyle=:dash)
save("link1.pdf",fig)

save("link2.pdf",mesh(linkmap(gamma.(t),sigma.(t)),normalnorm,colormap=:grays,rasterize=5))
save("link3.pdf",mesh(unit(linkmap(gamma.(t),sigma.(t))),normalnorm,colormap=:grays,rasterize=5))

# curvature

torus(x) = Chain(
    (2+0.5cos(x[1]))*cos(x[2]),
    (2+0.5cos(x[1]))*sin(x[2]),
    0.5sin(x[1]))
wobble(x) = (1+0.3sin(3x[1])+0.1cos(7x[2]))
wumble(x) = (3+0.5cos(x[2]))
wiggle(x) = Chain(
    (wumble(x)+wobble(x)*cos(x[1]))*cos(x[2]),
    (wumble(x)+wobble(x)*cos(x[1]))*sin(x[2]),
    wobble(x)*sin(x[1]))
tor = torus.(TorusParameter())
wig = wiggle.(TorusParameter())

save("tor1.pdf",mesh(tor,normalnorm,colormap=:grays,rasterize=6))
save("tor2.pdf",mesh(tor,gaussextrinsic,colormap=:grays,rasterize=6))
save("tor3.pdf",mesh(tor,gaussintrinsic,colormap=:grays,rasterize=6))
save("tor4.pdf",mesh(tor,gaussextrinsicnorm,colormap=:grays,rasterize=6))
save("tor5.pdf",mesh(tor,gausssign,colormap=:grays,rasterize=6))
#save("tor4.pdf",scaledarrows(tor,unitframe(tor),gridsize=(11,21)))
save("wig1.pdf",mesh(wig,normalnorm,colormap=:grays,rasterize=6))
save("wig2.pdf",mesh(wig,gaussextrinsic,colormap=:grays,rasterize=6))
save("wig3.pdf",mesh(wig,gaussintrinsic,colormap=:grays,rasterize=6))

# geodesic

tormet = surfacemetric(tor);
torcoef = secondkind(tormet);
solg = geosolve(torcoef,Chain(1.0,1.0),Chain(1.0,sqrt(2)),10π,7,1,4)
@basis MetricTensor([1 1; 1 1]) # abstract non-Euclidean V
solm = TensorField(tormet(solg),Chain{V}.(value.(fiber(solg))))
save("torgeo1.pdf",myplot(solg,3,rasterize=1))
save("torgeo2.pdf",myplot(torus.(solg),3,rasterize=3))
fig = myplot(arclength(solg),3,rasterize=2)
myplot!(arclength(solm),3,rasterize=2)
save("torgeo3.pdf",fig)

# klein

klein(x) = klein(x[1],x[2]/2)
function klein(v,u)
    x = cos(u)*(-2/15)*(3cos(v)-30sin(u)+90sin(u)*cos(u)^4-
        60sin(u)*cos(u)^6+5cos(u)*cos(v)*sin(u))
    y = sin(u)*(-1/15)*(3cos(v)-3cos(v)*cos(u)^2-
        48cos(v)*cos(u)^4+48cos(v)*cos(u)^6-
        60sin(u)+5cos(u)*cos(v)*sin(u)-
        5cos(v)*sin(u)*cos(u)^3-80cos(v)*sin(u)*cos(u)^5+
        80cos(v)*sin(u)*cos(u)^7)
    z = sin(v)*(2/15)*(3+5cos(u)*sin(u))
    Chain(x,y,z)
end

rescal3(x) = Chain(-x[1],-x[2],x[3])
kle = klein.(KleinParameter(100,100))
klemet = surfacemetric(kle);
klecoef = secondkind(klemet);
solg = geosolve(klecoef,Chain(1.0,1.0),Chain(1.0,2.0),2pi,7,1,4);
save("klegeo1.pdf",myplot(solg,3,rasterize=5))
save("klegeo2.pdf",myplot(rescal3.(klein.(solg)),3,rasterize=5))
save("klegeo22.pdf",myplot(klein.(solg),3,rasterize=5))

kle = klein.(KleinParameter(60,60))
klemet = surfacemetric(kle);
klecoef = secondkind(klemet);
save("kle.pdf",wireframe(rescal3.(kle),color=:black,linewidth=0.3))

# hyperbolic

halfplane(x) = TensorOperator(Chain(
    Chain(Chain(0.0,inv(x[2])),Chain(-inv(x[2]),0.0)),
    Chain(Chain(-inv(x[2]),0.0),Chain(0.0,-inv(x[2])))))
z1 = geosolve(halfplane,Chain(1.0,1.0),Chain(1.0,2.0),10pi,7)
z2 = geosolve(halfplane,Chain(1.0,0.1),Chain(1.0,2.0),10pi,7)
z3 = geosolve(halfplane,Chain(1.0,0.5),Chain(1.0,2.0),10pi,7)
z4 = geosolve(halfplane,Chain(1.0,1.0),Chain(1.0,1.0),10pi,7)
z5 = geosolve(halfplane,Chain(1.0,1.0),Chain(1.0,1.5),10pi,7)
fig = myplot(z1,rasterize=2); myplot!(z2,rasterize=2); myplot!(z3,rasterize=2); myplot!(z4,rasterize=2); myplot!(z5,rasterize=2)
save("halfplane.pdf",fig)

# hopf

stereohopf(x) = stereohopf(x[1],x[2],x[3])
function stereohopf(theta,phi,psi)
    a = cos(theta)*exp((im/2)*(psi-phi))
    b = sin(theta)*exp((im/2)*(psi+phi))
    Chain(imag(a),real(b),imag(b))/(1-real(a))
end
hs = stereohopf.(HopfParameter());
fig = alteration!(hs,wireframe,wireframe!)
ff,ax,pl = fig
zoom!(ax.scene,0.1)
save("hopf.pdf",fig)

# tangent spaces

save("stream1.pdf",streamplot(sph,vf2,colormap=:grays,rasterize=3))
save("stream2.pdf",streamplot(tor,vf3,colormap=:grays,rasterize=3))

# odesolve

#start0(x) = Chain(cos(x),sin(x),sin(3x))
start0(x) = Chain(cos(x),sin(x),sin(2x))
start0(x) = Chain(cos(x),sin(x),cos(1.5x)*sin(1.5x)/5)
x1 = start0.(TorusParameter(180));
darios(t,dt=tangent(fiber(t))) = hodge(wedge(dt,tangent(dt)))
sol1 = odesolve(darios,x1,1.0,11)
save("odesolve-example.pdf",mesh(sol1,normalnorm,colormap=:grays,rasterize=4))

# krylov

using MATLAB
pt,pe = initmesh(circleg,"hmax"=>0.1)
A,M = assemble(pt,1,1,0)
using KrylovKit # provides general eigsolve
yi,xi = geneigsolve((A,M),10,:SR;krylovdim=100) # maxiter=100
amp = TensorField.(Ref(pt),xi)
mode = TensorField.(graphbundle.(amp),xi)

fig = mesh(mode[9],colormap=:grays,rasterize=5)
wireframe!(pt,color=:black)
save("krylov9.pdf",fig)

# airfoil

pt,pe = initmesh(decsg(NACA"6511"),"hmax"=>0.1)
tf = solvepoisson(pt,pe,1,0,
    x->(x[2]>3.49 ? 1e6 : 0.0),0,x->(x[2]<-1.49 ? 1.0 : 0.0))
function kappa(z); x = base(z)
    if x[2]<-1.49 || sqrt((x[2]-0.5)^2+x[3]^2)<0.51
        1e6
    elseif x[2]>3.49
        fiber(z)[1]
    else
        0.0
    end
end
gtf = -gradient(tf)
kf = kappa.(gtf(immersion(pe)))
tf2 = solvetransportdiffusion(gtf,kf,0.01,1/50,
    x->(sqrt((x[2]-0.5)^2+x[3]^2)<0.7 ? 1.0 : 0.0))
save("air1.pdf",wireframe(pt,color=:black))
save("air2.pdf",streamplot(gtf,-0.3..1.3,-0.2..0.2,colormap=:grays,rasterize=2))
save("air3.pdf",mesh(-tf2,colormap=:grays,rasterize=3))

# sphere

using Grassmann, Cartan, FlowGeometry, MiniQhull, TetGen,Makie
pt,pe = FlowGeometry.spheremesh()
tf = solvepoisson(pt,pe,1,0,
    x->(x[2]>1.99 ? 1e6 : 0.0),0,x->(x[2]<-1.99 ? 1.0 : 0.0))
ps = sphere(sphere(∂(delaunay(PointCloud(sphere())))))
fig = wireframe(ps,color=:black)
streamplot!(gradient(tf),-1.1..1.1,-1.1..1.1,-1.1..1.1,gridsize=(10,10,10),colormap=:grays,rasterize=3)
save("sphereflow.pdf",fig)

# stokes

square = TensorField(ProductSpace(-3:0.03:3,-3:0.03:3))
cube = TensorField(ProductSpace(-4:0.1:4,-4:0.1:4,-1:0.2:10))
disk = (x->float(fiber(abs(x))<3)).(square)
paraboloid(x) = 9-x[1]*x[1]-x[2]*x[2]
para = graph(disk*paraboloid.(square))
F(x) = Chain(2x[3]-x[2],x[1]+x[3],3x[1]-2x[2])
vf = F.(cube)

fig = mesh(para,normalnorm,colormap=:grays,rasterize=4)
scaledarrows!(para,disk*unitnormal(para),gridsize=(22,22),rasterize=2)
streamplot!(vf,gridsize=(11,11,11),colormap=:grays,rasterize=3)
save("stokes1.pdf",fig)

