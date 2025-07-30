# github.com/chakravala

# Chapter 1
pt,pe = initmesh(0:1/5:1)
M,b = assemblemassload(pt,x->x[2]*sin(x[2]))
M,b = assemblemassload(pt,x->2x[2]*sin(2π*x[2])+3)
tf = TensorField(pt,M\b)

# Chapter 2
pt,pe = initmesh(2:0.1:8)
tf = solvepoisson(pt,pe,x->0.1*(5-0.6x[2]),x->0.03*(x[2]-6)^4,[1e6,0],[-1,0])

g,pt,pe = refinemesh(0:0.2:1)

g,pt,pe = refinemesh(0:0.25:1)
adaptpoisson(g,pt,pe,1,0,x->exp(-100abs2(x[2]-0.5)),1e6)
tf = solvepoisson(pt,pe,1,x->exp(-100abs2(x[2]-0.5)),1e6)

#Chatper 3
pt,pe = initmesh(Rectg(0,0,2π,2π),"hmax"=>0.1)
M,b = assemblemassload(pt,x->x[2]*x[3])
tf = TensorField(pt,x->x[2]*x[3])

# Chapter 4
pt,pe = initmesh(Airfoilg,"hmax"=>0.5)
tf = solvepoisson(pt,pe,1,0,x->(x[2]>29.99 ? 1e6 : 0.0),0,x->(x[2]<-29.99 ? 1.0 : 0.0))

# try NACA"3412"

using KrylovKit
pt,pe = initmesh(circleg,"hmax"=>0.1)
A,M = assemble(pt,1,1,0)
La,Xi = geneigsolve((A,M),5,:SR;krylovdim=60) # maxiter=100
Ξ = TensorField.(Ref(pt),Xi)
Ξ2 = TensorField.(graphbundle.(Ξ),fiber.(Ξ))

g,pt,pe = refinemesh(Rectg(0,0,1,1),"hmax"=>0.25)
tf = solvepoisson(pt,pe,1,f,0)

g,pt,pe = refinemesh(Rectg(0,0,1,1),"hmax"=>0.1)
rad2(x) = (x[2]-0.5)^2+(x[3]-0.5)^2
r2 = rad2.(pt);
rr2 = rad2.(FaceBundle(pt));
f = 4*400^2*(1-400r2)*exp(-400r2);
ff = 4*400^2*(1-400rr2)*exp(-400rr2);
tf = Adapode.poisson(pt,pe,1,0,f,0);
η = jumps(pt,1,0,ff,tf);
refinemesh!(g,pt,pe,select(η))

# Chapter 5

pt,pe = initmesh(0:0.01:1)
triangle(x) = 0.5-abs(0.5-x[2])
bt = solveheat(triangle.(pt),(x->2x[2]).(pt),(x->1e6).(pe),range(0,0.6,101))

pt,pe = initmesh(Dslitg,"hmax"=>0.02)
cone(x) = 0.5-sqrt((0.5-x[2])^2+(0.5-x[3])^2)
be = solveheat(cone.(pt),(x->0).(pt),(x->1e6).(pe),range(0,0.6,101))

pt,pe = initmesh(Dslitg,"hmax"=>0.02)
DoubleSlitBC(pe,Δt=0.005,T=2) = TimeParameter(pe,x->x[2]<-0.24999,0:Δt:T)
cn = solvewave(pt,0.1sin(8π*DoubleSlitBC(pe)))

pt,pe = initmesh(decsg(NACA"2414"),"hmax"=>0.1)
function WingWaveBC(pe,Δt=0.01,T=2)
    TimeParameter(pe,x->(atan(x[2]-0.7)<0.3)&&(abs(x[3])<0.3),0:Δt:T)
end
cn = solvewave(pt,0.1sin(4π*WingWaveBC(pe)))
variation(cn,0.01,mesh,mesh!)

# Chapter 6

function Arnoldi(A,q,m)
    n = size(A)[1]
    Q = zeros(n,m+1)
    H = zeros(m+1,m)
    Q[:,1] = q/norm(q)
    for i ∈ 1:m
        z = A*Q[:,j]
        for i ∈ 1:j
            Qi = Q[:,i]
            Hij = z⋅Qi
            H[i,j] = Hij
            z -= Hij*Qi
        end
        Hj1j = norm(z)
        H[j+1,j] = Hj1j
        iszero(Hj1j) && break
        Q[:,j+1] = z/Hj1j
    end
    return Q,H
end

function Richardson(A,x,b,maxit,α)
    for k ∈ 1:maxit
        x += α*(b-A*x)
    end
    return x
end
'
function TwoGrid(n=25)
    nf = 2*n-1 # number of fine nodes
    nc = Int((nf-1)/2) # coarse nodes
    h = 1/(nf+1)
    x = 0:h:1 # mesh
    A1 = -ones(nf-1)
    A = spdiagm(-1=>A1,0=>2ones(nf),-1=>A1) # fine stiffness matrix
    b = ones(nf)*h # load vector
    u = zeros(nf) # solution guess
    P = spzeros(nf,nc) # prolongation matrix
    for i ∈ 1:nc
        P[2i-1,i] = 0.5
        P[2i,i] = 1
        P[2i+1,i] = 0.5
    end
    R = P'/2 # prolongation matrix
    RAP = R*A*P # coarse stiffness matrix
    for k ∈ Values(1,2,3,4,5)
        u = Richardson(A,u,b,4,0.25h)
        r = R*(b-A*u) # residual
        e = RAP\r # correction
        u += P*e # solution update
    end
end

# Chapter 8

pt,pe = initmesh(squareg)
pt2 = LagrangeP2(pt)
A = assemblestiffnessP2(pt)

# interpCR
pt,pe = initmesh(Rectg(0,0,1,1),"hmax"=>0.25)
crfun(x) = 1+2sin(3x[2])
ed = edges(pt)
ei = edgesindices(pt,FaceBundle(ed))
iCR = Cartan.interpCR(pt,ed,crfun.(ei))
surface(iCR/5)
wireframe!(pt)

# Chapter 9

Afcn(u) = 0.125 + u^2
Ffcn(x) = 1
pt,pe = initmesh(Rectg(0,0,1,1),"hmax"=>0.05)
solvenonlinearpoisson(TensorField(FaceBundle(pt),1),pe,Afcn)

pt,pe = initmesh(Rectg(0,0,1,1),"hmax"=>0.025)
bistable_ic(x) = cos(2π*x[2]^2)*cos(2π*x[3]^2)
ξ = solvebistable(bistable_ic.(pt))

pt,pe = initmesh(Rectg(0,0,1,1),"hmax"=>0.02)
solvebistable_newton(bistable_ic.(pt))

# Chapter 10

pt,pe = initmesh(squareg,"hmax"=>0.05)
tf = solvetransport(pt,pe,iterable(immersion(pt),x->Chain{↓(Cartan.varmanifold(3))}(1,1)))

pt,pe = initmesh(decsg(NACA"2414"),"hmax"=>0.1)
tf = solvepoisson(pt,pe,1,0,x->(x[2]>3.49 ? 1e6 : 0.0),0,x->(x[2]<-1.49 ? 1.0 : 0.0))
gtf = -gradient(tf)
function κ(z)
    x = base(z)
    if x[2]<-1.49 || sqrt(x[2]^2+x[3]^2)<0.51
        1e6
    elseif x[2]>3.49
        fiber(z)[1]
    else
        0.0
    end
end
tf2 = solvetransportdiffusion(gtf,κ.(gtf(immersion(pe))),0.01,0.1,x->(sqrt((x[2]-0.5)^2+x[3]^2)<0.7 ? 1.0 : 0.0))

# combine Poisson + transport diffusion, GLS Galerkin Least Squares stabilized

linesegments(p(edges(pt)))

function κ(z)
    x = base(z)
    if x[2]<-1.49 || sqrt((x[2]-0.5)^2+x[3]^2)<0.51
        1e6
    elseif x[2]>3.49
        fiber(z)[1]
    else
        0.0
    end
end
tf2 = solvetransportdiffusion(gtf,κ.(gtf(e)),0.01,1/150,x->(sqrt((x[2]-0.5)^2+x[3]^2)<0.7 ? 1.0 : 0.0))

streamplot(-gradient(laplacian(tf)),-0.5..3.0,-0.3..0.3)
lines!(NACA"2414")
streamplot(gtf,-0.3..1.3,-0.2..0.2)

# Chapter 11

pt,pe = initmesh(Rectg(0,0,1,1),"hmax"=>0.1)

force(x) = force(x[2],x[3])
force(x,y) = Chain((35/13)*(y-y^2)+(10/13)*(x-x^2),(-25/26)*(2y-1)*(2x-1))

UV = solveelastic(force.(pt),TensorField(pe,0))

elementresiduals(UV)+edgeresiduals(UV)

# eigen modes
K,M = assembleelastic(f,μ,λ)
w,D = geneigsolve(K,M,10,:SR)

# Chapter 12

pt,pe = totalmesh(squareg)
lid_bc(x) = Chain(x[3]>0.99 ? 1 : 0,0)
bc = lid_bc.(FaceBundle(pe))
ncv,ncp = solvestokes(pt,bc)

pt,pe = initmesh(DFGg,"hmax"=>0.25)
function channel_bc(pe) # inflow and outflow boundary conditions
    V = Cartan.varmanifold(3)(2,3)
    out = TensorField(pe(findall(x->x[2]>2.199999,fullpoints(pe))),1e6) # nodes on outflow
    inp = pe(findall(x->x[2]<0.000001,fullpoints(pe)))
    y = getindex.(points(inp),3); Umax = 0.3 # maximum inflow velocity
    ins = TensorField(inp,Chain{V}.((4Umax/0.41^2)*y.*(0.41.-y),0))
    return ins,out
end
ns = solvenavierstokes(pt,pe,channel_bc(pe)...,0.001,range(0,1,101))

# Chapter 13

pt,pe,sd = totalmeshes(Scatterg,"hmax"=>0.1)
scatter_ic_bc(sd,pe,ω=2π/1,ϵ=1) = (scatter_ic(sd,ω,ϵ),scatter_bc(pe,ω))
function scatter_ic(sd,ω=2π/1,ϵ=1)
    nt = length(sd)
    κ = Vector{Complex{Float64}}(undef,nt)
    p = points(sd)
    for i ∈ 1:nt
        _,x,y = value(p[i])
        sdi = fiber(sd)[i]
        σ = if sdi == 7 # cavity
            0
        elseif sd == 2 || sd == 9 # up down pml
            (abs(y)-3)^4/2^4
        elseif sd == 6 || sd == 8 # left right pml
            (abs(x)-3)^4/2^4
        else
            ((abs(x)-3)^4+(abs(y)-3)^4)/2^4
        end
        κ[i] = sqrt(Complex(-ϵ*(ω*ω),σ*ω))
    end
    TensorField(base(sd),κ)
end
function scatter_bc(in,ω=2π/1)
    e = fullimmersion(in)
    pe = FaceBundle(fullcoordinates(in),e)
    ne = length(e)
    fixed = Int[] # fixed nodes
    gvals = Complex{Float64}[] # nodal values of g
    for i ∈ 1:ne # loop over edges
        _,x,y = value(base(pe[i])) # coordinate of edge mid-point
        if sqrt(x*x+y*y)<1.001 # cylinder
            normal = -Chain(x,y)/sqrt(x*x+y*y)
            push!(fixed,i) # fix velocity nodes
            push!(gvals,normal[2]*exp(ω*y*im)) # bc values (u,v)
        elseif abs(x)>4.999 || abs(y)>4.999 # pml
            push!(fixed,i)
            push!(gvals,0)
        end
    end
    TensorField(pe(e[fixed]),gvals)
end
κ = scatter_ic(sd)
ef = Cartan.edgefacets(FaceBundle(pt))
bc = scatter_bc(ef)
ξ = solvemaxwell(κ,bc)
ξ = solvemaxwell(scatter_ic_bc(sd,pe)...)

# Chapter 14

pt,pe = initmesh(Rectg(0,0,1,1),"hmax"=>0.25)
force(x) = 2π*π*sin(π*x[2])*sin(π*y[3])
forse = interp(discontinuous(force.(FaceBundle(pt))))
solvepoissonDIPG(forse,-1,9) # α = SIPG parameter, β = penalty parameter
#dpt = discontinuous(pt)
#forse = interp(discontinuous(means(force.(pt))))

