# github.com/chakravala
# FEM: MacCormack

# Chapter 5

# 5.2.1 Hyperbolic explicit

icfun(x) = x[1] ≤ 1/2 ? 1.0 : 1/2 #ic = icfun.(x)
Δt(x,CFL=0.9,c=1) = CFL/c*step(points(x))

function D0(u,alg=FirstDifference{1})
    mat = Tridiagonal(u,alg)
    mat.d[1] = 0
    mat.du[1] = 0
    return mat
end
Dp(u) = D0(u,FirstDifference{0})
Dm(u) = D0(u,FirstDifference{-1})

x = TensorField(range(0,2,41))
D0x = D0(x)
dt = Δt(x,0.9,1)

# backward, forward

function explicitwave(fun,dt,f=icfun.(x))
    ic = InitialCondition(Flow(fun,9dt),f)
    odesolve(ic,ExplicitIntegrator{1}(dt))
end

back(x) = -gradient_back(localfiber(x))
forw(x) = -gradient_forw(localfiber(x))

back(x) = -Dmx*localfiber(x)
forw(x) = -Dpx*localfiber(x)

surface(explicitwave(back,dt))
surface(explicitwave(forw,dt))

# central difference

cent(x) = -gradient_fast(localfiber(x))
cent2(x) = -D0x*localfiber(x)

surface(explicitwave(cent,dt))
surface(explicitwave(cent2,dt))

# Jameson

function jamesonwave(fun,dt,f=icfun.(x))
    ic = InitialCondition(Flow(fun,9dt),f)
    odesolve(ic,ExplicitIntegrator{5}(dt))
end

surface(jamesonwave(cent2,dt))
surface(jamesonwave(cent2,Δt(x,2)))

# Lax

avg(x,i) = isone(i)||length(x)==i ? x[i] : (x[i-1]+x[i+1])/2
avg(x) = TensorField(x,[avg(fiber(x),i) for i ∈ 1:length(x)])

function laxwave(ic,dt,c=1)
    D0dt = (c*dt)*D0(ic)
    fixedpoint(x -> avg(x) - D0dt*x,ic,range(0,9dt,10))
end

surface(laxwave(icfun.(x),dt))

# Lax-Wendroff

function laxwendroffwave(ic,dt,c=1)
    D0dt = (dt*c)*D0(ic)
    fun(x) = x - D0dt*x + ((c*dt)^2)/2*gradient_forw(gradient_back(x))
    fixedpoint(fun,ic,range(0,9dt,10))
end

surface(laxwendroffwave(icfun.(x),dt))

# MacCormack

function maccormackwave(ic,dt,c=1)
    D0dt = (dt*c)*D0(ic)
    function fun(x)
        du = (-c*dt)*gradient_forw(x)
        x+(du+(-c*dt)*gradient_back(x+du))/2
    end
    fixedpoint(fun,ic,range(0,9dt,10))
end

surface(maccormackwave(icfun.(x),dt))

# Warming-Beam

function warmingbeamwave(ic,dt,c=1)
    dx = step(points(ic))
    D0dt = (dt*c)*D0(ic)
    function fun(x)
        du = ((dx-c*dt)/2)*gradient_back(x)
        x+(-c*dt)*gradient_back(x+du)
    end
    fixedpoint(fun,ic,range(0,9dt,10))
end

surface(warmingbeamwave(icfun.(x),dt))
surface(warmingbeamwave(icfun.(x),Δt(x,2)))

# upwind

upwind(x,i,c=1) = isone(i)||length(x)==i ? x[i] : (c*(x[i+1]+x[i])-abs(c)*(x[i+1]-x[i]))/2
upwind(x,c=1) = TensorField(x,[upwind(fiber(x),i,c) for i ∈ 1:length(x)])
fun(x) = -gradient_back(upwind(localfiber(x)))

surface(explicitwave(fun,dt))

# 5.2.2 implicit

function implicitwave(ic,dt)
    ID0 = I+dt*D0(ic)
    fixedpoint(x -> ID0\x,ic,range(0,9dt,10))
end
function implicitwave(ic,dt,α)
    D0x = D0(ic)
    ID0 = I+(α*dt)*D0x
    D0a = D0x*((1-α)*dt)
    fixedpoint(x -> ID0\(x-D0a*x),ic,range(0,9dt,10))
end

# implicit delta form

function implicitwave(u,dt)
    ID0 = dt*D0(u)
    ic = InitialCondition(Flow(x -> ((I+ID0)\-(ID0*x))/dt,9dt),u)
    odesolve(ic,ExplicitIntegrator{1}(dt))
end
function implicitwave(u,dt,α)
    D0x = D0(u)
    ID0 = (α*dt)*D0x
    D0a = ID0 + D0x*((1-α)*dt)
    ic = InitialCondition(Flow(x -> ((I+ID0)\-(D0a*x))/dt,9dt),u)
    odesolve(ic,ExplicitIntegrator{1}(dt))
end

# 5.2.3 implicit

surface(implicitwave(icfun.(x),Δt(x,0.9)))
surface(implicitwave(icfun.(x),Δt(x,2.0)))
surface(implicitwave(icfun.(x),Δt(x,0.9),1/2))
surface(implicitwave(icfun.(x),Δt(x,2.0),1/2))

# 5.2.4 nonlinear wave

back2(x) = -gradient_back((localfiber(x)^2)/2)

surface(explicitwave(back2,dt))

function D0u(u,d=Cartan.centraldiff_fast_points(u))
    b = fiber(u)[2:end-1]./d[2:end-1]
    c = pushfirst!(copy(b),0)
    bc = inv(points(x)[end]-points(x)[end-1])
    Tridiagonal(-push!(b,bc*fiber(u)[end-1]),push!(zeros(length(u)-1),bc*fiber(u)[end]),c)
end

function implicitwave2(ic,dt)
    d = Cartan.centraldiff_fast_points(ic)
    D0x = D0(ic)
    fixedpoint(x -> (I+D0u(dt*x,d))\(x-D0x*((dt*x/2)^2)),ic,range(0,9dt,10))
end

surface(implicitwave2(icfun.(x),Δt(x,1.8)))

# delta form

function implicitwave2(u,dt)
    d = Cartan.centraldiff_fast_points(u)
    D0x = D0(u)
    function fun(x)
        A = D0u(dt*localfiber(x),d)
        ((I+A)\-(A*localfiber(x)+D0x*((dt*localfiber(x)/2)^2)))/dt
    end
    ic = IC(Flow(fun,9dt),u)
    odesolve(ic,ExplicitIntegrator{1}(dt))
end


# 5.3 Elliptic

XY = TensorField(rakichproduct())
x,y = split(points(XY))

p∞,ρ∞,γ = 1,1,1.4
a∞ = sqrt(γ*p∞/ρ∞)
M∞ = 0.5
V∞ = M∞*a∞
A = 1-M∞^2
ϕ = V∞*getindex.(XY,1)

velocity(M∞=0.5,p∞=1,ρ∞=1,γ=1.4) = M∞*sqrt(γ*p∞/ρ∞)

struct Freestream
    M∞::Float64
    p∞::Float64
    ρ∞::Float64
    γ::Float64
    V∞::Float64
    a∞::Float64
    A::Float64
    function Freestream(M∞=0.5,p∞=1,ρ∞=1,γ=1.4,V∞=velocity(M∞,p∞,ρ∞,γ))
        new(M∞,p∞,ρ∞,γ,V∞,sqrt(γ*p∞/ρ∞),1-M∞^2)
    end
end

mach(f::Freestream) = f.M∞
pressure(f::Freestream) = f.p∞
density(f::Freestream) = f.ρ∞
velocity(f::Freestream) = f.V∞
celerity(f::Freestream) = f.a∞

function flowfield(ϕ::AbstractMatrix,δx=0.4,δy=1+4δx)
    x,y = split(points(ϕ))
    i = findall(x->abs(x-0.5)<0.5+δx,x)
    j = findall(y->abs(y)<δy,y)
    ϕ.(TensorField(ProductSpace(x[i],y[j])))
end
function flowfield(ϕ,δx=0.4,δy=3+4δx,δz=1+4δx)
    x,y,z = split(points(ϕ))
    i = findall(x->abs(x-0.5)<0.5+δx,x)
    j = findall(y->abs(y)<δy,y)
    k = findall(z->abs(z)<δz,z)
    ϕ.(TensorField(ProductSpace(x[i],y[j],z[k])))
end


function airstrip(ϕ,δx=0.3)
    x = split(points(ϕ))[1]
    i = findall(x->abs(x-0.5)<0.5+δx,x)
    ϕ[:,1].(TensorField(x[i]))
end

function initbc(u::AbstractMatrix)
    nx,ny = size(u)
    out = Array{fibertype(u),2}(undef,size(u)...)
    out[1,:] .= view(fiber(u),1,:)
    out[end,:] .= view(fiber(u),nx,:)
    out[:,end] .= view(fiber(u),:,ny)
    return out
end
function initbc(u)
    nx,ny,nz = size(u)
    out = Array{fibertype(u),3}(undef,size(u)...)
    out[1,:,:] .= view(fiber(u),1,:,:)
    out[end,:,:] .= view(fiber(u),nx,:,:)
    out[:,end] .= view(fiber(u),:,ny)
    return out
end


Dxx(x) = Tridiagonal(x,SecondDifference{1})

function elliptictri(y,A=1)
    mat = (-A)*Dxx(y)
    mat.d[1] = -1/(y[2]-y[1])
    mat.d[end] = 1
    mat.dl[end] = 0
    mat.du[1] = 1/(y[2]-y[1])
    return mat
end

function ellipticpressure(ϕ;M∞=0.5,p∞=1,ρ∞=1,γ=1.4)
    V∞ = velocity(M∞,p∞,ρ∞,γ)
    p∞*(1-(M∞^2*(γ-1)/2)*(Real(abs2(gradient(ϕ)))/V∞^2-1))^(γ/(γ-1))
end

function ellipticcpressure(ϕ;M∞=0.5,p∞=1,ρ∞=1,γ=1.4)
    V∞ = velocity(M∞,p∞,ρ∞,γ)
    (ellipticpressure(ϕ;M∞,p∞,ρ∞,γ)-p∞)*(2/(ρ∞*V∞^2))
end

function solveelliptic_pointsweep(U,out,matx,maty,bc)
    u = fiber(U)
    nx,ny = size(U)
    for i ∈ 2:nx-1
        out[i,1] = u[i,2]-bc[i]
        out[i,end] = u[i,end]
        dx1,dx2,dx3 = matx.dl[i-1],matx.d[i],matx.du[i]
        for j ∈ 2:ny-1
            out[i,j] = -(u[i-1,j]*dx1+u[i+1,j]*dx3 + u[i,j-1]*maty.dl[j-1]+u[i,j+1]*maty.du[j])/(dx2+maty.d[j])
        end
    end
    return TensorField(U,out)
end

function elliptic_vertical_rhs(u,i,matx,bc1,bc2)
    f = -view(fiber(u),i-1,:)*matx.dl[i-1]-view(fiber(u),i+1,:)*matx.du[i]
    f[1] = bc1[i]
    f[end] = bc2[i]
    return f
end

function elliptic_vertical_lhs(i,matx,maty)
    xmat = matx.d[i]*I+maty
    xmat.d[1] = maty.d[1]
    xmat.d[end] = maty.d[end]
    xmat.dl[end] = maty.dl[end]
    xmat.du[1] = maty.du[1]
    return xmat
end

function solveelliptic_vertical_sweep(u,out,matx,maty,bc1,bc2)
    for i ∈ 2:size(u)[1]-1
        f = elliptic_vertical_rhs(u,i,matx,bc1,bc2)
        xmat = elliptic_vertical_lhs(i,matx,maty)
        #out[i,:] .= xmat\f
        out[i,:] .= xmat\(f-xmat*view(fiber(u),i,:)) # delta form
    end
    return TensorField(u,out)
end

function elliptic_horizontal_rhs(u,j,maty,bc1,bc2)
    f = -view(fiber(u),:,j-1)*maty.dl[j-1]-view(fiber(u),:,j+1)*maty.du[j]
    f[1] = bc1
    f[end] = bc2
    return f
end

function elliptic_horizontal_lhs(j,matx,maty)
    ymat = maty.d[j]*I+matx
    ymat.d[1] = 1
    ymat.d[end] = 1
    ymad.dl[end] = 0
    ymat.du[1] = 0
    return ymat
end

function solveelliptic_horizontal_sweep(u,out,matx,maty,bc,bc1,bc2)
    out[:,1] .= view(fiber(u),:,2)-bc
    for j ∈ 2:size(u)[2]-1
        f = elliptic_horizontal_rhs(u,j,matx,bc1,bc2)
        ymat = elliptic_horizontal_lhs(j,matx,maty)
        #out[:,j] .= ymat\f
        out[:,j] .= ymat\(f-ymat*view(fiber(u),:,j)) # delta form
    end
    return TensorField(u,out)
end

function solveelliptic(P=CircularArc{6,21}(),k=1000;M∞=0.5,p∞=1,ρ∞=1,γ=1.4,n=51,initbc=fiber,adi=false,line=true)
    XY,ps = TensorField(rakichproduct(P,50,n)),profileslope(P)
    V∞,A = velocity(M∞,p∞,ρ∞,γ),1-M∞^2
    x,y = split(points(XY))
    matx,maty = elliptictri(x,A),elliptictri(y)
    ϕ,bc = V∞*getindex.(XY,1),V∞*ps.(x)
    if !line
        fun1(u) = solveelliptic_pointsweep(u,initbc(u),matx,maty,bc)
        fixedpoint(fun1,ϕ,k)
    elseif !adi
        fun2(u) = solveelliptic_vertical_sweep(u,initbc(u),matx,maty,bc,V∞*x)
        fixedpoint(fun2,ϕ,k)
    else
        function fun3(U)
            n,u = point(U),fiber(U)
            if iseven(n)
                solveelliptic_vertical_sweep(u,initbc(u),matx,maty,bc,V∞*x)
            else
                solveelliptic_horizontal_sweep(u,initbc(u),matx,maty,bc,V∞*x[1],V∞*x[end])
            end
            return (n+1) ↦ u
        end
        fiber(fixedpoint(fun3,0↦ϕ,k))
    end
end
function solveelliptic_point_jacobi(P=CircularArc{6,21}(),k=1000;M∞=0.5,p∞=1,ρ∞=1,γ=1.4,n=51)
    solveelliptic(P,k;M∞=M∞,p∞=p∞,ρ∞=ρ∞,γ=γ,n=n,initbc=initbc,line=false)
end
function solveelliptic_point_gauss_seidel(P=CircularArc{6,21}(),k=1000;M∞=0.5,p∞=1,ρ∞=1,γ=1.4,n=51)
    solveelliptic(P,k;M∞=M∞,p∞=p∞,ρ∞=ρ∞,γ=γ,n=n,initbc=fiber,line=false)
end
function solveelliptic_line_jacobi(P=CircularArc{6,21}(),k=1000;M∞=0.5,p∞=1,ρ∞=1,γ=1.4,n=51,adi=false)
    solveelliptic(P,k;M∞=M∞,p∞=p∞,ρ∞=1,γ=γ,n=n,initbc=initbc,adi=adi)
end
function solveelliptic_gauss_seidel(P=CircularArc{6,21}(),k=1000;M∞=0.5,p∞=1,ρ∞=1,γ=1.4,n=51,adi=false)
    solveelliptic(P,k;M∞=M∞,p∞=p∞,ρ∞=1,γ=γ,n=n,initbc=fiber,adi=adi)
end
function solveelliptic_line_adi(P=CircularArc{6,21}(),k=1000;M∞=0.5,p∞=1,ρ∞=1,γ=1.4,n=51,initbc=fiber)
    solveelliptic(P,k;M∞=M∞,p∞=p∞,ρ∞=1,γ=γ,n=n,initbc=initbc,adi=true)
end

cp = ellipticcpressure(solveelliptic(CircularArc{6,21}();n=51))
scatter(-airstrip(cp,0.3))

cp = ellipticcpressure(solveelliptic(CircularArc{6,41}();n=101))
scatter(-airstrip(cp,0.3))

cp0 = ellipticcpressure(solveelliptic(CircularArc{6,41}();n=101,M∞=0))
scatter(-airstrip(cp0,0.3)/sqrt(A))

# 3D Elliptic

function elliptic_vertical_rhs(u,i,j,matx,maty,bc1,bc2)
    f = if isone(j)
        -view(fiber(u),i-1,j,:)*matx.dl[i-1]-view(fiber(u),i+1,j,:)*matx.du[i]#-2view(fiber(u),i,j+1,:)*maty.du[j]
    else
        -view(fiber(u),i-1,j,:)*matx.dl[i-1]-view(fiber(u),i+1,j,:)*matx.du[i]-view(fiber(u),i,j-1,:)*maty.dl[j-1]-view(fiber(u),i,j+1,:)*maty.du[j]
    end
    f[1] = bc1[i,j]
    f[end] = bc2[i]
    return f
end

function elliptic_vertical_lhs(i,j,matx,maty,matz)
    xymat = if isone(j)
        matx.d[i]*I+matz
    else
        (matx.d[i]+maty.d[j])*I+matz
    end
    xymat.d[1] = matz.d[1]
    xymat.d[end] = matz.d[end]
    xymat.dl[end] = matz.dl[end]
    xymat.du[1] = matz.du[1]
    return xymat
end

function solveelliptic_vertical_sweep(u,out,matx,maty,matz,bc1,bc2)
    for j ∈ 1:size(u)[2]-1
        for i ∈ 2:size(u)[1]-1
            f = elliptic_vertical_rhs(u,i,j,matx,maty,bc1,bc2)
            xymat = elliptic_vertical_lhs(i,j,matx,maty,matz)
            #out[i,j,:] .= xymat\f
            out[i,j,:] .= xymat\(f-xymat*view(fiber(u),i,j,:)) # delta form
        end
    end
    return TensorField(u,out)
end

function solveelliptic3(P=CircularArc{6,21}(),k=1000;M∞=0.5,p∞=1,ρ∞=1,γ=1.4,n=51,initbc=fiber,adi=false,line=true)
    XYZ,ps = TensorField(rakichproduct3(P,50,n)),profileslope(P)
    V∞,A = velocity(M∞,p∞,ρ∞,γ),1-M∞^2
    x,y,z = split(points(XYZ))
    matx,maty,matz = elliptictri(x,A),elliptictri(y),elliptictri(z)
    ϕ,bc = V∞*getindex.(XYZ,1),fiber(fiberproduct(TensorField(x,V∞*ps.(x)),TensorField(y,y.≤3),*))
    if !line
        fun1(u) = solveelliptic_pointsweep(u,initbc(u),matx,maty,matz,bc)
        fixedpoint(fun1,ϕ,k)
    else
        fun2(u) = solveelliptic_vertical_sweep(u,initbc(u),matx,maty,matz,bc,V∞*x)
        fixedpoint(fun2,ϕ,k)
    end
end

# 5.5 Parabolic

Δt(x,ν=1) = (step(points(x))^2)/(2ν)

function bc!(x)
    x[1] = 0
    x[end] = 1
    return x
end

function Dpm(u) # with some boundary conditions
    mat = Dxx(u)
    mat.d[1] = 0
    mat.d[end] = 0
    mat.dl[end] = 0
    mat.du[1] = 0
    return mat
end

fun(x) = gradient_forw(gradient_back(localfiber(x)))
fun(x) = bc!(Dpmx*localfiber(x))

x = TensorField(range(0,1,41))
Dpmx = Dpm(x)
dt = Δt(x)
ic = InitialCondition(Flow(fun,100dt),hat(x,1))
out = odesolve(ic,ExplicitIntegrator{1}(dt))
surface(out)

# implicit

function implicitparabolic(ic,dt)
    IDpm = I-dt*Dpm(ic)
    fixedpoint(x -> IDpm\x,ic,0:dt:10dt)
end
function implicitparabolic(ic,dt,α)
    Dpmx = Dpm(ic)
    IDpm = I-Dpmx*(α*dt)
    Dpma = I+Dpmx*((1-α)*dt)
    fixedpoint(x -> IDpm\(Dpma*x),ic,0:dt:10dt)
end

mp = implicitparabolic(hat(x,1),1e9)
lines(mp[:,1])
lines!(mp[:,2])

surface(implicitparabolic(hat(x,1),1/8))
surface(implicitparabolic(hat(x,1),1/8,1/2))

# delta form

function implicitparabolic(u,dt)
    IDpm = I-dt*Dpm(u)
    ic = IC(Flow(x -> (IDpm\((I-IDpm)*x))/dt,10dt),u)
    odesolve(ic,ExplicitIntegrator{1}(dt))
end
function implicitparabolic(u,dt,α)
    Dpmx = Dpm(u)
    IDpm = I-Dpmx*(α*dt)
    Dpma = I+Dpmx*((1-α)*dt)-IDpm
    ic = IC(Flow(x -> (IDpm\(Dpma*x))/dt,10dt),u)
    odesolve(ic,ExplicitIntegrator{1}(dt))
end

# Chapter 6

XY = TensorField(rakichproduct())
x,y = split(points(XY))

V∞,ρ∞,γ = 1,1,1.4
M∞ = 0.908
a∞ = M∞/V∞
p∞ = ρ∞*a∞^2/γ
A = 1-M∞^2
ϕ = TensorField(XY,1.0)

pressure(M∞=0.908,V∞=1,ρ∞=1,γ=1.4) = (V∞/M∞)^2*ρ∞/γ
Freestream2(M∞=0.5,V∞=1,ρ∞=1,γ=1.4,p∞=pressure(M∞,V∞,ρ∞,γ)) = Freestream(M∞,p∞,ρ∞,γ,V∞)

function tsdmach(ϕ;M∞=0.908,V∞=1,ρ∞=1,γ=1.4)
    #p∞ = pressure(M∞,V∞,ρ∞,γ)
    #sqrt((γ*p∞/ρ∞)*(1-(M∞^2*(γ-1)/2)*(Real(abs2(gradient(ϕ)+Chain(V∞,0)))/V∞^2-1))^γ)
    sqrt(M∞^2+((γ+1)*(M∞^2)/V∞)*gradient(ϕ,1))
end

function tsdpressure(ϕ::AbstractMatrix;M∞=0.908,V∞=1,ρ∞=1,γ=1.4)
    p∞ = pressure(M∞,V∞,ρ∞,γ)
    p∞*(1-(M∞^2*(γ-1)/2)*(Real(abs2(gradient(ϕ)+Chain(V∞,0)))/V∞^2-1))^(γ/(γ-1))
end
function tsdpressure(ϕ;M∞=0.908,V∞=1,ρ∞=1,γ=1.4)
    p∞ = pressure(M∞,V∞,ρ∞,γ)
    p∞*(1-(M∞^2*(γ-1)/2)*(Real(abs2(gradient(ϕ)+Chain(V∞,0,0)))/V∞^2-1))^(γ/(γ-1))
end


function tsdcpressure(ϕ;M∞=0.908,V∞=1,ρ∞=1,γ=1.4)
    p∞ = pressure(M∞,V∞,ρ∞,γ)
    (tsdpressure(ϕ;M∞=M∞,V∞=V∞,ρ∞=ρ∞,γ=γ)-p∞)*(2/(ρ∞*V∞^2))
end

function solvetsd_pointsweep(U,out,matx,maty,bc,A,dx)
    u = fiber(U)
    nx,ny = size(U)
    A1 = Vector{fibertype(u)}(undef,size(u)[2])
    A2 = zeros(size(u)[2])
    μ1,μ2 = A2.<0,A2.<0
    for i ∈ 2:nx-1
        out[i,1] = u[i,2]-bc[i]
        out[i,end] = u[i,end]
        A1 = A.-(view(fiber(u),i+1,:)-view(fiber(u),i-1,:))/(dx[i-1]+dx[i])
        μ1 = A1.<0
        A1μ,A2μ = (1.0.-μ1).*A1,μ2.*A2
        dx1,dx2,dx3 = matx.dl[i-1],matx.d[i],matx.du[i]
        if i == 2
             for j ∈ 2:ny-1
                 out[i,j] = -((u[i-1,j]*dx1+u[i+1,j]*dx3)*A1μ[j] + u[i,j-1]*maty.dl[j-1]+u[i,j+1]*maty.du[j])/(dx2*A1μ[j]+maty.d[j])
            end
        else
            dz1,dz2,dz3 = matx.dl[i-2],matx.d[i-1],matx.du[i-1]
            for j ∈ 2:ny-1
                out[i,j] = -((u[i-1,j]*dx1+u[i+1,j]*dx3)*A1μ[j] + (u[i-1,j]*dz2+u[i-2,j]*dz3)*A2μ[j] + u[i,j-1]*maty.dl[j-1]+u[i,j+1]*maty.du[j])/(dx2*A1μ[j]+dz1*A2μ[j]+maty.d[j])
            end
        end
        A2,μ2 = A1,μ1
    end
    return TensorField(U,out)
end

function tsd_vertical_rhs(u,i,matx,bc1,bc2,A1μ,A2μ)
    f = if i == 2
        (-view(fiber(u),i-1,:).*matx.dl[i-1]-view(fiber(u),i+1,:).*matx.du[i]).*A1μ
    else
        (-view(fiber(u),i-1,:).*matx.dl[i-1]-view(fiber(u),i+1,:).*matx.du[i]).*A1μ .+(-view(fiber(u),i-1,:).*matx.d[i-1]-view(fiber(u),i-2,:).*matx.du[i-1]).*A2μ
    end
    f[1] = bc1[i]
    f[end] = bc2[i]
    return f
end

function tsd_vertical_lhs(i,matx,maty,A1μ,A2μ)
    xmat = if i == 2
        Diagonal(matx.d[i]*A1μ)+maty
    else
        Diagonal(matx.d[i]*A1μ+matx.dl[i-2]*A2μ)+maty
    end
    xmat.d[1] = maty.d[1]
    xmat.d[end] = maty.d[end]
    xmat.dl[end] = maty.dl[end]
    xmat.du[1] = maty.du[1]
    return xmat
end

function solvetsd_vertical_sweep(u,out,matx,maty,bc1,bc2,A,dx)
    A1 = Vector{fibertype(u)}(undef,size(u)[2])
    A2 = zeros(size(u)[2])
    μ1,μ2 = A2.<0,A2.<0
    for i ∈ 2:size(u)[1]-1
        A1 = A.-(view(fiber(u),i+1,:)-view(fiber(u),i-1,:))/(dx[i-1]+dx[i])
        μ1 = A1.<0
        A1μ,A2μ = (1.0.-μ1).*A1,μ2.*A2
        f = tsd_vertical_rhs(u,i,matx,bc1,bc2,A1μ,A2μ)
        xmat = tsd_vertical_lhs(i,matx,maty,A1μ,A2μ)
        #out[i,:] .= xmat\f
        out[i,:] .+= xmat\(f-xmat*view(fiber(u),i,:)) # delta form
        A2,μ2 = A1,μ1
    end
    return TensorField(u,out)
end

function solvetsd(P=CircularArc{6,21}(),k=400;M∞=0.908,V∞=1,ρ∞=1,γ=1.4,n=51,initbc=fiber,line=true)
    XY,ps = TensorField(rakichproduct(P,50,n)),profileslope(P)
    p∞,A = pressure(M∞,V∞,ρ∞,γ),1-M∞^2
    x,y = split(points(XY))
    dx,bc2 = diff(x)/((γ+1)*M∞^2/V∞),ones(length(x))
    matx,maty = elliptictri(x),elliptictri(y)
    ϕ,bc1 = TensorField(XY,1.0),V∞*ps.(x)
    if !line
        fun1(u) = solvetsd_pointsweep(u,initbc(u),matx,maty,bc1,A,dx)
        fixedpoint(fun1,ϕ,k)
    else
        fun2(u) = solvetsd_vertical_sweep(u,initbc(u),matx,maty,bc1,bc2,A,dx)
        fixedpoint(fun2,ϕ,k)
    end
end

function solvetsd_point_jacobi(P=CircularArc{6,21}(),k=400;M∞=0.908,V∞=1,ρ∞=1,γ=1.4,n=51)
    solvetsd(P,k;M∞=M∞,V∞=V∞,ρ∞=ρ∞,γ=γ,n=n,initbc=initbc,line=false)
end
function solvetsd_point_gauss_seidel(P=CircularArc{6,21}(),k=400;M∞=0.908,V∞=1,ρ∞=1,γ=1.4,n=51)
    solvetsd(P,k;M∞=M∞,V∞=V∞,ρ∞=ρ∞,γ=γ,n=n,initbc=fiber,line=false)
end
function solvetsd_line_jacobi(P=CircularArc{6,21}(),k=400;M∞=0.908,V∞=1,ρ∞=1,γ=1.4,n=51)
    solvetsd(P,k;M∞=M∞,V∞=V∞,ρ∞=ρ∞,γ=γ,n=n,initbc=initbc,line=true)
end
function solvetsd_gauss_seidel(P=CircularArc{6,21}(),k=400;M∞=0.908,V∞=1,ρ∞=1,γ=1.4,n=51)
    solvetsd(P,k;M∞=M∞,V∞=V∞,ρ∞=ρ∞,γ=γ,n=n,initbc=fiber,line=true)
end

cp = cpressure(gauss_seidel(CircularArc{6,21}();n=51))
scatter(-airstrip(cp,0.3))

out = gauss_seidel(CircularArc{6,41}();n=101,M∞=0.908)
cp = cpressure(gauss_seidel(CircularArc{6,41}();n=101,M∞=0.908),M∞=0.908)
scatter(-airstrip(cp,0.3))
contour(flowfield(mach(out)),levels=[1,1.1,1.2,1.3])
streamplot!(flowfield(gradient(out)+Chain(1,0)))

cp0 = cpressure(gauss_seidel(CircularArc{6,41}();n=101,M∞=0))
scatter(-airstrip(cp0,0.3)/sqrt(A))

cp = cpressure(gauss_seidel(CircularArc{6,41}();n=101,M∞=0.908),M=0.908)
scatter(-airstrip(cp,0.3))

# 3D TSD

function tsd_vertical_rhs(u,i,j,matx,maty,bc1,bc2,A1μ,A2μ)
    f = if isone(j) && i == 2
        (-view(fiber(u),i-1,j,:).*matx.dl[i-1]-view(fiber(u),i+1,j,:).*matx.du[i]).*A1μ
    elseif isone(j)
        (-view(fiber(u),i-1,j,:).*matx.dl[i-1]-view(fiber(u),i+1,j,:).*matx.du[i]).*A1μ .+(-view(fiber(u),i-1,j,:).*matx.d[i-1]-view(fiber(u),i-2,j,:).*matx.du[i-1]).*A2μ
    elseif i == 2
        (-view(fiber(u),i-1,j,:).*matx.dl[i-1]-view(fiber(u),i+1,j,:).*matx.du[i]).*A1μ-view(fiber(u),i,j-1,:)*maty.dl[j-1]-view(fiber(u),i,j+1,:)*maty.du[j]
    else
        (-view(fiber(u),i-1,j,:).*matx.dl[i-1]-view(fiber(u),i+1,j,:).*matx.du[i]).*A1μ .+(-view(fiber(u),i-1,j,:).*matx.d[i-1]-view(fiber(u),i-2,j,:).*matx.du[i-1]).*A2μ-view(fiber(u),i,j-1,:)*maty.dl[j-1]-view(fiber(u),i,j+1,:)*maty.du[j]
    end
    f[1] = bc1[i,j]
    f[end] = bc2[i]
    return f
end

function tsd_vertical_lhs(i,j,matx,maty,matz,A1μ,A2μ)
    xymat = if isone(j) && i == 2
        Diagonal(matx.d[i]*A1μ)+matz
    elseif isone(j)
        Diagonal(matx.d[i]*A1μ+matx.dl[i-2]*A2μ)+matz
    elseif i == 2
        (Diagonal(matx.d[i]*A1μ)+I*maty.d[j])+matz
    else
        (Diagonal(matx.d[i]*A1μ+matx.dl[i-2]*A2μ)+I*maty.d[j])+matz
    end
    xymat.d[1] = matz.d[1]
    xymat.d[end] = matz.d[end]
    xymat.dl[end] = matz.dl[end]
    xymat.du[1] = matz.du[1]
    return xymat
end

function solvetsd_vertical_sweep(u,out,matx,maty,matz,bc1,bc2,A,dx)
    A1 = Vector{fibertype(u)}(undef,size(u)[2])
    A2 = zeros(size(u)[2])
    μ1,μ2 = A2.<0,A2.<0
    for j ∈ 1:size(u)[2]-1
        for i ∈ 2:size(u)[1]-1
            A1 = A.-(view(fiber(u),i+1,j,:)-view(fiber(u),i-1,j,:))/(dx[i-1]+dx[i])
            μ1 = A1.<0
            A1μ,A2μ = (1.0.-μ1).*A1,μ2.*A2
            f = tsd_vertical_rhs(u,i,j,matx,matx,bc1,bc2,A1μ,A2μ)
            xymat = tsd_vertical_lhs(i,j,matx,maty,matz,A1μ,A2μ)
            #out[i,j,:] .= xymat\f
            out[i,j,:] .+= xymat\(f-xymat*view(fiber(u),i,j,:)) # delta form
            A2,μ2 = A1,μ1
        end
    end
    return TensorField(u,out)
end

function solvetsd3(P=CircularArc{6,21}(),k=400;M∞=0.908,V∞=1,ρ∞=1,γ=1.4,n=51,initbc=fiber,line=true)
    XYZ,ps = TensorField(rakichproduct3(P,50,n)),profileslope(P)
    p∞,A = pressure(M∞,V∞,ρ∞,γ),1-M∞^2
    x,y,z = split(points(XYZ))
    dx,bc2 = diff(x)/((γ+1)*M∞^2/V∞),ones(length(x))
    matx,maty,matz = elliptictri(x),elliptictri(y),elliptictri(z)
    ϕ,bc1 = TensorField(XYZ,1.0),fiber(fiberproduct(TensorField(x,V∞*ps.(x)),TensorField(y,y.≤3),*))
    if !line
        fun1(u) = solvetsd_pointsweep(u,initbc(u),matx,maty,matz,bc1,A,dx)
        fixedpoint(fun1,ϕ,k)
    else
        fun2(u) = solvetsd_vertical_sweep(u,initbc(u),matx,maty,matz,bc1,bc2,A,dx)
        fixedpoint(fun2,ϕ,k)
    end
end

# Chapter 8

function flowfield2(ϕ::AbstractMatrix,δx=0.4,δy=1+4δx)
    x,y = split(points(ϕ))
    i = findall(x->abs(x-0.5)<0.5+δx,x)
    j = findall(y->abs(y)<δy,y)
    ϕ.(TensorField(ProductSpace(x[i],y[j])))
end

function airstrip(ϕ,δx=0.3)
    x = split(points(ϕ))[1]
    i = findall(x->abs(x-0.5)<0.5+δx,x)
    ϕ[:,1].(TensorField(x[i]))
end


function fppressure(ϕ;M∞=0.5,p∞=1,ρ∞=1,γ=1.4)
    V∞ = velocity(M∞,p∞,ρ∞,γ)
    p∞*(1-(M∞^2*(γ-1)/2)*(Real(abs2(gradient(ϕ)))/V∞^2-1))^(γ/(γ-1))
end

function fpcpressure(ϕ;M∞=0.5,p∞=1,ρ∞=1,γ=1.4)
    V∞ = velocity(M∞,p∞,ρ∞,γ)
    (ellipticpressure(ϕ;M∞,p∞,ρ∞,γ)-p∞)*(2/(ρ∞*V∞^2))
end

function fpdensity(ϕ,iJ,dxy=centraldiff(points(ϕ));M∞=0.5,p∞=1,ρ∞=1,γ=1.4)
    V∞ = velocity(M∞,p∞,ρ∞,γ)
    uv = iJ > gradient(ϕ,dxy)
    ρ∞*(1-(M∞^2*(γ-1)/2)*(Real(abs2(uv))/V∞^2-1))^(1/(γ-1))
end

function fp_center_rhs(u,i,ρim,ρip,ρjm,ρjp,T12,T21)
    nx,ny = size(u)
    (ρip.*view(T12,i+1,2:ny-1) .+ ρjp.*view(T21,i,3:ny)).*view(fiber(u),i+1,3:ny) .-
    (ρip.*view(T12,i+1,2:ny-1) .+ ρjm.*view(T21,i,1:ny-2)).*view(fiber(u),i+1,1:ny-2) .-
    (ρim.*view(T12,i-1,2:ny-1) .+ ρjp.*view(T21,i,3:ny)).*view(fiber(u),i-1,3:ny) .+
    (ρim.*view(T12,i-1,2:ny-1) .+ ρjm.*view(T21,i,1:ny-2)).*view(fiber(u),i-1,1:ny-2)
end

function fp_vertical_rhs(u,i,ρim,ρip,ρjm,ρjp,T11,T12,T21,T22,bc)
    nx,ny = size(u)
    f = (ρim.*(view(T11,i,1:ny-2).+view(T11,i,2:ny-1))./2 .+ (ρjm.-ρjp).*view(T21,i,2:ny-1)./4).*view(fiber(u),i-1,2:ny-1) .+
    (ρip.*(view(T11,i,2:ny-1).+view(T11,i,3:ny))./2 .+ (ρjp.-ρjm).*view(T21,i,2:ny-1)./4).*view(fiber(u),i+1,2:ny-1) .+
    fp_center_rhs(u,i,ρim,ρip,ρjm,ρjp,T12,T21)
    push!(pushfirst!(f,(fiber(u)[i+1,1]-fiber(u)[i-1,1])*T21[i,1]/T22[i,1]/2),bc[i])
end

function fp_vertical_lhs(u,i,ρim,ρip,ρjm,ρjp,T11,T12,T21,T22)
    nx,ny = size(u)
    T22m = ρjm.*(view(T22,i,1:ny-2).+view(T22,i,2:ny-1))./2
    T22p = ρjp.*(view(T22,i,2:ny-1).+view(T22,i,3:ny))./2
    b = (ρip.-ρim).*view(T12,i,2:ny-1)./4 .- T22m
    a = ρim.*(view(T11,i,1:ny-2).+view(T11,i,2:ny-1))./2 .+ ρip.*(view(T11,i,2:ny-1).+view(T11,i,3:ny))./2 .+ T22m .+ T22p
    c = (ρim.-ρip).*view(T12,i,2:ny-1)./4 .- T22p
    Tridiagonal(push!(b,0),pushfirst!(push!(a,1),-1),pushfirst!(c,1))
end

function solvefp_vertical_sweep(u,out,iJ,dxy,bc,T11,T12,T21,T22)
    nx,ny = size(u)
    ρ = fpdensity(u,iJ,dxy)#;M∞=M∞,p∞=p∞,ρ∞=ρ∞,γ=γ)
    for i ∈ 2:size(u)[1]-1
        ρim = (view(fiber(ρ),i-1,2:ny-1).+view(fiber(ρ),i,2:ny-1))./2
        ρip = (view(fiber(ρ),i,2:ny-1).+view(fiber(ρ),i+1,2:ny-1))./2
        ρjm = (view(fiber(ρ),i,1:ny-2).+view(fiber(ρ),i,2:ny-1))./2
        ρjp = (view(fiber(ρ),i,2:ny-1).+view(fiber(ρ),i,3:ny))./2
        f = fp_vertical_rhs(u,i,ρim,ρip,ρjm,ρjp,T11,T12,T21,T22,bc)
        xmat = fp_vertical_lhs(u,i,ρim,ρip,ρjm,ρjp,T11,T12,T21,T22)
        #out[i,:] .= xmat\f
        out[i,:] .+= xmat\(f-xmat*view(fiber(u),i,:)) # delta form
    end
    return TensorField(u,out)
end

function solvefp(P=CircularArc{6,21}(),k=400;M∞=0.5,p∞=1,ρ∞=1,γ=1.4,n=51,initbc=fiber,line=true,adi=false)
    XY = rakichfield(P,50,n)
    V∞ = velocity(M∞,p∞,ρ∞,γ)
    dxy = centraldiff(points(XY))
    J = jacobian_slow(XY)
    iJ = inv(J)
    T = cofactor(J)*iJ
    T11,T21 = fiber.(split(getindex.(T,1)))
    T12,T22 = fiber.(split(getindex.(T,2)))
    ϕ = V∞*getindex.(XY,1)
    bc = fiber(ϕ)[:,end]
    if !line
        fun1(u) = solvefp_pointsweep(u,initbc(u),matx,maty,bc1,A,dx)
        fixedpoint(fun1,ϕ,k)
    elseif !adi
        fun2(u) = solvefp_vertical_sweep(u,initbc(u),iJ,dxy,bc,T11,T12,T21,T22)
        fixedpoint(fun2,ϕ,k)
    else
        function fun3(U)
            n,u = point(U),fiber(U)
            if iseven(n)
                solvefp_vertical_sweep(u,initbc(u),matx,maty,bc,V∞*x)
            else
                solvefp_horizontal_sweep(u,initbc(u),matx,maty,bc,V∞*x[1],V∞*x[end])
            end
            return (n+1) ↦ u
        end
        fiber(fixedpoint(fun3,0↦ϕ,k))
    end
end

# Chapter 9

function euler_rfft(U,γ=1.4)
    ρ,ρu,e = split(localfiber(U))
    u = ρu/ρ
    ε = (e/ρ) - (u^2)/2
    p = (γ-1)*ρ*ε
    -Chain.(gradient_rfft(ρu),gradient_rfft((ρu)^2/ρ+p),gradient_rfft((e+p)*u))
end



function eulerF(U,γ=1.4)
    ρ,ρu,e = split(localfiber(U))
    u = ρu/ρ
    ε = (e/ρ) - (u^2)/2
    p = (γ-1)*ρ*ε
    Chain.(ρu,(ρu)^2/ρ+p,(e+p)*u)
end

euler_ic(ρ,p,u,γ=1.4) = Chain.(ρ,ρ*u,p/(γ-1) + (u^2)*ρ/2)
euler1(U,γ=1.4) = gradient_forw(eulerF(localfiber(U),γ))
Δt(x,CFL=0.9,c=1) = CFL/c*step(points(x))

x = TensorField(range(0,2,81))
ρ0 = heaviside(x-1)+1
p0 = heaviside(x-1)+1


# Lax-Wendroff

function laxwendroffwave(ic,dt,c=1)
    D0dt = (dt*c)*D0(ic)
    fun(x) = x - D0dt*x + ((c*dt)^2)/2*gradient_forw(gradient_back(x))
    fixedpoint(fun,ic,range(0,9dt,10))
end

surface(laxwendroffwave(icfun.(x),dt))

# MacCormack

function maccormackeuler1(ic,dt)
    function fun(x)
        du = (-dt)*gradient_forw(eulerF(x))
        x+(du+(-dt)*gradient_back(eulerF(x+du)))/2
    end
    fixedpoint(fun,ic,0:dt:100dt)
end

surface(maccormackeuler1(euler_ic(ρ0,p0,0x),dt))

# Warming-Beam

function warmingbeamwave(ic,dt,c=1)
    dx = step(points(ic))
    D0dt = (dt*c)*D0(ic)
    function fun(x)
        du = ((dx-c*dt)/2)*gradient_back(x)
        x+(-c*dt)*gradient_back(x+du)
    end
    fixedpoint(fun,ic,range(0,9dt,10))
end

surface(warmingbeamwave(icfun.(x),dt))
surface(warmingbeamwave(icfun.(x),Δt(x,2)))

# Jameson

function rffteuler(fun,dt,f=euler_ic(ρ0,p0,0*x))
    ic = InitialCondition(Flow(fun,100dt),f)
    odesolve(ic,ExplicitIntegrator{1}(dt))
end


function jamesoneuler(fun,dt,f=euler_ic(ρ0,p0,0*x))
    ic = InitialCondition(Flow(fun,100dt),f)
    odesolve(ic,ExplicitIntegrator{5}(dt))
end



