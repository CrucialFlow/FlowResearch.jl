# github.com/chakravala
# Phys: Classical Mechanics, Hasbun

# Chapter 1

# ho1.m

function init_mass_spring(;tmax=1,k=1000,m=5,C=0,E=0,ω=0,x0=0.1,v0=0,t0=0)
    k # 1000 # spring constant
    m # 5.0 # bob mass
    C # 0.0 # damping coefficient
    E # 0.0 # driving force magnitude
    ω # 0.0 # driving force frequency
    x0 # 0.1 # initial position
    v0 # 0.0 # initial velocity
    t0 # 0.0 # initial time
    F0 = -k*x0-C*v0+E*sin(ω*t0) # initial force
    a0 = F0/m # initial acceleration
    InitialCondition(mass_spring(k,C,F0,ω,m),t0↦Chain(x0,v0),tmax)
end

mass_spring(k,C,F0,ω,m) = x -> Chain(x[2],(F0*sin(ω*point(x))-k*x[1]-C*x[2])/m)
function solve_mass_spring(;n=100,tmax=1,k=1000,m=5,C=0,E=0,ω=0,x0=0.1,v0=0,t0=0)
    ic = init_mass_spring(;tmax,k,m,C,E,ω,x0,v0,t0)
    getindex.(odesolve(ic,ExplicitIntegrator{4}(tmax/n)),1)
end

solve_mass_spring(n=100,tmax=1,k=1000,m=5,C=0,E=0,ω=0,x0=0.1,v0=0,t0=0)
solve_mass_spring(n=100,tmax=10,k=1,m=1,C=0,E=0,ω=0,x0=1,v0=0,t0=0)
solve_mass_spring(n=200,tmax=20,k=1,m=1,C=0.5,E=0,ω=0,x0=1,v0=0,t0=0)
solve_mass_spring(n=200,tmax=20,k=1,m=1,C=0.5,E=0.1,ω=0.8,x0=1,v0=0,t0=0)

# ho2.m

function init_free_fall(;tmax=20,g=9.8,m=1,C=0.05,y0=0,v0=110,flag=false)
    g # gravity
    m # mass
    C # drag coefficient
    y0 # initial height
    v0 # initial velocity
    t0 = 0 # initial time
    F = flag ? -m*g-C*v0*abs(v0) : -m*g-C*v0 # initial force
    vt = flag ? sqrt(m*g/C) : abs(m*g/C) # terminal velocity
    a0 = F/m # initial acceleration
    InitialCondition(free_fall(g,m,C,flag),float(t0)↦Chain(float(y0),float(v0)),tmax)
end

function free_fall(g,m,C,flag)
    flag ? x-> Chain(x[2],-C*x[2]*abs(x[2])/m-g) : x -> Chain(x[2],-C*x[2]/m-g)
end
function solve_free_fall(;n=200,tmax=20,g=9.8,m=1,C=0.05,y0=0,v0=110,flag=false)
    ic = init_free_fall(;tmax,g,m,C,y0,v0,flag)
    getindex.(odesolve(ic,ExplicitIntegrator{4}(tmax/n)))
end

solve_free_fall(n=200,tmax=20,g=9.8,m=1,C=0.05,y0=0,v0=110,flag=false)
solve_free_fall(n=200,tmax=10,g=9.8,m=1,C=0.05,y0=0,v0=110,flag=true)

# Chapter 2

# foft.m

function init_acceleration(;m=1,F0=1,ω=3)
    t = TensorField(0:0.1:10)
    a = F0*cos(ω*t)/m
end

function solve_position(;m=1,F0=1,ω=3,x0=0,v0=0.05)
    x0 + integral(v0+integral(init_acceleration(;m,F0,ω)))
end

solve_position(m=1,F0=1,ω=3,x0=0,v0=0.05)

# fofx.m

function fofx(;m=1,k=0.01,x0=0,v0=0.5)
    m # mass
    k # spring constant
    ω = sqrt(k/m) # natural frequency
    x0 # initial position
    v0 # initial velocity
    t = TensorField(0:0.05:2π/ω)
    E0 = (m/2)*v0^2 # total initial energy
    x = sqrt(2*E0/k)*sin(ω*t) # position vs time
    v = v0*cos(ω*t) # velocity vs time
    # a = (-k/m)*x # acceleration vs time
    PE = (k/2)*x^2 # potential energy
    KE = (m/2)*v^2 # kinetic energy
    E = PE + KE
    F = -k*x # force
    display(lines(Chain.(x,PE)))
    lines!(Chain.(x,KE))
    lines!(Chain.(x,E))
    lines!(Chain.(x,F))
end

# fofv.m

function fofv(;g=9.8,m=1,C=0.05,y0=10,v0=20,NPTS=100,tmax=4.5)
    g # gravity
    m # mass
    C # drag coefficient
    C < 1e-3 && (C=1e-3)
    y0 # initial height
    v0 # initial velocity
    NPTS # number of  points
    tmax
    tz = v0/g + sqrt((v0/g)^2+2*y0/g)
    f(x) = y0-m*(g*t+(m*g/C+v0)*(exp(-C*t/m)-1))/C
    t = TensorField(0:tmax/NPTS:tmax)
    y = y0-m*(g*t+(m*g/C+v0)*(exp(-C*t/m)-1))/C
    v = (m*g/C+v0)*exp(-C*t/m)-m*g/C
    a = -g-C*v/m
    display(lines(y))
    lines!(v)
    lines!(a)
end

# Chapter 3

# xoft.m

function xoft(;n=100,τ=2,r=2)
    n # number of points
    τ # period
    r # radius
    t = TensorField(0:2/n:2τ)
    x = r*cos(2π*t/τ)
end

lines(xoft(n=100,τ=2,r=2))

# Example 3.1

lines(tangent(xoft(n=100,τ=2,r=2)))

# Example 3.2

function v_and_f(;xb=3/2,vmin=-4/27,xmin=0.5,n=100)
    xmax = 5xb
    x = TensorField(xmin:2/n:2xmax)
    V = 1/x^3 - 1/x^2
    F = 3/x^4 - 2/x^3
    st = lines(V)
    fig,ax,plt = st
    lines!(F)
    ax.limits = ((1,7),(-0.2,0.2))
    st
end

function v_and_f(;xb=3/2,vmin=-4/27,xmin=0.5,n=100)
    xmax = 5xb
    x = TensorField(xmin:2/n:2xmax)
    V = 1/x^3 - 1/x^2
    F = 3/x^4 - 2/x^3
    return V,F
end

lines(bound(V,0.2))
lines!(bound(F,0.2))

# over_crit_damp

function over_crit_damp(;m=0.05,k=1,c=0.5,x0=1,v0=5,tmax=2,n=100)
    γ = c/2/m
    desc = γ^2-k/m
    desc ≤ 0 && throw("γ needs to be smaller")
    γ1 = γ+sqrt(desc)
    γ2 = γ-sqrt(desc)
    Bo = (v0+γ1*x0)/(γ1-γ2)
    Ao = x0-Bo
    Ac = x0
    Bc = v0+γ*x0
    t = TensorField(0:tmax/n:tmax) # time
    xo = Ao*exp(-γ1*t) + Bo*exp(-γ2*t) # overdamped
    xc = Ac*exp(-γ*t) + Bc*t*exp(-γ*t) # critically damped
    return (xo,xc)
end

xo,xc = over_crit_damp(m=0.05,k=1,c=0.5,x0=1,v0=5,tmax=2,n=100)
lines(xo)
lines!(xc)

# under_damp

function under_damp(;m=0.05,k=1,c=0.08,x0=1,v0=5,tmax=5,n=100)
    γ = c/2/m
    ω0 = sqrt(k/m)
    desc = ω0^2-γ^2
    desc ≤ 0 && throw("γ needs to be smaller")
    ω = sqrt(desc)
    B = sqrt(x0^2+(v0+γ*x0)^2/ω^2)
    θ = atan(ω*x0/(v0+γ*x0))
    t = TensorField(0:tmax/n:tmax)
    xe = B*exp(-γ*t)
    x = xe*sin(ω*t+θ)
end

under_damp(m=0.05,k=1,c=0.08,x0=1,v0=5,tmax=5,n=100)

# dive_amp

function drive_amp(;m=0.5,k=0.5,F0=0.5,ωmin=0.1,ωmax=2,n=200,cmin=0.2)
    ωo = sqrt(k/m) # SHO natural frequency
    ω = TensorField(ωmin:ωmax/n:ωmax)
    cmax = m*ωo*(2/sqrt(2))
    cstep = (cmax-cmin)/5
    c = cmin:cstep:cmax
    for i ∈ 1:length(c) # loop over drag coefficient
        γ = c[i]/2/m
        desc = (2γ*ω)^2+(ωo^2-ω^2)^2
        A = F0/m/sqrt(desc) # driven ho amplitude
        out = Chain.(ω+0,A)
        isone(i) ? display(lines(out)) : lines!(out)
        om_res = sqrt(ωo^2-2γ^2) # resonant frequency
        Amax = F0/2/m/γ/sqrt(ωo^2-γ^2) # maximum amplitude at om_res
        scatter!([Chain(om_res,Amax)])
    end
end

drive_amp(m=0.5,k=0.5,F0=0.5,ωmin=0.1,ωmax=2,n=200,cmin=0.2)

# drive_phase

function drive_phase(;m=0.5,k=0.5,ωmin=0,n=2.5,cmin=0.01,cmax=1)
    ωo = sqrt(k/m)
    ωmax = n*ωo
    NPTS = Int(round(33n+1))
    ωstep = (ωmax-ωmin)/NPTS
    cstep = (cmax-cmin)/5
    ω = ωmin:ωstep:ωmax
    ϕ = Vector{Float64}(undef,length(ω))
    for c ∈ cmin:cstep:cmax
        γ = c/2/m
        for i ∈ 1:length(ω)
            ωi = ω[i]
            den = ωo^2-ωi^2
            ϕi = atan((2γ/den)*ωi)
            ϕ[i] = ωi≤ωo ? ϕi : π+ϕi
        end
        out = TensorField(ω,Chain.(ω,ϕ))
        c==cmin ? display(lines(out)) : lines!(out)
    end
end

drive_phase(m=0.5,k=0.5,ωmin=0,n=2.5,cmin=0.01,cmax=1)

# drive_sol

function drive_sol(ω,c;x0=1,v0=5,m=0.5,k=0.5,F0=0.5,θ=0,dt=0.05)
    ωo = sqrt(k/m)
    τ = 2π/ω
    tmax = 5τ
    n = tmax/dt
    γ = c/2/m
    desc = (2γ*ω)^2+(ωo^2-ω^2)^2
    A = F0/m/sqrt(desc)
    den = ωo^2-ω^2
    iszero(den) && (den = 1e-3)
    ϕ = ω≤ωo ? atan((2γ/den)*ω) : π+atan((2γ/den)*ω)
    δ = θ-ϕ # forced solutions phase
    t = TensorField(0:tmax/n:tmax)
    xf = A*cos(ω*t+δ)
    F = F0*cos(ω*t+θ)
    # homoegeneous + forced solution
    desc = ωo-γ^2
    desc ≤ 0 && throw("γ needs to be smaller")
    ωu = sqrt(desc)
    th = atan(ωu*x0/(v0+γ*x0))
    B = sqrt(x0^2+(v0+γ*x0)^2/ωu^2)
    xh = B*exp(-γ*t)*sin(ωu*t+th)
    return xf,F,xf+xh
end

xf,F,xh = drive_sol(0.1,0.1;x0=1,v0=5,m=0.5,k=0.5,F0=0.5,θ=0,dt=0.05)
lines(xf)
lines!(F)
lines(xh)

# drive_power

function drive_power(;m=0.5,k=0.5,F0=0.5,ωmin=0.01,ωmax=3,n=200,cmin=0.2,cmax=1)
    ωo = sqrt(k/m)
    dw = (ωmax-ωmin)/(n-1)
    cstep = (cmax-cmin)/3
    c = cmin:cstep:cmax
    ω = TensorField(ωmin:dw:ωmax)
    power = Vector{Float64}(undef,length(ω))
    for k ∈ 1:length(c)
        γ = c[k]/2/m
        for i ∈ 1:n
            ωi = fiber(ω)[i]
            desc = (2γ*ωi)^2+(ωo^2-ωi^2)^2
            A = F0/m/sqrt(desc) # driven ho amplitude
            den = ωo^2-ωi^2
            iszero(den) && (den = 1e-3)
            ϕ = ωi≤ωo ? atan((2γ/den)*ωi) : π+atan((2γ/den)*ωi)
            power[i] = 0.5F0*A*ωi*sin(ϕ)
        end
        out = Chain.(ω,TensorField(ω,power))
        isone(k) ? display(lines(out)) : lines!(out)
    end
end

drive_power(m=0.5,k=0.5,F0=0.5,ωmin=0.01,ωmax=3,n=200,cmin=0.2,cmax=1)

# Chapter 5

# gradient_ex

function gradient_ex(;vmax=2,vs=0.1,dv=0.1)
    xmax,ymax,zmax = vmax,vmax,vmax
    xs,ys,zs = vs,vs,vs
    N = Int(round(2vmax/vs))
    dx,dy,dz = dv,dv,dv
    m = round(N/2+5)
    zm = -zmax+(m-1)*zs # value of z at which we plot f(x,y,z)
    x,y,z = -xmax:xs:xmax,-ymax:ys:ymax,-zmax:zs:zmax
    xyz = TensorField(ProductSpace(x,y,z))
end

fun(x) = x[1]*exp(-x[1]^2-x[2]^2-x[3]^2)
f = fun.(gradient_ex(vmax=2,vs=0.1,dv=0.1))
df = gradient(f)
dx,dy,dz = getindex.(df,1),getindex.(df,2),getindex.(df,3)
dxy = Chain.(dx,dy)

streamplot(df)
surface(leaf(f,0.4,3))
streamplot(leaf(dxy,0.4,3))
contour!(leaf(f,0.4,3),levels=-0.3:0.1:0.3)

# divergence_ex

function divergence_ex(;vmax=2,vs=0.1,dv=0.1)
    xmax,ymax,zmax = vmax,vmax,vmax
    xs,ys,zs = vs,vs,vs
    N = Int(round(2vmax/vs))
    dx,dy,dz = dv,dv,dv
    m = round(N/2+5)
    zm = -zmax+(m-1)*zs # value of z at which we plot f(x,y,z)
    x,y,z = -xmax:xs:xmax,-ymax:ys:ymax,-zmax:zs:zmax
    xyz = TensorField(ProductSpace(x,y,z))
end

fun(x) = x[1]*exp(-x[1]^2-x[2]^2-x[3]^2)
f = fun.(gradient_ex(vmax=2,vs=0.1,dv=0.1))
df = gradient(f)
dx,dy,dz = getindex.(df,1),getindex.(df,2),getindex.(df,3)
dxy = Chain.(dx,dy)

streamplot(df)
surface(leaf(f,0.4,3))
streamplot(leaf(dxy,0.4,3))
contour!(leaf(f,0.4,3),levels=-0.3:0.1:0.3)


