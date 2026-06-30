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

# Chapter 4

# inter_spr1

function inter_spr1(;m1=1,m2=2,k0=0.5,x10=1,x20=-1,v10=0.02,v20=0.04)
    μ = m1*m2/(m1+m2) # reduced mass
    xcm0 = (m1*x10+m2*x20)/(m1+m2) # initial center of mass
    vcm = (m1*v10+m2*v20)/(m1+m2) # center of mass speed
    xr0,vr0 = x20-x10,v20-v10 # relative coordinate / speed
    om = sqrt(k0/μ) # frequency
    τ = 2pi/om # period
    A,B = vr0/om,xr0 # amplitudes
    tmax = 2τ # time range
    t = TensorField(0:τ/50:tmax)
    xr = A*sin(om*t)+B*cos(om*t) # solution
    xcm = xcm0 + vcm*t # cm position vs time
    x1 = xcm-m2*xr/(m1+m2) # mass position vs time
    x2 = xcm+m1*xr/(m1+m2)
    lines([xcm,x1,x2])
end

# eigen

M = TensorOperator(Chain(1,-1),Chain(-1,1))
P = eigvecs(M)
L = eigvals(M)
P/P
P\M*P

# inter_spr2

function inter_spr2(;m=1,k0=1.0,k=10.0,x10=1,x20=0)
    xs,xd = (x10+x20)/2,(x10-x20)/2
    om1,om2 = sqrt(k/m),sqrt((k+2k0)/m)
    om = min(om1,om2)
    τ = 2pi/om # period
    tmax = 20τ # time range
    t = TensorField(0:τ/50:tmax)
    x1 = xs*cos(om1*t)+xd*cos(om2*t)
    x2 = xs*cos(om1*t)-xd*cos(om2*t)
    #tom1,tom2 = t*(om1+om2)/2,t*(om1-om2)/2
    #x3 = x10*cos(tom1)*cos(tom2) + x20*sin(tom1)*sin(tom2)
    #x4 = x10*sin(tom1)*sin(tom2) + x20*cos(tom1)*cos(tom2)
    lines([x1,x2])
end

# pend0

function pend0(;w0=1,th0=0,thmax=90,imax=10,tol=1e-5,N=25)
    a3 = w0^2/6
    dth = (thmax-th0)/N
    th = [th0+(j-1)*dth for j ∈ 1:N]
    A1 = zeros(N)
    for j ∈ 1:N
        x = th[j]-1 # initial guess
        xn = 999
        f = 999
        i = 0
        while (abs(xn-x) ≥tol) & (f != 0) & (i <imax)
            x = xn
            f = x + a3*x^3/(27a3*x^2-32w0^2)-th[j]
            fp = 1+3a3*x^2/(27a3*x^2-32w0^2) - 54a3^2*x^4/((27a3*x^2-32w0^2)^2)
            xn = x - f/fp
            i += 1
        end
        A1[j] = xn
    end
    TensorField(th,A1)
end

# pend1

function pend1(;thmax=90,N=100)
    dth = thmax/N
    th = TensorField(0:dth:thmax)
    m = sin(th*2pi/360/2)^2
    y1 = 1/sqrt(1-(th*2pi/360)^2/8)
    y2 = 2ellipke(m)/pi
    lines([y1,y2])
end

# pend2

function pend2(;w0=1)
    cf = 2pi/360
    τ0 = 2pi/w0
    tmax = 4τ0
    th = 90 # initial amplitude in degrees
    thr = th*cf # initial angle in radians
    ic1 = Chain(thr,0.0)
    ic = InitialCondition(Flow(pend2_der(w0),tmax),ic1)
    th2 = odesolve(ic,ExplicitIntegrator{4}(1e-4))
    om = w0*sqrt(1-thr^2/8)
    a3 = w0^2/6
    A1 = thr
    A3 = a3*A1^3/(27*a3*A1^2-32*w0^2)
    t = TensorField(points(th2))
    th1 = thr*cos(om*t) + A3*cos(3om*t)
    th0 = thr*cos(w0*t) # the SHO case
    lines([th0,th1,getindex.(th2,1)]/cf)
end

pend2_der(w0) = x -> Chain(x[2],-w0^2*sin(x[1]))

# molec

function molec(;tmax=2,xb=3/2,xi=0.2)
    ic1 = Chain(xb+xi,0.0)
    ic = InitialCondition(Flow(molec_der(),tmax),ic1)
    x2 = odesolve(ic,ExplicitIntegrator{4}(1e-4))
    t = TensorField(points(x2))
    A1 = (-1+sqrt(1+xi*32/9))*9/16
    x1 = xb+A1*cos(2pi*t)+4A1^2*(3-cos(4pi*t))/9
    x0 =xb+xi*cos(2pi*t)
    lines([x0,x1,getindex.(x2,1)])
end

molec_der() = x -> Chain(x[2],81pi^2*(3/(x[1]^4)-2/x[1]^3)/8)

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

# curl_ex

function curl_ex(n,j=1;vmin=-1,vmax=1)
    xmax,ymax,zmax = vmax,vmax,vmax
    xmin,ymin,zmin = vmin,vmin,vmin
    vs = 0.3
    xs,ys,zs = vs,vs,vs
    N = (vmax-vmin)/vs
    m = Int(round(N/2+1))
    zm = zmin+(m-1)*zs
    XYZ = TensorField(ProductSpace{S"3"}(xmin:xs:xmax,ymin:ys:ymax,zmin:zs:zmax))
    x,y,z = split(XYZ)
    f = if n == 1
        Chain.(-y,x,0)
    elseif n == 2
        Chain.(x*z,-y^2,2x^2*y)
    elseif n == 3
        Chain.(y^2,2*x*y+z^2,2y*z)
    elseif n == 4
        Chain.(x*y,y*z,z*x)
    elseif n == 5
        Chain.(x^2*y,y^2*z,z^2*x)
    end
    if j == 1
        streamplot(f)
    elseif j == 2
        Chain.(x[:,:,m],y[:,:,m])
    elseif j == 3
        cav2 = curl(Chain.(getindex.(f[:,:,m],1),getindex.(f[:,:,m],2)))
    elseif j == 4
        streamplot(curl(f))
    end
end

# Chapter 6

# parabola

function parabola(;x0=1,z0=1,xmin=-3,xmax=3,xs=0.1)
    x = TensorField(xmin+x0:xs:xmax+x0)
    zt,zb = 0,0
    out = nothing
    for e ∈ -0.5:-0.5:-1.5
        z = z0+(x-x0)^2/4/e
        e == -0.5 ? (out = lines(z)) : lines!(z)
    end
    return out
end

# projectile

function projectile()
end



