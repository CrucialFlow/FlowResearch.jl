# github.com/chakravala
# Phys: Classical Mechanics, Hasbun

# Chapter 1

# ho1.m

function init_mass_spring(;tmax=1,k=1000,m=5,C=0,E=0,ω=0,x0=0.1,v0=0)
    k # 1000 # spring constant
    m # 5.0 # bob mass
    C # 0.0 # damping coefficient
    E # 0.0 # driving force magnitude
    ω # 0.0 # driving force frequency
    x0 # 0.1 # initial position
    v0 # 0.0 # initial velocity
    t0 = 0.0 # initial time
    F0 = -k*x0-C*v0+E*sin(ω*t0) # initial force
    a0 = F0/m # initial acceleration
    InitialCondition(mass_spring(k,C,F0,ω,m),t0↦Chain(x0,v0),tmax)
end

mass_spring(k,C,F0,ω,m) = x -> Chain(x[2],(F0*sin(ω*point(x))-k*x[1]-C*x[2])/m)
function solve_mass_spring(;n=100,tmax=1,k=1000,m=5,C=0,E=0,ω=0,x0=0.1,v0=0)
    ic = init_mass_spring(;tmax,k,m,C,E,ω,x0,v0)
    getindex.(odesolve(ic,ExplicitIntegrator{4}(tmax/n)),1)
end

solve_mass_spring(n=100,tmax=1,k=1000,m=5,C=0,E=0,ω=0,x0=0.1,v0=0)
solve_mass_spring(n=100,tmax=10,k=1,m=1,C=0,E=0,ω=0,x0=1,v0=0)
solve_mass_spring(n=200,tmax=20,k=1,m=1,C=0.5,E=0,ω=0,x0=1,v0=0)
solve_mass_spring(n=200,tmax=20,k=1,m=1,C=0.5,E=0.1,ω=0.8,x0=1,v0=0)

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
