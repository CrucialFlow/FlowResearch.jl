# https://github.com/chakravala

using Grassmann, Cartan
using Makie, GLMakie
using FFTW, ToeplitzMatrices, SparseArrays

# Program 1 # fourth order difference convergence

fun(N) = exp(sin(pts(N)))
pts(N=24,h=2π/N) = TensorField((h*(1:N)).-π)

function program1(N)
    D = spdiagm(1=>2ones(N-1)/3,2=>-ones(N-2)/12,N-1=>-2ones(1)/3,N-2=>ones(2)/12)
    (D-D')/(2π/N)
end

dprime(N) = program1(N)*fun(N)
lines(dprime(240))
lines(dprime(240)-cos(pts(240))*fun(240))

error(N) = fiber(maximum(abs(dprime(N)-cos(pts(N))*fun(N))))
errorlog(x) = TensorField(log10.(x),log10.(error.(x)))
scatter(errorlog((2).^(3:12)))


fun(N) = exp(sin(pts(N)))
pts(N=24,h=2π/N) = TensorField((h*(1:N)).-π)
error(N) = abs(gradient(fun(N))-cos(pts(N))*fun(N))
error_fft(N) = abs(gradient(fun(N))-cos(pts(N))*fun(N))
errorlog(x) = TensorField(log10.(x),log10.(fiber.(maximum.(error.(x)))))
errorlog_fft(x) = TensorField(log10.(x),log10.(fiber.(maximum.(error_fft.(x)))))

scatter(errorlog((2).^(3:12)))

fun(N) = exp(sin(pts(N)))
pts(N=24,h=2pi/N) = TensorField((h*(1:N)).-pi)
error(N) = abs(gradient(fun(N))-cos(pts(N))*fun(N))
error_fft(N) = abs(gradient_fft(fun(N))-cos(pts(N))*fun(N))
errorlog(x,y) = TensorField(log10.(x),log10.(fiber.(maximum.(y))))
graylines(errorlog((2).^(2:12),error.((2).^(2:12))))
graylines!(errorlog((2).^(0:12),error_fft.((2).^(0:12))))


# infinite-like operator

cols(N) = vcat(0,((-1).^(1:N-1))./(1:N-1))
infinite(N) = Toeplitz(cols(N),-cols(N))

# Program 2 # periodic spectral convergence

fun(N) = exp(sin(pts(N)))
pts(N=24,h=2π/N) = TensorField((h*(1:N)).-π)

cols(N,h=2π/N) = vcat(0,0.5*(-1).^(1:N-1).*cot.((1:N-1)*h/2))
spectral(N) = Toeplitz(cols(N),-cols(N))

dprime(N,fun) = spectral(N)*fun(N)
lines(dprime(240,fun))
lines(dprime(240,fun)-cos(pts(240))*fun(240))

error(N) = spectral(N,h)*exp(sin(pts(N))) - cos(pts(N))*exp(sin(pts(N)))
errorlog(x) = TensorField(log10.(x),log10.(error.(maximum.(x))))
scatter(errorlog(2:2:100))


error(N) = fiber(maximum(abs(dprime(N,fun)-cos(pts(N))*fun(N))))
errorlog(x) = TensorField(log10.(x),log10.(error.(x)))
scatter(errorlog(2:2:100))

# Program 3 # band-limited interpolation

function program3(fun,h=1,xmax=10)
    x = -xmax:h:xmax # computational grid
    xx = TensorField(-xmax-h/20:h/10:xmax+h/20) # plot grid
    v = fun.(x)
    p = zeros(length(xx))
    for i ∈ 1:length(x)
        p += v[i]*fiber(sin(π*(xx-x[i])/h)/(π*(xx-x[i])/h))
    end
    return TensorField(xx,p)
end

program3(x->iszero(x[1]))
program3(x->abs(x[1])≤3)
program3(x->max(0,1-abs(x[1])/3))

lines(program3(x->iszero(x[1])))
scatter!(resample(program3(x->iszero(x[1])),21)) # n = 2xmax/h+1

lines(program3(x->abs(x[1])≤3))
scatter!(resample(program3(x->abs(x[1])≤3),21)) # n = 2xmax/h+1

lines(program3(x->max(0,1-abs(x[1])/3)))
scatter!(resample(program3(x->max(0,1-abs(x[1])/3)),21)) # n = 2xmax/h+1

t = TensorField(-10:1:10)
tt = TensorField(-10.5:0.01:10.5)
lines(sin(π*tt)/(π*tt))
scatter!(spike(t))


# Program 4 # periodic spectral differentiation

pts(N=24,h=2π/N) = TensorField(h*(1:N))

cols(N,h=2π/N) = vcat(0,0.5*(-1).^(1:N-1).*cot.((1:N-1)*h/2))
spectral(N) = Toeplitz(cols(N),-cols(N))
dprime(N,fun) = derivetoeplitz(N)*fun(N)

hat(N) = max(0,1-abs(pts(N)-π)/2)

lines(hat(24))
scatter!(hat(24))
lines(dprime(24,hat))
scatter!(dprime(24,hat))

fun(N) = exp(sin(pts(N)))

lines(fun(24))
scatter!(fun(24))
lines(dprime(24,fun))
scatter!(dprime(24,fun))

# Program 5

pts(N=24,h=2π/N) = TensorField(h*(1:N))
hat(N) = max(0,1-abs(pts(N)-π)/2)

diffop(N::Int) = im*vcat(0:N/2-1,0,-N/2+1:-1)
function program5(v)
    v_hat = fft(v)
    w_hat = diffop(length(v)).*v_hat
    w = real(ifft(w_hat))
end

lines(hat(24))
scatter!(hat(24))

lines(gradient_fft(hat(24)))
scatter!(gradient_fft(hat(24)))

fun(N) = exp(sin(pts(N)))

lines(fun(24))
scatter!(fun(24))

lines(gradient_fft(fun(24)))
scatter!(gradient_fft(fun(24)))

# Program 6 # variable coefficient wave equation

pts(N=128,h=2π/N) = TensorField(h*(1:N))


function leapfrog(fun,N=128,tmax=8,tplot=0.15,t=0)
    h = 2pi/N
    dt = h/4
    x = TensorField(h*(1:N))
    c = 0.2 + sin(x-1)^2
    fprime(v) = c*gradient_rfft(localfiber(v))
    vold = LocalTensor(t-dt,fun(x-0.2dt))
    v = LocalTensor(t,fun(x))
    lc = LeapCondition(fprime,vold,v,tmax)
    odesolve(lc,LeapIntegrator{1}(Int(round(tplot/dt))))
end

fun(x) = exp(-100*(x-1)^2)
linegraph(leapfrog(fun,128,8,0.15))

alteration(leapfrog(fun),0.01,lines,lines!)
linegraph(leapfrog(fun))



function my_leapfrog(fun;N=128,tmax=8,tplot=0.15,t=0,o=3)
    x = pts(N)
    dt = (2π/N)/4
    c = 0.2 + sin(x-1)^2
    plotgap = Int(round(tplot/dt))
    f(x) = (-c)*gradient_fft(localfiber(x))
    ic = InitialCondition(f,t↦fun(x),tmax)
    #odesolve(ic,MultistepIntegrator{o}(tplot/plotgap))
    odesolve(ic,ExplicitIntegrator{o}(tplot/plotgap,plotgap))
end

wireframe(graph(my_leapfrog(fun)))

# Project 7 # accuracy of periodic spectral differentiation

pts(N=128,h=2π/N) = TensorField(h*(1:N))
cols(N,h=2π/N) = vcat(0,0.5*(-1).^(1:N-1).*cot.((1:N-1)*h/2))
spectral(N) = Toeplitz(cols(N),-cols(N))
dprime(N,fun::Function) = derivetoeplitz(N)*fun(N)
dprime(N,fun) = derivetoeplitz(N)*fun
error(N,v,vprime) = fiber(maximum(norm(dprime(N,v)-vprime)))

function fun(N,o)
    h = 2π/N
    x = pts(N)
    if o == 1
        v = abs(sin(x))^3
        vprime = 3sin(x)*cos(x)*abs(sin(x))
        error(N,v,vprime)
    elseif o == 2
        v = exp(-sin(x/2)^-2)
        vprime = 0.5v*sin(x)/sin(x/2)^4
        error(N,v,vprime)
    elseif o == 3
        v = 1/(1+sin(x/2)^2)
        vprime = -sin(x/2)*cos(x/2)*v^2
        error(N,v,vprime)
    elseif o == 4
        v = sin(10x)
        vprime = 10cos(10x)
        error(N,v,vprime)
    end
end

scatter(TensorField(6:2:50,log.(fun.(6:2:50,1))))
scatter(TensorField(6:2:50,log.(fun.(6:2:50,2))))
scatter(TensorField(6:2:50,log.(fun.(6:2:50,3))))
scatter(TensorField(6:2:50,log.(fun.(6:2:50,4))))

# Program 8 # eigenvalues of harmonic oscillator

pts(N=128,h=2π/N) = TensorField(h*(1:N))

cols(N,h=2π/N) = vcat(-π^2/(3*h^2)-1/6,-0.5*(-1).^(2:N)./sin.((1:N-1)*h/2).^2)
spectral(N) = Toeplitz(cols(N),cols(N))

function program8(N;L=8)
    x = TensorField(fiber((L/π)*(pts(N)-π)))
    D2 = (π/L)^2*derivetoeplitz2(N)
    eigvals(Diagonal(fiber(x^2))-D2)
end

program8(6)
program8(12)
program8(18)
program8(24)
program8(30)
program8(36)

# Program 9 # Chebyshev interpolation

equispace(N) = TensorField(2(0:N)/N.-1)
chebyshev(N) = TensorField(cos.(π*(0:N)/N))

function program9(x,N=16)
    xx = TensorField(-1.01:0.005:1.01)
    u = 1/(1+16x^2)
    uu = 1/(1+16xx^2)
    # todo
end

# cheb.m

ChebyshevMatrix(1)
ChebyshevMatrix(2)
ChebyshevMatrix(3)
ChebyshevMatrix(4)
ChebyshevMatrix(5)

# Program 11 # Chebyshev differentiation of smooth function

xx = TensorField(-1.01:0.01:1.01)
lines(exp(xx)*sin(5xx))

function program11(N)
    x = Chebyshev(N)
    error = ChebyshevMatrix(x)*(exp(x)*sin(5x)) - exp(x)*(sin(5x)+5cos(5x))
end

lines(program11(10))
lines(program11(20))

x = TensorField(Chebyshev(20))
y = exp(x)*sin(5x)
lines(gradient_chebyshev(y))
lines(gradient_chebyshev(y) - gradient_chebfft(y))
lines(gradient_chebyshevfft(y))

# Program 12 # accuracy of Chebyshev spectral differentiation (compare Program 7)

pts(N=128,h=2π/N) = TensorField(h*(1:N))
dprime(N,fun::Function) = ChebyshevMatrix(N)*fun(N)
dprime(N,fun) = ChebyshevMatrix(N)*fun
error(N,v,vprime) = fiber(maximum(norm(dprime(N,v)-vprime)))

function fun(N,o)
    h = 2π/N
    x = Chebyshev(N+1)
    if o == 1
        v = abs(x)^3
        vprime = 3x*abs(x)
        error(N,v,vprime)
    elseif o == 2
        v = exp(-x^-2)
        vprime = 2v/x^3
        error(N,v,vprime)
    elseif o == 3
        v = 1/(1+x^2)
        vprime = -2x*v^2
        error(N,v,vprime)
    elseif o == 4
        v = x^10
        vprime = 10x^9
        error(N,v,vprime)
    end
end

scatter(TensorField(1:50,log.(fun.(1:50,1))))
scatter(TensorField(1:50,log.(fun.(1:50,2))))
scatter(TensorField(1:50,log.(fun.(1:50,3))))
scatter(TensorField(1:50,log.(fun.(1:50,4))))

# Program 13

x = TensorField(Chebyshev(17)) # domain
lines(solvedirichlet(helmholtz(x),exp(4x))) # Dirichlet invert

# Program 14

u0 = TensorField(Chebyshev(17),zeros(17)) # initialize
lines(solveiteration(helmholtz(u0),exp,u0,solvedirichlet))

# Program 15

function solve_eig(N=36)
    x = Chebyshev(N)
    D2 = (ChebyshevMatrix(x)^2)[2:N-1,2:N-1]
    U = eigvecs(D2)
    [TensorField(x,vcat(0,U[:,j],0)) for j ∈ 1:N-2]
end

lines(se[30])

# Program 16

fun(x) = 10sin(8x[1]*(x[2]-1))
x = Chebyshev(50)
XY = TensorField(ProductSpace{2}(x,x))
u = solvedirichlet(helmholtz(XY),fun.(XY))
linegraph(graph(u),u)

# Program 17

x = Chebyshev(50)
XY = TensorField(ProductSpace{2}(x,x))

fun(x) = exp(-10*((x[2]-1)^2+(x[1]-0.5)^2))
contour(solvedirichlet(helmholtz(XY,9),fun.(XY)))

function solvehelmholtzwave(u0,k,t)
    data = zeros(Complex{Float64},size(u0)...,length(t))
    data[:,:,1] = fiber(u0)
    for i ∈ 2:length(t)
        data[:,:,i] = fiber(exp(k*im*fiber(t)[i])*u0)
    end
    TensorField(base(u0)⊕base(TensorField(t)),data)
end

hw = solvehelmholtzwave(hh,9,0:0.01:1)

# chebfft

function gradient_chebfft(v)
    N = length(v)-1
    ii = 0:N-1
    U = real.(fft(vcat(fiber(v),reverse(fiber(v)[2:N]))))
    N = length(U)÷2
    ii = 0:N-1
    W = real.(ifft((im*vcat(ii,0,1-N:-1)).*U))
    w = zeros(N+1)
    w[2:N] = -W[2:N]./sqrt.(1.0.-cos.(π*(1:N-1)/N).^2)
    w[1] = sum((ii.^2).*U[ii.+1])/N .+ (0.5N)*U[N+1]
    w[N+1] = sum((-1).^(ii.+1).*(ii.^2).*U[ii.+1])/N .+ 0.5N*(-1)^(N+1)*U[N+1]
    return TensorField(v,w)
end

function gradient_chebrfft(v)
    N = length(v)-1
    x = points(v)
    iszero(N) && (return 0)
    ii = 0:N-1
    V = vcat(fiber(v),reverse(fiber(v)[2:N])) # transform x -> θ
    U = real.(rfft(V))
    W = real.(irfft((im*vcat(ii,N)).*U,2N+1))
    w = zeros(N+1)
    w[2:N] = -W[2:N]./sqrt.(1.0.-x[2:N].^2)
    w[1] = sum((ii.^2).*U[ii.+1])/N .+ (0.5N)*U[N+1]
    w[N+1] = sum((-1).^(ii.+1).*(ii.^2).*U[ii.+1])/N .+ 0.5N*(-1)^(N+1)*U[N+1]
    return TensorField(x,w)
end

# Program 19

fprime(v) = dirichlet!(gradient_chebyshevfft(gradient_chebyshevfft(v)))

function leapfrog(u0,N=81,tmax=4,tplot=0.075)
    dt = 8/(N-1)^2
    x = TensorField(Chebyshev(N))
    lc = LeapCondition(fprime,u0(x-dt),u0(x),dt,tmax)
    odesolve(lc,LeapIntegrator{2}(Int(round(tplot/dt))))
end

fun(x) = exp(-200x^2)
lf = leapfrog(fun,81,4,0.075)
linegraph(lf)

# Program 20


function fprime(v)
    uxx = gradient_chebyshevfft(gradient_chebyshevfft(v,1),1)
    uyy = gradient_chebyshevfft(gradient_chebyshevfft(v,2),2)
    dirichlet!(uxx + uyy)
end

fprime(v) = dirichlet!(laplacian_chebyshevfft(v))

function leapfrog(u0,N=25,tmax=1,tplot=1/3)
    dt = 6/(N-1)^2
    lc = LeapCondition(fprime,u0,u0,dt,tmax)
    odesolve(lc,LeapIntegrator{2}(Int(round(tplot/dt))))
end

x = Chebyshev(41)
XY = TensorField(ProductSpace{2}(x,x))
X,Y = split(XY)
ex = exp(-40((X-0.4)^2+Y^2))
lf = leapfrog(ex,41,4,1/30)
linegraph(lf[:,:,2])


x = Chebyshev(31)
XYZ = TensorField(ProductSpace{3}(x,x,x))
X,Y,Z = split(XYZ)
ex = exp(-40((X-0.4)^2+Y^2+Z^2))
lf = leapfrog(ex,41,4,1/30)
contour(lf[:,:,2],alpha=0.1)



# Program 22

function eigairy(N)
    x = Chebyshev(N)
    D2 = (ChebyshevMatrix(x)^2)[2:N-1,2:N-1]
    yi,xi = geneigsolve((D2,Diagonal(x[2:N-1])),10,:SR)
end

# Program 23

x = Chebyshev(17)
XY = TensorField(ProductSpace{2}(x,x))
L = helmholtz(XY)
V = eigvecs(-L)
contour(reshapedirichlet(XY,V[:,1]))

y = x[2:end-1]
X,Y = split(TensorField(ProductSpace(y,y)))
perturb = Diagonal(vec(fiber(exp(20*(X-Y-1)))))


# Program 27

function soliton(N=256,A=25,B=16)
    x = TensorField((2pi/N)*(-N/2:N/2-1))
    3A^2*sech(0.5(A*(x+2))).^2 + 3B^2*sech(0.5*(B*(x+1))).^2
end

function solvekdv(u0,tmax=0.006,o=4)
    N = length(u0)
    dt = 0.4/N^2
    tplot = tmax/25
    plotgap = Int(floor(tplot/dt))
    nplots = Int(round(tmax/tplot))
    v = rfft(u0)
    u = irfft(v)
    k = 0:length(v)-1#vcat(0:N/2-1,0,-N/2+1:-1)
    ik2 = (im/2)*k
    ik3 = im*k.^3
    g = dt.*ik2
    E = exp.(dt.*ik3./2)
    E2 = E.^2
    data = zeros(N,nplots+1)
    data[:,1] = fiber(u0)
    for i ∈ 1:nplots
        for n ∈ 1:plotgap
            a = g.*rfft(irfft(v).^2)
            b = g.*rfft(irfft(E.*(v.+a./2)).^2)
            c = g.*rfft(irfft(E.*v.+b./2).^2)
            d = g.*rfft(irfft(E2.*v.+E.*c).^2)
            v = E2.*v.+(E2.*a.+(2).*E.*(b.+c).+d)./6
            #u = irfft(v)
        end
        data[:,i+1] = fiber(irfft(v))
    end
    TensorField(base(u0)⊕(0:dt*plotgap:tmax),data)
end

function solvekdv(u0,tmax=0.006,o=4)
    N = Cartan.invdim(points(u0))
    dt = 0.4/N^2
    tplot = tmax/25
    plotgap = Int(floor(tplot/dt))
    k = 0:length(u0)-1#vcat(0:N/2-1,0,-N/2+1:-1)
    ik2 = (im/2)*k
    ik3 = im*k.^3
    function kdv(x)
        eik = exp.(point(x)*ik3)
        (ik2.*conj.(eik)).*rfft(irfft(eik.*localfiber(x))^2)
    end
    ic = InitialCondition(kdv,u0,tmax)
    odesolve(ic,ExplicitIntegrator{o}(tplot/plotgap,plotgap))
end

# Program 28


toeplitz2(N,h=2π/N) = vcat(-π^2/(3*h^2)-1/6,0.5*(-1).^(2:N)./sin.((1:N-1)*h/2).^2)
spectral(N,h=2pi/N,c=toeplitz2(N,h)) = Toeplitz(c,c)


function polarlaplacian(N=13,M=21)
    r = Chebyshev(2N)
    D = ChebyshevMatrix(r); D2 = D^2
    D1,D2 = D2[2:N,2:N],D2[2:N,2N-1:-1:N+1]
    E1,E2 = D[2:N,2:N],D[2:N,2N-1:-1:N+1]
    D2t,R = derivetoeplitz2(M-1),Diagonal(inv.(r[2:N]))
    M2 = Int((M-1)/2); Z = zeros(M2,M2)
    kron(I(M-1),D1+R*E1)+kron([Z I;I Z],D2+R*E2)+kron(D2t,R^2)
end


function reshapedisk(Vi,N,M)
    t = TorusParameter(M)
    r = TensorField(Chebyshev(2N)[N+1:2N])
    disk = unitdisk(t,r)
    x,y = fiber.(split(disk))
    u = reshapepolardirichlet(points(disk),real.(Vi))
    TensorField(disk,Chain.(x,y,u/norm(u,Inf)))
end

V = eigvecs(-polarlaplacian(13,21))
third(x) = getindex.(x,3)
mesh(reshapedisk(V[:,1],13,21),third)

# Program 29

diskparam = base(TorusParameter(41))⊕base(TensorField(Chebyshev(32)[17:end]))
fun(x) = -x[2]^2*sin(x[1]/2)^4 + sin(6x[1])*cos(x[1]/2)^2
pd = solvepolardirichlet(L,fun.(diskparam))


# clencurt

function clenshawcurtis(N::Int)
    θ = π*(0:N)/N
    #x = cos.(θ)
    w = zeros(N)
    ii = 2:N
    v = ones(N-1)
    if iszero(N%2)
        w[1] = 1/(N^2-1) #w[N+1] = w[1]
        for k ∈ 1:Int(N/2)-1
            v .-= 2cos.(2k*θ[ii])/(4k^2-1)
        end
        v .-= cos.(N*θ[ii])/(N^2-1)
    else
        w[1] = 1/N^2 # w[N+1] = 1/N^2
        for k ∈ 1:Int((N-1)/2)
            v .-= 2cos.(2k*θ[ii])/(4k^2-1)
        end
    end
    w[ii] .= (2/N).*v
    return w
end

# Gauss

function legendre(N)
    β = 0.5./sqrt.(1.0.-(2*(1:N-1)).^-2)
    T = spdiagm(1=>β,-1=>β)
    eigvals(Matrix(T))
end

function gausslegendre(N)
    β = 0.5./sqrt.(1.0.-(2*(1:N-1)).^-2)
    T = spdiagm(1=>β,-1=>β)
    2eigvecs(Matrix(T))[1,:].^2
end

# Program 31

function integrate_haar(f,z,θ=range(-pi,pi,70))
    N = length(θ)
    out = 0z
    for i ∈ 1:N
        out .+= f(fiber(θ)[i],z)
    end
    out/N
end

function integrate_haar(f,z,θ=range(-pi,pi,70))
    N = length(θ)
    out = zeros(typeof(z),N)
    out[1] = f(fiber(θ)[1],z)
    for i ∈ 2:N
        out[i] = out[i-1]+f(fiber(θ)[i],z)
    end
    TensorField(θ,out/N)
end

XY = TensorField(ProductSpace{2}(-3.5:0.1:4,-2.5:0.1:2.5))
zz = Complex(complexify(XY))
function gamcirc(c=-11,r=16)
    function fun(θ,z)
        t = c + r*exp(im*θ)
        (exp(t)*(t-c))*t^-z
    end
end

gam = inv(integrate_haar(gamcirc(-11,16),zz))
surface(bound(abs(gam),5))

# Program 33

function allencahn(u0,eps=0.01,tplot=2,tmax=50,n=1)
    x = TensorField(points(u0))
    dt = min(0.01,50length(u0)^-4/eps)
    D2 = ChebyshevMatrix(u0)^2
    plotgap = round(tplot/dt)
    allencahn(v) = (eps*D2)*(v-x) + v - v^3
    ic = InitialCondition(allencahn,u0,tmax)
    out = odesolve(ic,ExplicitIntegrator{n}(tplot/plotgap,plotgap))
    x,t = split(points(out))
    TensorField(ProductSpace(x,t/10),fiber(out)/10)
end

x = TensorField(Chebyshev(21))
u0 = 0.53x + 0.47sin(-1.5pi*x)
surface(bound(allencahn(u0),0.2))

# Program 38

x = TensorField(Chebyshev(24))
lines(solvedirichlet(biharmonic(x),exp(-0.1x^2)))

# Program 39


x = Chebyshev(18)
XY = TensorField(ProductSpace{2}(x,x))
B = eigvecs(biharmonic(XY))
contour(reshapedirichlet(XY,B[:,1]))

# Program 40

function orrsommerfeld(v,R=5772)
    N = length(v)
    x = points(v)
    D = ChebyshevMatrix(x)
    D2 = (D^2)[2:N-1,2:N-1]
    D4 = biharmonic(v,D)
    A = (D4-2D2+I(N-2))/R - 2im*I(N-2) - im*Diagonal(1.0.-x[2:N-1].^2)*(D2-I(N-2))
    B = D2 - I(N-2)
    return A,B
end


