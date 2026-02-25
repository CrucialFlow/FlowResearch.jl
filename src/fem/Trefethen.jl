# https://github.com/chakravala

using Grassmann, Cartan
using Makie, GLMakie
using FFTW, ToeplitzMatrices, SparseArrays

#= Chapter 1 =#

# Program 1 # fourth order difference convergence

fun(N) = exp(sin(pts(N)))
pts(N=24,h=2π/N) = TensorField((h*(1:N)).-π)

function program1(N)
    D = spdiagm(1=>2ones(N-1)/3,2=>-ones(N-2)/12,N-1=>-2ones(1)/3,N-2=>ones(2)/12)
    (D-D')/(2π/N)
end

dprime(x) = program1(length(x))*x
lines(dprime(fun(240)))

error(N,grad) = grad(fun(N[1])) - cos(pts(N[1]))*fun(N[1])
lines(error(240,dprime))
lines(error(24,gradient))
lines!(error(24,gradient_rfft))

errorlog(x,grad) = TensorField(log10.(x),log10.(fiber.(maximum.(abs.(error.(x,grad))))))
scatter(errorlog((2).^(2:12),gradient))
scatter(errorlog((2).^(0:12),gradient_rfft))
scatter(errorlog((2).^(1:12),dprime))

# infinite-like operator

cols(N) = vcat(0,((-1).^(1:N-1))./(1:N-1))
infinite(N) = Toeplitz(cols(N),-cols(N))

# Program 2 # periodic spectral convergence

fun(N) = exp(sin(pts(N)))
pts(N=24,h=2π/N) = TensorField((h*(1:N)).-π)

dprime(N,fun) = gradient_toeplitz(fun(N))
lines(dprime(240,fun))

error(N,grad) = grad(fun(N[1])) - cos(pts(N[1]))*fun(N[1])
lines(error(240,gradient_toeplitz))

errorlog(x,grad) = TensorField(log10.(x),log10.(fiber.(maximum.(abs.(error.(x,grad))))))
scatter(errorlog(2:2:100,gradient_toeplitz))

#= Chapter 2 =#

# Program 3 # band-limited interpolation

t = TensorField(-11:1.0:11)

lines(resample_sinc(iszero(t),221))
scatter!(iszero(t))

lines(resample_sinc((x->abs(x[1])≤3).(t),221))
scatter!((x->abs(x[1])≤3).(t))

lines(resample_sinc(max(0,1-abs(t)/3)))
scatter!(max(0,1-abs(t)/3))

#= Chapter 3 =#

t = TensorField(range(-pi,3pi,2*24+1)) # coarse grid
tt = TensorField(range(-pi,3pi,1000)) # fine grid
lines(mysinc(tt)) # mysinc interpolation
scatter!(mysinc(t)) # periodic sinc hat
lines(mysinc(t)) # line interpolation
scatter!(mysinc(t)) # periodic sinc hat

# Program 4 # periodic spectral differentiation

dprime(N,fun) = derivetoeplitz(N)*fun(N)
pts(N=24,h=2π/N) = TensorField(h*(1:N))

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

#= Chapter 4 =#

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

#= Chapter 5 =#

# Program 9 # Chebyshev interpolation

x = TensorField(range(-1,1,17))
xx = TensorField(Chebyshev(17))

lines(resample_lagrange(1/(1+16x^2),401))
scatter!(1/(1+16x^2))

lines(resample_lagrange(1/(1+16xx^2),401))
scatter!(1/(1+16xx^2))

tt = TensorField(range(-1,1,401));
maximum(abs.(fiber(1/(1+16tt^2)) - fiber(resample_lagrange(1/(1+16xx^2),401))))

# Program 10

x = TensorField(LagrangeWeights(range(-1,1,17)))
xx = TensorField(Chebyshev(17))
Z = Complex(complexify(TensorField(ProductSpace{2}(-1.4:0.02:1.4,-1.12:0.02:1.12))))

lines(resample_roots(x,401))
scatter!(0x)

lines(resample_roots(xx,401))
scatter!(0xx)

contour(rootspolynomial(x,Z),levels=10.0.^(-4:0))
scatter!(x+0im)

contour(rootspolynomial(xx,Z),levels=10.0.^(-4:0))
scatter!(xx+0im)

#= Chapter 6 =#

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

#= Chapter 7 =#

# Program 13

x = TensorField(Chebyshev(17)) # domain
u = solvedirichlet(helmholtz(x),exp(4x)) # Dirichlet invert
lines(resample_lagrange(u,201))
scatter!(u)

# Program 14

u0 = TensorField(Chebyshev(17),zeros(17)) # initialize
u = solveiteration(helmholtz(u0),exp,u0,solvedirichlet)
lines(resample_lagrange(u,201))
scatter!(u)

# Program 15

function solve_eig(N=36)
    x = Chebyshev(N)
    D2 = (ChebyshevMatrix(x)^2)[2:N-1,2:N-1]
    U = eigvecs(D2)
    [TensorField(x,vcat(0,U[:,j],0)) for j ∈ 1:N-2]
end

function show_eig(x,n=201)
    out = lines(resample_lagrange(x,n))
    scatter!(x)
    return out
end

se = solve_eig(36)
show_eig(se[30],201)
show_eig(se[25],201)
show_eig(se[20],201)
show_eig(se[15],201)
show_eig(se[10],201)
show_eig(se[5],201)

# Program 16

fun(x) = 10sin(8x[1]*(x[2]-1))
x = Chebyshev(25)
XY = TensorField(ProductSpace{2}(x,x))
u = solvedirichlet(helmholtz(XY),fun.(XY))
linegraph(graph(u),u)
ur = resample_lagrange(u,51,51)
linegraph(graph(ur),ur)

# Program 17

x = Chebyshev(25)
XY = TensorField(ProductSpace{2}(x,x))

fun(x) = exp(-10*((x[2]-1)^2+(x[1]-0.5)^2))
u = solvedirichlet(helmholtz(XY,9),fun.(XY))
contour(u)
contour(resample_lagrange(u,61,61))

function solvehelmholtzwave(u0,k,t)
    data = zeros(Complex{Float64},size(u0)...,length(t))
    data[:,:,1] = fiber(u0)
    for i ∈ 2:length(t)
        data[:,:,i] = fiber(exp(k*im*fiber(t)[i])*u0)
    end
    TensorField(base(u0)⊕base(TensorField(t)),data)
end

hw = solvehelmholtzwave(hh,9,0:0.01:1)

#= Chapter 8 =#

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

# Program 18

xx = TensorField(-1.01:0.01:1.01)
lines(exp(xx)*sin(5xx))

function program18(N)
    x = TensorField(Chebyshev(N))
    error = gradient_chebyshevfft(exp(x)*sin(5x)) - exp(x)*(sin(5x)+5cos(5x))
end

lines(program18(10))
lines(program18(20))

x = TensorField(Chebyshev(20))
y = exp(x)*sin(5x)
lines(gradient_chebyshev(y))
lines(gradient_chebyshev(y) - gradient_chebyshevfft(y))
lines(gradient_chebyshevfft(y))

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
linegraph(resample_lagrange(lf[:,:,2],33,33))

x = Chebyshev(31)
XYZ = TensorField(ProductSpace{3}(x,x,x))
X,Y,Z = split(XYZ)
ex = exp(-40((X-0.4)^2+Y^2+Z^2))
lf = leapfrog(ex,41,4,1/30)
contour(lf[:,:,2],alpha=0.1)

#= Chapter 9 =#

# Program 21

x = TensorField((2π/42)*(1:42))
D2 = derivetoeplitz2(42)
evs(q) = eigvals(-D2 + 2q*Diagonal(fiber(cos(2x))))
out = hcat(evs.(0:0.2:15)...)
lines(out[1,:])
[lines!(out[i,:]) for i ∈ 2:11]

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
contour(resample_lagrange(reshapedirichlet(XY,V[:,1]),101,101))

y = x[2:end-1]
X,Y = split(TensorField(ProductSpace(y,y)))
perturb = Diagonal(vec(fiber(exp(20*(X-Y-1)))))

# Program 24

function daviesosc(c,L=6,N=71)
    x = Chebyshev(range(-L,L,N))
    D = ChebyshevMatrix(x)
    (-D^2)[2:end-1,2:end-1] + c*Diagonal(x[2:end-1].^2)
end
function pseudospectra(A,Z)
    In = I(size(A)[1])
    sigmin = zeros(size(Z)...)
    for j ∈ 1:size(Z)[2], i ∈ 1:size(Z)[1]
        sigmin[i,j] = minimum(svd(fiber(Z)[i,j]*In-A).S)
    end
    return TensorField(Z,sigmin)
end

A = daviesosc(1+3im,6,71)
Z = Complex(complexify(TensorField(ProductSpace(0:2:50,0:2:40))))
contour(pseudospectra(A,Z),levels=10.0.^(-4:0.5:-0.5))

#= Chapter 10 =#

# Program 25

t = TensorField(range(0,2pi,201))
z = exp(im*t)
r = z-1

# Adams-Bashforth
lines(r/1)
lines!(2r/(3-1/z))
lines!(12r/(23-16/z+5/z^2))

# Adams-Moulton
lines(12r/(5z+8-1/z))
lines!(24r/(9z+19-5/z+1/z^2))
lines!(720r/(251z+646-264/z+106/z^2-19/z^3))
dd = 1-1/z
lines!(dd/(1-dd/2-dd^2/12-dd^3/24-19dd^4/720-3dd^5/160))

# Backward differentiation
r = dd
lines(r)
for i ∈ 2:6
    r += (dd^i)/i
    lines!(r)
end

# Runge-Kutta

W = TensorField(z,zeros(Complex{Float64},length(z)))
for i ∈ 2:length(z)
    W[i] = W[i-1]-(1+W[i-1]-fiber(z)[i])
end
lines(W)
W = TensorField(z,zeros(Complex{Float64},length(z)))
for i ∈ 2:length(z)
    W[i] = W[i-1]-(1+W[i-1]+0.5W[i-1]^2-fiber(z)[i]^2)/(1+W[i-1])
end
lines!(W)
W = TensorField(z,zeros(Complex{Float64},length(z)))
for i ∈ 2:length(z)
    W[i] = W[i-1]-(1+W[i-1]+0.5W[i-1]^2+W[i-1]^3/6-fiber(z)[i]^3)/(1+W[i-1]+W[i-1]^2/2)
end
lines!(W)
W = TensorField(z,zeros(Complex{Float64},length(z)))
for i ∈ 2:length(z)
    W[i] = W[i-1]-(1+W[i-1]+0.5W[i-1]^2+W[i-1]^3/6+W[i-1]^4/24-fiber(z)[i]^4)/(1+W[i-1]+W[i-1]^2/2+W[i-1]^3/6)
end
lines!(W)

# Program 26

x = TensorField(Chebyshev(61))
D2 = (ChebyshevMatrix(x)^2)[2:end-1,2:end-1]
Lam,V = eigen(D2)
lines(TensorField(log10.(1:59),log10.(-reverse(Lam))))

v4N = TensorField(x,[0; V[:,end-Int(60/4)+2]; 0])
lines(resample_lagrange(-v4N,201))
scatter!(-v4N)

vN = TensorField(x,[0;V[:,2];0])
scatter(abs(vN))

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

#= Chapter 11 =#

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

#= Chapter 12 =#

# Program 30

fun1(x) = abs(integrate_chebyshev(abs(x)^3) - 0.5)
fun2(x) = abs(integrate_chebyshev(exp(-x^-2)) -2*(exp(-1)+sqrt(π)*(erf(1)-1)))
fun3(x) = abs(integrate_chebyshev(1/(1+x^2)) - π/2)
fun4(x) = abs(integrate_chebyshev(x^10) - 2/11)

myfun(fun,N=51) = TensorField(2:N,[fun(TensorField(Chebyshev(i))) for i ∈ 2:N])

lines(log(myfun(fun1)))
lines(log(myfun(fun2)))
lines(log(myfun(fun3)))
lines(bound(log(myfun(fun4)),50))

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

#= Chapter 13 =#

# Program 32

x = TensorField(Chebyshev(17)) # domain
u = solvedirichlet(helmholtz(x),exp(4x))+(x+1)/2 # u(1) = 1
lines(resample_lagrange(u,201))
scatter!(u)

# Program 33

function myoperator(x)
    D = ChebyshevMatrix(x)
    D2 = D^2
    D2[1,:] = D[1,:]
    D2[1:end-1,1:end-1]
end

leftneumann(L,u) = TensorField(u,vcat(L\vcat(0,fiber(u)[2:end-1]),0))

x = TensorField(Chebyshev(17)) # domain
u = leftneumann(myoperator(x),exp(4x))
lines(resample_lagrange(u,201))
scatter!(u)
exact = (exp(4x)-4exp(-4)*(x-1) - exp(4))/16
maximum(abs(u-exact))

# Program 34

function allencahn(u0,eps=0.01,tplot=2,tmax=50,n=1)
    x = TensorField(points(u0))
    dt = min(0.01,50length(u0)^-4/eps)
    D2 = ChebyshevMatrix(u0)^2
    D2[[1,end],:] .= 0
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

# Program 35

function acbc(x)
    fiber(x)[1] = -1
    fiber(x)[end] = 1 + sin(point(x)/5)^2
    return x
end
function allencahn(u0,eps=0.01,tplot=2,tmax=50,n=1)
    x = TensorField(points(u0))
    dt = min(0.01,50length(u0)^-4/eps)
    D2 = ChebyshevMatrix(u0)^2
    plotgap = round(tplot/dt)
    allencahn(v) = (eps*D2)*(v-x) + v - v^3
    ic = InitialCondition(allencahn,u0,tmax)
    out = odesolve(ic,ExplicitIntegrator{n}(tplot/plotgap,plotgap),acbc)
    x,t = split(points(out))
    TensorField(ProductSpace(x,t/10),fiber(out)/10)
end

x = TensorField(Chebyshev(21))
u0 = 0.53x + 0.47sin(-1.5pi*x)
surface(bound(allencahn(u0),0.2))

# Program 36

x = Chebyshev(25)
XY = TensorField(ProductSpace{2}(x,x))
X,Y = split(XY)
D2 = ChebyshevMatrix(x)^2
L = kron(I(25),D2) + kron(D2,I(25))
b = findall(x -> isone(abs(x[1])) | isone(abs(x[2])),vec(fiber(XY)))
L[b,:] .= 0
L[b,b] .= I(length(b))
bx = findall(x -> isone(x[1]),vec(fiber(XY)))
by = findall(x -> isone(x[2]) & (x[1]<0),vec(fiber(XY)))
rhs = zeros(25^2)
rhs[bx] = 0.2.*sin.(3π.*fiber(Y)[bx])
rhs[by] = sin.(π.*fiber(X)[by]).^4
u = TensorField(XY,reshape(L\vec(fiber(rhs)),size(XY)))
surface(u)
u(0,0)
wireframe(resample_lagrange(u,51,51))

# Program 37

x = TensorField(range(-3,3,50))
D2x = (π/3)^2*derivetoeplitz2(50)
y = Chebyshev(16)
Dy = ChebyshevMatrix(y)
D2y = Dy^2
BC = Dy[[1,16],[1,16]]\Dy[[1,16],2:15]
XY = TensorField(ProductSpace{2}(points(x),y))
X,Y = split(XY)
dt = 5/(50+15^2)
plotgap = Int(round(2/dt))
dt = 2/plotgap
vv = exp(-8((X+1.5)^2+Y^2))
vvold = exp(-8((X+dt+1.5)^2+Y^2))

fprime(v) = TensorField(base(fiber(v)),D2x*fiber(fiber(v)) + fiber(fiber(v))*D2y)
function wavebc(v)
    fiber(fiber(v))[:,[1,16]] .= fiber(fiber(v))[:,2:end-1]*transpose(BC)
    return v
end

function leapfrog(vv,vvold,tmax=4)
    lc = LeapCondition(fprime,vvold,vv,2/plotgap,tmax)
    odesolve(lc,LeapIntegrator{2}(plotgap),wavebc)
end

lf = leapfrog(vv,vvold,4)
surface(lf[:,:,1])

#= Chapter 14 =#

# Program 38

x = TensorField(Chebyshev(24))
u = solvedirichlet(biharmonic(x),exp(-0.1x^2))
lines(resample_lagrange(u,201))
scatter!(u)

# Program 39

x = Chebyshev(18)
XY = TensorField(ProductSpace{2}(x,x))
B = eigvecs(biharmonic(XY))
contour(reshapedirichlet(XY,B[:,1]))
contour(resample_lagrange(reshapedirichlet(XY,B[:,1]),201,201))

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

