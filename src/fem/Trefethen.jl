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

# Program 4 # periodic spectral differentiation

pts(N=24,h=2π/N) = TensorField(h*(1:N))

cols(N,h=2π/N) = vcat(0,0.5*(-1).^(1:N-1).*cot.((1:N-1)*h/2))
spectral(N) = Toeplitz(cols(N),-cols(N))
dprime(N,fun) = spectral(N)*fun(N)

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
fun(N) = exp(-100*(pts(N)-1)^2)

function leapfrog(fun;N=128,tmax=8,tplot=0.15,t=0)
    x = pts(N)
    dt = (2π/N)/4
    c = 0.2 + sin(x-1)^2
    v = fun(x)
    vold = fun(x-0.2dt)
    plotgap = Int(round(tplot/dt))
    nplots = Int(round(tmax/tplot))
    dt = tplot/plotgap
    data = zeros(N,nplots+1)
    data[:,1] = fiber(v)
    for i ∈ 1:nplots
        for n ∈ 1:plotgap
            vold,v = v,vold-(2dt*c)*gradient_fft(v)
        end
        data[:,i+1] = fiber(v)
    end
    TensorField(base(v)⊕(0:dt*plotgap:tmax),data)
end

alteration(leapfrog(fun),0.01,lines,lines!)
linegraph(leapfrog(fun))

function my_leapfrog(fun;N=128,tmax=8,tplot=0.15,t=0,o=3)
    x = pts(N)
    dt = (2π/N)/4
    c = 0.2 + sin(x-1)^2
    plotgap = Int(round(tplot/dt))
    f(x) = (-c)*gradient_fft(localfiber(x))
    ic = InitialCondition(f,t↦fun(N),tmax)
    #odesolve(ic,MultistepIntegrator{o}(tplot/plotgap))
    odesolve(ic,ExplicitIntegrator{o}(tplot/plotgap,plotgap))
end

wireframe(graph(my_leapfrog(fun)))

# Project 7 # accuracy of periodic spectral differentiation

pts(N=128,h=2π/N) = TensorField(h*(1:N))
cols(N,h=2π/N) = vcat(0,0.5*(-1).^(1:N-1).*cot.((1:N-1)*h/2))
spectral(N) = Toeplitz(cols(N),-cols(N))
dprime(N,fun::Function) = spectral(N)*fun(N)
dprime(N,fun) = spectral(N)*fun
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

cols(N,h=2π/N) = vcat(-π^2/(3*h^2)-1/6,-0.5*(-1).^(1:N-1)./sin.((1:N-1)*h/2).^2)
spectral(N) = Toeplitz(cols(N),cols(N))

function program8(N;L=8)
    x = TensorField(fiber((L/π)*(pts(N)-π)))
    D2 = (π/L)^2*spectral(N)
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

chebyshev(N) = TensorField(cos.(π*(0:N)/N))

function cheb(N::Int) # differentiation matrix
    iszero(N) && (return [0;;])
    x = chebyshev(N) # grid
    c = vcat(2,ones(N-1),2).*(-1).^(0:N)
    X = repeat(fiber(x),1,N+1)
    D = (c*inv.(c)')./((X-X')+I) # off-diagonal entries
    D-Diagonal(vec(sum(D,dims=2))) # diagonal entries
end

cheb(1)
cheb(2)
cheb(3)
cheb(4)
cheb(5)

# Program 11 # Chebyshev differentiation of smooth function

xx = TensorField(-1.01:0.01:1.01)
lines(exp(xx)*sin(5xx))

function program11(N)
    x = chebyshev(N)
    error = cheb(N)*(exp(x)*sin(5x)) - exp(x)*(sin(5x)+5cos(5x))
end

lines(program11(10))
lines(program11(20))

# Program 12 # accuracy of Chebyshev spectral differentiation (compare Program 7)

pts(N=128,h=2π/N) = TensorField(h*(1:N))
dprime(N,fun::Function) = cheb(N)*fun(N)
dprime(N,fun) = cheb(N)*fun
error(N,v,vprime) = fiber(maximum(norm(dprime(N,v)-vprime)))

function fun(N,o)
    h = 2π/N
    x = chebyshev(N)
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

