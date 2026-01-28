# Fourier

function FourierSine(N,f)
    x = TensorField(base(f))
    L = fiber(x[end]-x[1]) # T/2
    ω = Cartan.FourierSpace((π/L)*(1:N),points(x))
    TensorField(ω,[(2/L)*integrate(f*sin(nω*x)) for nω ∈ ω])
end

x = TensorField(0:0.01:5)
b = FourierSine(8,1+x^2)
save("fouriersine1.pdf",graylines(b,3,rasterize=2))

fsum(b) = sum([fiber(b)*sin(x*point(b)) for b in b])
out = graylines(1+x^2,3,rasterize=2)
graylines!(fsum(b),3,rasterize=2)
save("fouriersine2.pdf",out)

lines(1+x^2)
for i ∈ 1:3
    lines!(sum(b.*sin.(((π/5)*(1:2^(2+i))).*Ref(x))))
end

b = FourierSine(2^(2+3),1+x^2)
save("fouriersine3.pdf",graylines(b,3,rasterize=2))

out = graylines(1+x^2,3,rasterize=2)
graylines!(fsum(b),3,rasterize=2)
save("fouriersine4.pdf",out)

function FourierSine(N,M,f)
    xy = TensorField(base(f))
    x,y = getindex.(xy,1),getindex.(xy,2)
    L1 = points(x).v[1][end]-points(x).v[1][1]
    L2 = points(x).v[2][end]-points(x).v[2][1]
    ω1 = Cartan.FourierSpace((π/L1)*(1:N),points(x).v[1])
    ω2 = Cartan.FourierSpace((π/L2)*(1:N),points(x).v[1])
    TensorField(ProductSpace(ω1,ω2),
        [(4/(L1*L2))*integrate(f*sin(n*x)*sin(m*y)) for n ∈ ω1, m ∈ ω2])
end

xy = TensorField(ProductSpace(0:0.01:2,0:0.01:3))
x,y = getindex.(xy,1),getindex.(xy,2)
f = exp(-(1-x)^2*(1-y)^2)*sin(π*x)*sin(2π*y)+exp(-(1.5-x)^2*(2.3-y)^2)*sin(2π*x)*sin(π*y)
b = FourierSine(10,10,f)
S = sum([fiber(b)[n,m]*(sin((n*π/2)*x)*sin((m*π/3)*y)) for n ∈ 1:10, m ∈ 1:10])

save("fouriersine5.pdf",contour(b,colormap=:grays,linewidth=3))
save("fouriersine6.pdf",contour(S,colormap=:grays),rasterize=2,linewidth=3)
save("fouriersine7.pdf",contour(norm(f-S),colormap=:grays),rasterize=2,linewidth=3)

t = TensorField(-10:1:10)
tt = TensorField(-10.5:0.01:10.5)
out = graylines(sinc(tt),3,rasterize=2)
scatter!(hat(t),color=:black)
save("spike1.pdf",out)

out = graylines(hat(t),3,rasterize=2)
scatter!(hat(t),color=:black)
save("spike2.pdf",out)

save("freq1.pdf",graylines(fftspace(0:0.01:10),3,rasterize=2,color=:gray))
save("freq2.pdf",graylines(rfftspace(0:0.01:10),3,rasterize=2,color=:gray))

t = TensorField(0:0.01:5);
saw(t,T) = ((t+T/2)%T)-T/2
save("fft1.pdf",graylines(saw(t,1/2),3,rasterize=2))
save("fft2.pdf",graylines(abs(fft(saw(t,1/2))),3,rasterize=2))

freq = TensorField(fftspace(0:0.01:10))
save("fft3.pdf",graylines(hat(freq,30pi/5)+hat(freq,10pi),3,rasterize=2))
save("fft4.pdf",graylines(ifft(hat(freq,3pi/5)+hat(freq,pi)),3,rasterize=2))

t = TensorField(0:0.01:5);
save("rfft1.pdf",graylines(saw(t,1/2),3,rasterize=2))
save("rfft2.pdf",graylines(abs(rfft(saw(t,1/2))),3,rasterize=2))
save("rfft3.pdf",graylines(angle(rfft(saw(t,1/2))),3,rasterize=2))
save("rfft4.pdf",graylines(rfft(saw(t,1/2)),3,rasterize=2))

freq = TensorField(rfftspace(0:0.01:10))
save("rfft5.pdf",graylines(hat(freq,30pi/5)+hat(freq,10pi),3,rasterize=2))
save("rfft6.pdf",graylines(irfft(hat(freq,3pi/5)+hat(freq,pi)),3,rasterize=2))

save("rfft7.pdf",graylines(exp(-0.003(freq-100)^2),3,rasterize=2))
save("rfft8.pdf",graylines(irfft(exp(-3(freq-10)^2)),3,rasterize=2))

t = TensorField(0:0.01:5);
save("r2r1.pdf",graylines(saw(t,1/2),3,rasterize=2))
save("r2r2.pdf",graylines(FFTW.r2r(saw(t,1/2),FFTW.RODFT10)/2length(t),3,rasterize=2))

t = TensorField(0:0.01:5);
save("dct1.pdf",graylines(saw(t,1/2),3,rasterize=2))
save("dct2.pdf",graylines(dct(saw(t,1/2)),3,rasterize=2))

t = TensorField(0:0.01:5);
save("fftshift1.pdf",graylines(abs(fft(saw(t,1/2))),3,rasterize=2))
save("fftshift2.pdf",graylines(fftshift(abs(fft(saw(t,1/2)))),3,rasterize=2))

t = TensorField(-1:0.001:1)
out = graylines(heaviside(t),3,rasterize=2)
save("heaviside1.pdf",out)
save("heavisidefilter.pdf",graylines(irfft(rfft(saw(t+1,1))*heaviside(-rfftspace(t),-100)),3,rasterize=2))

XY = TensorField(ProductSpace{2}(-1:0.1:1,-1:0.1:1))
h3 = hat(XY,0,0)+hat(XY,0.5,-0.5)+hat(XY,-0.5,0.5)
save("spike3.pdf",mesh(graph(h3),-h3,colormap=:grays))
h2 = heaviside(-XY,0,0)
save("heaviside2.pdf",mesh(graph(h2),-h2,colormap=:grays))

t = TensorField(0:0.01:5)
sp = hat(t,1)+0.5hat(t,2)
save("convolve1.pdf",graylines(sp,3,rasterize=2))
save("convolve2.pdf",graylines(heaviside(t)-heaviside(t,2),3,rasterize=2))
save("convolve3.pdf",graylines(convolve(heaviside(t)-heaviside(t,2),sp),3,rasterize=2))
save("convolve4.pdf",graylines(cos(t),3,rasterize=2))
save("convolve5.pdf",graylines(convolve(cos(t),sp),3,rasterize=2))

t = TensorField(range(-pi,3pi,2*24+1)) # coarse grid
tt = TensorField(range(-pi,3pi,1000)) # fine grid
mysinc(t,h) = iszero(t[1]) ? one(t[1]) : sin(pi*t[1]/h)/(2pi/h*tan(t[1]/2))

out = graylines(mysinc.(tt,2pi/24),3,rasterize=2)
scatter!(mysinc.(t,2pi/24),color=:black)
save("mysinc1.pdf",out)

out = graylines(mysinc.(t,2pi/24),3,rasterize=2)
scatter!(mysinc.(t,2pi/24),color=:black)
save("mysinc2.pdf",out)

mysincprime(N,h) = vcat(0,0.5*(-1).^(1:N-1).*cot.((1:N-1)*h/2))
spectral(N,h) = Toeplitz(mysincprime(N,h),-mysincprime(N,h))
pts(N=24,h=2pi/N) = TensorField(h*(1:N))
save("mysinc3.pdf",graylines(spectral(100,2pi/100)*exp(sin(pts(100))),3,rasterize=2))
save("mysinc4.pdf",graylines(spectral(100,2pi/100)*exp(sin(pts(100)))-gradient(exp(sin(pts(100)))),3,rasterize=2))

pts(N=24,h=2pi/N) = TensorField((h*(1:N)).-pi)
save("gradfft1.pdf",graylines(gradient_impulse(pts(24)),3,rasterize=2))
fun(N) = exp(sin(pts(N)))
save("gradfft2.pdf",graylines(gradient_rfft(fun(24)),3,rasterize=2))
error(N,grad) = grad(fun(N[1])) - cos(pts(N[1]))*fun(N[1])
out = graylines(error(24,gradient),3,rasterize=2)
graylines!(error(24,gradient_rfft),3,rasterize=2)
save("gradfft3.pdf",out)

maxerror(N,grad) = fiber(maximum(abs(error(N[1],grad))))
loglog(x) = TensorField(log10.(points(x)),log10.(fiber(x)))
out = graylines(loglog(maxerror.(TensorField((2).^(2:12)),gradient)),3,rasterize=2)
graylines!(loglog(maxerror.(TensorField((2).^(0:12)),gradient_rfft)),3,rasterize=2)
save("gradfft4.pdf",out)

save("integralfft1.pdf",graylines(integral_impulse(pts(24)),3,rasterize=2))
save("integralfft2.pdf",graylines(integral_rfft(fun(24)),3,rasterize=2))

# Chebyshev

cheb = Chebyshev(-2:0.2:5)
out = scatter(TensorField(cheb,sin.(angle(cheb))),color=:gray)
scatter!(TensorField(cheb,0),color=:black)
save("chebyshev1.pdf",out)

x = TensorField(Chebyshev(20))
y = exp(x)*sin(5x)
save("chebyshev2.pdf",graylines(gradient_chebyshev(y),3,rasterize=2))
save("chebyshev3.pdf",graylines(gradient_chebyshev(y) - gradient_chebfft(y),3,rasterize=2))
save("chebyshev4.pdf",graylines(gradient_chebfft(y),3,rasterize=2))

x = TensorField(range(-pi,pi,25))
xc = TensorField(Chebyshev(range(-pi,pi,25)))
xl = TensorField(GaussLegendre(range(-pi,pi,25)))
out = scatter(exp(-xc^2),color=:gray)
graylines!(integral_clenshawcurtis(exp(-xc^2)),3,rasterize=2)
save("chebyshev5.pdf",out)

# Phasor

basis"2"
t = TorusParameter(7)
hp = Phasor(sqrt(3),v12)
hexphase = hp(t,pi/6)
hexvec = vectorize(complexify(hexphase))
out = lines(hexvec,color=:black) # lines(complexify(hp(t,pi/6)))
graylines!(sqrt(3)*unitcircle(60),3,rasterize=2)
save("hexphase1.pdf",out)

k = TensorField(ProductSpace{Submanifold(2)}(-5:5,-5:5))
ts = Endomorphism{Manifold(k)}([3 3/2; 0 sqrt(3)*3/2])
save("hexphase2.pdf",lines(vec(Ref(hexvec).+fiber(ts > k)),color=:black))

pt,pe = initmesh("circleg","hmax"=>0.5)
save("simplexhat1.pdf",wireframe(graphbundle(hat(pt,0,0)),color=:black))

pt,pe = initmesh("circleg","hmax"=>0.3)
save("simplexhat2.pdf",wireframe(graphbundle(hat(pt,0.3,-0.3)),color=:black))

