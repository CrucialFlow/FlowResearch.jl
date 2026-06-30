# https://github.com/chakravala
# dsp: Chui

# Chapter 2

gaussian(t,α) = exp(t^2/(-4α))/(2sqrt(π*α))

t = TensorField(range(-π,π,1001))
lines(gaussian(t,1))
lines!(gaussian(t,1/4))
lines!(gaussian(t,1/16))

# Chapter 3

gabor(t,α,b,ω) = exp((im*ω)*t)*gaussian(t-b,α)
lines(-real(gabor(t,0.2925,0,2π)))
lines(imag(gabor(t,0.23,0,2π)))

x = TensorField(0:0.01:5)
lines(ChebyshevFirst(1+x^2))
b = FourierSine(1+x^2,15)
lines(b)
lines(FourierSine(b))
lines!(1+x^2)

XY = TensorField(ProductSpace{2}(0:0.01:2,0:0.01:3))
X,Y = split(XY)
g = exp(-(1-X)^2*(1-Y)^2)*sin(pi*X)*sin(2pi*Y) +
exp(-(1.5-X)^2*(2.3-Y)^2)*sin(2pi*X)*sin(pi*Y)
b = FourierSine(g,10,10)
contour(b)
contour(FourierSine(b))
contour(norm(g-FourierSine(b)))


