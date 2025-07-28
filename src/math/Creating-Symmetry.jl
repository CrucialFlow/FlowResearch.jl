# github.com/chakravala
# Math: Creating Symmetry, Farris

# 1 Going in Circles

t = TensorField(0:0.01:2pi)

# circle
y(t) = Chain(cos(t),sin(t))
exp(t*im)


# center (h,k) and radius R
Chain(h,k) + R*y(t)
h+k*im + R*exp(t*im)

# rolling quarter
y(t) = Chain(2sin(t)-sin(2t),2cos(t)-cos(2t))
roll = 2exp((-t+π/2)*im)+exp((-2t-π/2)*im)
arclength(y.(t))
arclength(roll)

# epicycloid
Chain((a+1)*sin(t) - sin((a+1)*t),(a+1)*cos(t)-cos((a+1)*t))
(a+1)*exp(-t+π/2)+exp((-(a+1)*t-π/2)*im)

# mystery curve

μ(t) = Chain(cos(t)+cos(6t)/2+sin(14t)/3,sin(t)+sin(6t)/2+cos(14t)/3)
exp(t*im)+exp(t*6im)/2 + exp(t*(-14im))*(im/3)

# ancient mathematics

y(t) = Chain(1-t[1]^2,2t[1])/(1+t[1]^2)

x[1]*exp(t*im)+x[2]*exp(t*6im) + exp(t*(-14im))*(im/3)

function fourier(a::Vector,t=TensorField(0:0.01:2pi))
    out = exp(t*im)
end

