# https://github.com/chakravala
# Math: Signal Processing First,McClellan, Schafer, Yoder

# Figure 9-11

t = TensorField(0:0.01:5)
hh = heaviside(t)-heaviside(t,2)
sp = spike(t,1)+0.5spike(t,2)
lines(convolve(hh,sp))


