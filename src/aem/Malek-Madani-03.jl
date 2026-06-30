# https://github.com/chakravala
# AEM: Malek-Madani

# Project A

dt,dx = 0.3,0.2
t = TensorField(0:dt:3)
x = TensorField(-2:dt:2)

T = ones(length(x))'*points(t)
X = points(x)'*ones(length(t))
M = xprime(T,X)
Scale = 0.1./(2sqrt(1+M.*M))
Tleft,Tright = T-Scale,T+Scale
newScale=M.*Scale
Xleft=X-newScale
Xright=X+newScale
newTleft = Tfleft


# Project C

picardintegral(f) = x -> fiber(x)[1] + integral(f(x))
picardintegral(f,t) = x -> fiber(x)[1] + integral(f(x,t))

t = TensorField(0:0.01:1)
fprime(x,t) = sin(t)-2x
pic = picardintegral(fprime,t)
lines(fixedpoint(pic,1+8t,3))
lines(fixedpoint(pic,1+8t,6))
lines(fixedpoint(pic,1+8t,9))
lines(fixedpoint(pic,1+8t,12))
lines(fixedpoint(pic,1+8t,16))
lines(fixedpoint(pic,1+8t))

fprime(x) = -x*x
pic = picardintegral(fprime)
lines(fixedpoint(pic,1+0t))


