# https://github.com/chakravala
# AEM: Malek-Madani

# Example 3.1.3

solve313(t,c1=1,c2=1,k=1,m=1,g=1) = c1*exp(-k*t/m)-(m*g/k)*t+c2

function verify313(t,c1=1,c2=1,k=1,m=1,g=1)
    y = solve313(t,c1,c2,k,m,g)
    gradient(gradient(y)) + k*gradient(y)+m*g
end

# Problem 1

t = TensorField(0:0.01:2)
xt = TensorField(ProductSpace{2}(0:0.01:2,0:0.01:2))

# 1(a)
y = 3exp(-2t)
lines(y'+2y)

# 1(b)
y = t^5/60 + t^3/6 + t^2/2 + t + 1
lines(y''' - t^2-1)

# 1(c)
y = -7cos(2t)
lines(y'' + 4y)

# 1(d)
y = t*exp(-2t)
lines(y''+4y'+4y)

# 1(e)
y = exp(-2t)*sin(t)
lines(y''+4y'+5y)

# 1(f)
u(xt,a=1) = exp(-a^2*xt[2])*sin(a*xt[1])
surface(gradient(u.(xt),2) - gradient(gradient(u.(xt),1),1))

# 1(g)
u(xt,a=1) = exp(-9a^2*xt[2])*cos(a*xt[1])
surface(gradient(u.(xt),2) - 9gradient(gradient(u.(xt),1),1))

# Example 3.1.7

mixic(V0=40,ri=3,ci=10,ro=2.9) = IC(Flow(x->ri*ci-fiber(x)/(V0+(ri-ro)*point(x)),2pi),0.0)
odesolve(mixic(40,3,10,2.9),ExplicitIntegrator{4}(2e-4))

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

picardintegral(f) = (x -> fiber(x)[1] + integral(f(x)))∘localfiber
picardintegral(f,t) = (x -> fiber(x)[1] + integral(f(x,t)))∘localfiber

t = TensorField(0:0.01:1)
fprime(x,t) = sin(t)-2x
pic = picardintegral(fprime,t)
lines(orbit(pic,1+8t,3))
lines(orbit(pic,1+8t,6))
lines(orbit(pic,1+8t,9))
lines(orbit(pic,1+8t,12))
lines(orbit(pic,1+8t,16))
lines(orbit(pic,1+8t))

fprime(x) = -x*x
pic = picardintegral(fprime)
lines(orbit(pic,1+0t))


