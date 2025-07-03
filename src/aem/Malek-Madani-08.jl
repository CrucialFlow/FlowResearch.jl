# github.com/chakravala
# AEM: Malek-Madani
# 08 Introduction to Vector Calculus

# 8.13 Project B

f(x) = 1-x[1]^2-x[2]^2

ff = TensorField(OpenParameter(-1:0.1:1,-1:0.1:1),f)

surface(ff)
streamplot(graph(ff),-gradient(ff))

# 8.14 Project C

t = TensorField(0:0.01:2pi);
helix(x) = Chain(cos(x[1]),sin(x[1]),x[1])/sqrt(2)
lines(helix.(t))
scaledarrows!(helix.(t),frame(helix.(t)),gridsize=11)

curve(x) = Chain(sin(2x[1])*cos(3x[1]),sin(2x[1])*sin(3x[1]),cos(x[1]))
lines(curve.(t))
scaledarrows!(curve.(t),frame(curve.(t)),gridsize=11)

# 8.15 Project D

vfa(x) = Chain(x[2],-x[1])
vfb(x) = vfa(x)/norm(fiber(x))
vfc(x) = vfa(x)/abs2(fiber(x))

ini = Chain.(-0.4:0.1:-0.1,0)
lines(vfa,ini)
lines(vfb,ini)
lines(vfc,ini)

ini = TensorField(0:5,Chain.(-2:0.1:-1.5,0))
lines(Flow(vfa,0.2),ini,7)
lines(Flow(vfb,0.2),ini,7)
lines(Flow(vfc,0.2),ini,7)

t = TensorField(LinRange(0,2pi,20))
t = TensorField(0:0.01:2pi)
circ = Chain.(0.2cos(t)+1,0.2sin(t))
lines(Flow(vfa,0.2),circ,7)
lines(Flow(vfb,0.2),circ,7)
lines(Flow(vfc,0.2),circ,7)

d4a(x) = Chain(x[2],-sin(x[1]))
d4b(x) = Chain(x[2],-0.1x[2]-sin(x[1]))
d4c(x) = Chain(x[2],-(x[1]^2-1)*x[2]-x[1])

circ = Chain.(0.1cos(t)+π/4,0.1sin(t))
lines(Flow(d4a,0.2,circ,7))

circ = Chain.(0.1cos(t)+π/4,0.1sin(t))
lines(Flow(d4b,0.2,circ,7))

circ = Chain.(0.1cos(t)+0.5,0.1sin(t))
lines(Flow(d4c,0.2,circ,7))

van = d4c.(OpenParameter(-3:0.1:3,-3:0.1:3))
streamplot!(van)

# 8.16 Project E

sys1(λ) = x->Chain(λ+x[1]^2,-x[2])
sys2(λ) = x->Chain(-λ*x[1]-x[1]^3,-x[2])
sys3(λ) = x->Chain(x[2]+x[1]*(λ-x[1]^2-x[2]^2),-x[1]+x[2]*(λ-x[1]^2-x[2]^2))

dom = OpenParameter(-3:0.1:3,-2:0.1:2)
streamplot(sys1(-1).(dom))
streamplot(sys1(0).(dom))
streamplot(sys1(1).(dom))

streamplot(sys2(-1).(dom))
streamplot(sys2(0).(dom))
streamplot(sys2(1).(dom))

streamplot(sys3(-0.1).(dom))
streamplot(sys3(0).(dom))
streamplot(sys3(0.1).(dom))

ini = [Chain(-3.,-1.),Chain(-3.,-0.5),Chain(-3.,0.5),Chain(-3.,1.)]
lines(sys1(-1),ini)
lines!(sys1(0),ini)
lines!(sys1(1),ini)

lines(sys2(-1),ini)
lines!(sys2(0),ini)
lines!(sys2(1),ini)

ini = [Chain(-0.5,-1),Chain(-0.5,-0.5)]
lines(sys3(-0.1),ini)
lines!(sys3(0),ini)
lines!(sys3(0.1),ini)

# 8.17 Project F

ψ(x) = sin(π*x[1])*sin(π*x[2])
dom = OpenParameter(0:0.1:3,0:0.1:1)

contour(ψ.(dom))
streamplot(!gradient(ψ.(dom)))

ini = [Chain(0.5,0.1),Chain(0.5,0.2),Chain(0.5,0.3),
       Chain(1.5,0.1),Chain(1.5,0.2),Chain(1.5,0.3),
       Chain(2.5,0.1),Chain(2.5,0.2),Chain(2.5,0.3)]
lines(!gradient(ψ.(dom)),ini)

myfun!(x;args...) = streamplot!(!gradient(x);args...)
stream(ϵ=0.1,ω=0.1) = x->ψ(x)+ϵ*cos(ω*x[3])*cos(x[1])*sin(π*x[2])
dom3 = OpenParameter(0:0.1:3,0:0.02:1,0:0.1:5)

sd = stream(0.1,0.1).(dom3)
variation(sd,0.1,contourf,contourf!)
variation(sd,0.1,myfun,myfun!)

sd = stream(0.1,1.0).(dom3)
variation(sd,0.1,contourf,contourf!)
variation(sd,0.1,myfun,myfun!)

sd = stream(0.1,10.0).(dom3)
variation(sd,0.1,contourf,contourf!)
variation(sd,0.1,myfun,myfun!)



