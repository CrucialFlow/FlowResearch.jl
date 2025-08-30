# https://github.com/chakravala

# hexphase

basis"2"

t = TorusParameter(7)
hp = Phasor(sqrt(3),v12)
hexphase = hp(t,π/6)
hexvec = vectorize(complexify(hexphase))
lines(hexvec)

lines(complexify(hp(TorusParameter(7),π/6)))

k = TensorField(ProductSpace{Submanifold(2)}(-5:5,-5:5))
ts = Endomorphism{Manifold(k)}([3 3/2; 0 sqrt(3)*3/2])
lines(vec(Ref(hexvec).+fiber(ts⋅k)))

vectorize(sqrt(3)*exp(t*v12+π/6*v12))
vectorize(complexify(Phasor(sqrt(3),v12)(t,π/6)))

Phasor(sqrt(3),v12)(t,π/6) # |> complexify |> vectorize
Phasor(sqrt(3),v12)(π/6,t) # |> complexify |> vectorize

# Riemann sphere

xy(x) = Chain(x[1],x[2])/(1-x[3])
uvw(x) = Chain(x[1],x[2],x[1]^2+x[2]^2)/(x[1]^2+x[2]^2+1)

riemannline(x,n=11) = resample(TensorField(0:2,[Chain(x[1],x[2],0.0),Chain(x[1],x[2],x[1]^2+x[2]^2)/(x[1]^2+x[2]^2+1),Chain(0,0,1.0)]),n)

rs(n=61,r=55) = sectorize(TensorField(exp.(LinRange(-exp(1),exp(1),r)),unitcircle(n)))

wireframe(uvw.(rs()))
wireframe!(rs())
lines(riemannline(Chain(0.75,0.75)),color=:black)

