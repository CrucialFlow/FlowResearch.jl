# github.com/chakravala

# Chapter 5

# Hyperbolic explicit

c = 1

x = TensorField(range(0,2,41))
icfun(x) = x[1] ≤ 1/2 ? 1.0 : 1/2

ic = icfun.(x)

Δt(x,CFL=0.9,c=1) = CFL/c*step(points(x))

function xt(x,CFL=0.9,c=1)
    dt = Δt(x,CFL,c)
    TensorField(points(x)⊕range(0,9dt,10))
end


