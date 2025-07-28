# github.com/chakravala
# DGA: Differentiable Manifolds, Conlon

# 4.5.16

f(x) = (1-x[1]^2-x[2]^2)*exp(x[3])
set = fol.(OpenParameter(-1.5:0.1:1.5,-1.5:0.1:1.5,0:0.1:2))
contour(set)

