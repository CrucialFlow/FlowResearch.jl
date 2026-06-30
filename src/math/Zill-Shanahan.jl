# https://github.com/chakravala
# Math: A First Course in Complex Analysis with Applications, Zill-Shanahan

z = complexify(TensorField(ProductSpace(-π/2:0.1:3π/2,-π/2:0.1:π/2)))
exp(z)

riemannsurface(z) = Chain(z[1],z[2],z[1])

