# https://github.com/chakravala
# math: Duren

# Chapter 1

Naturals
sum(Naturals)
PositiveRationals

rn = FunctionVector(^,100)
rn(0.1)
rn(0.9)
rn(1.1)

isconverging(rn(0.1))
isconverging(rn(1.1))

fun(x,k) = x^inv(k)
seq = FunctionVector(fun,100)
seq(0.1)
seq(0.2)
isconverging(seq(0.1))

# Chapter 6

# Interpolation

chebyshevpoly(x,n) = chebyshevpoly(x,Val(n))
chebyshevpoly(x,::Val{0}) = 1
chebyshevpoly(x,::Val{1}) = x
chebyshevpoly(x,::Val{N}) where N = 2x*chebyshevpoly(x,Val(N-1))-chebyshevpoly(x,Val(N-2))

chebpoly(x,n) = cos(n*acos(x))




