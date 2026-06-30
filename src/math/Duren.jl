# https://github.com/chakravala
# math: Duren

chebyshevpoly(x,n) = chebyshevpoly(x,Val(n))
chebyshevpoly(x,::Val{0}) = 1
chebyshevpoly(x,::Val{1}) = x
chebyshevpoly(x,::Val{N}) where N = 2x*chebyshevpoly(x,Val(N-1))-chebyshevpoly(x,Val(N-2))

chebpoly(x,n) = cos(n*acos(x))




