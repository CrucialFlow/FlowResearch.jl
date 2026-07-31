# https://github.com/chakravala
# math: Wade

# Chapter 1

# Example 1.23

2.0^(1-Naturals)
(2Naturals-1)/(2Naturals)

# Chapter 2

2.0^(1-Naturals)
(-1)^Naturals
Naturals

# Example 2.2

1/Naturals

# Exercises

3+1/Naturals
2*(1-1/Naturals)
CountableVector(n->(5+n)/n^2)
float(π)-3/sqrt(Naturals)

CountableVector(n->n+n*(-1)^(3n))

# Example 2.10

CountableVector(n->cos(n^3-n^2+n-13)*2.0^-n)

# Example 2.13

CountableVector(n->(n^3+n^2-1)/(1-3n^3))

# Exercises

limit(CountableVector(n->sin((n^4+n+1)/(n^2+1))/n),1e-7)
limit(CountableVector(n->n/(n^2+1)),1e-7)
limit(CountableVector(n->(sqrt(2n)+1)/(n+1)),1e-7)
limit(CountableVector(n->n/2.0^n),1e-7)

limit(CountableVector(n->(1+n-3n^2)/(3-2n+n^2)),1e-7)
limit(CountableVector(n->(n^3+n-5)/(5n^3+n-1)),1e-7)
limit(CountableVector(n->sqrt(2n^2-1)/(n+1)),1e-7)
limit(CountableVector(n->sqrt(n+1)-sqrt(n)),1e-7)

orbit(x->sqrt(2+x),0.1,eps())
collect(ans)

# Exercises

lines(CountableVector(n->(n^2+20n+35)*sin(n^3)/(n^2+n+1)))

orbit(x->1-sqrt(1-x),0.1,1e-7)
residuals(ans,/)

orbit(x->sqrt(2x+3),0.1,1e-7)
collect(ans)
isincreasing(ans)

orbit(x->1+sqrt(x-1),3.1,1e-7)
collect(ans)
isdecreasing(ans)

orbit(x->1+sqrt(x-1),1.7,1e-7)
collect(ans)
isincreasing(ans)

xk(x,n) = x^(1/(2n-1))
xkv = FunctionVector(xk)
xkv(0.1)
xkv(0)
xkv(-0.1+0im)

orbit(x->(1+x)/2,0.1,1e-7)

orbit(z->((z[1]+z[2])/2,sqrt(z[1]*z[2])),(2.0,1.0),100)

orbit(z->(z[1]+2z[2],z[1]+z[2]),(1.0,0.0),100)
/(last(ans)...)
sqrt(2)-ans

tupnorm1(x,y) = supnorm(x[1],y[1])
tupnorm2(x,y) = supnorm(x[2],y[2])
orbit(z->(x = (2z[1]*z[2])/(z[1]+z[2]); (x,sqrt(x*z[2]))),(2sqrt(3),3.0),100,tupnorm1)

# Remark 2.31

log(Naturals)
residuals(log(Naturals))

# Exercises

limit(sum(CountableVector(k->(-1)^k/k)),1e-5)

surface(10log(TensorField(ProductSpace(1:100,1:100),residualproduct(cumsum(CountableVector(k->(-1)^k/k))))))

sq = supseq((-1)^Naturals)
iq = infseq((-1)^Naturals)

limsup((-1)^Naturals)
liminf((-1)^Naturals)

sq = supseq(1+1/Naturals)
iq = infseq(1+1/Naturals)

limsup(1+1/Naturals)
liminf(1+1/Naturals)

surface(10log(TensorField(ProductSpace(1:100,1:100),residualproduct(1+1/Naturals))))

# Exercies

limsup(3-(-1)^Naturals)
liminf(3-(-1)^Naturals)

limsup(cos(Naturals*(π/2)))
liminf(cos(Naturals*(π/2)))

limsup(CountableVector(n->(-1)^(n+1)+(-1)^n/n))
liminf(CountableVector(n->(-1)^(n+1)+(-1)^n/n))

limsup(CountableVector(n->sqrt(1+n^2)/(2n-5)))
liminf(CountableVector(n->sqrt(1+n^2)/(2n-5)))

limsup(CountableVector(n->n*(1+(-1)^n)+((-1)^n-1)/n))
liminf(CountableVector(n->n*(1+(-1)^n)+((-1)^n-1)/n))

limsup(CountableVector(n->(n^3+n^2-n+1)/(n^2+2n+5)))
liminf(CountableVector(n->(n^3+n^2-n+1)/(n^2+2n+5)))

# Chapter 3




