# https://github.com/chakravala
# math: Rogawski-Adams

# Chapter 5

# Example 5.1

fib(f,n) = n<3 ? 1 : f[n-1] + f[n-2]
fs = SequenceVector([1,1],fib)
fs[10]
fs

# Example 5.2

rec(a,n) = (a[n-1]+2/a[n-1])/2
rs1 = SequenceVector([1//1],rec)
rs1[6]
rs1

rs2 = SequenceVector([1.0],rec)
rs2[7]
rs2
lines(rs2)

# Example 5.3

seq(n) = (n+4)/(n+1)
an = CountableVector(seq,15)
lines(an)

lines(CountableVector(cos,15))

# Example 5.4

seq(n) = (n^2-2)/n^2
an = CountableVector(seq)
lines(an)

# Example 5.5

seq(n) = (n+log(n))/n^2
an = CountableVector(seq)

# Example 5.6

seq(n) = 365.5n^2/(n^2-4)
an = CountableVector(seq)

# Example 5.7

funs(r,n) = r^n
fs = FunctionVector(funs,100)




