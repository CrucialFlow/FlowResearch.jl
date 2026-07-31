# https://github.com/chakravala
# math: Gallian

# Chapter 2

Z(n,f=(+),g=groupinverse(f)) = Semimagma(Mod{n}.(0:n-1),f,g)
U(n,Zn=Z(n,*)) = Semimagma(Zn[findall(is_invertible,Zn)])

# Example 2.1

Integers
Rationals

# Example 2.3

Semimagma([1,-1,im,-im])

# Example 2.4

PositiveRationals

# Example 2.7

Z(10)

# Example 2.11

U(10)

value.(cayley(U(10)))

# Example 2.12

isgroup(Z(4,*))

# Exercise 2.4

isgroup(Semimagma(Mod{16}.([0,4,8,12])))
isgroup(Semimagma(Mod{15}.([0,4,8,12])))
isgroup(Semimagma(Mod{15}.([1,4,7,13])))
isgroup(Semimagma(Mod{9}.([1,4,5,7])))

# Exercise 2.9

isgroup(Semimagma(Mod{4}.(1:3)))
isgroup(Semimagma(Mod{5}.(1:4)))

# Exercise 2.16

isgroup(Semimagma(Mod{40}.([5,15,25,35])))
issemigroup(Semimagma(Mod{40}.([5,15,25,35])))
U(8)

# Exercise 2.21

group(Semimagma(Mod{91}.([1,9,16,22,53,74,79,81])))
ans[end]

# Exercise 2.29

group(Semimagma(Mod{56}.([5,15])))

# Exercise 2.32

value.(cayely(U(12)))

# Chapter 3

# Example 3.1

U(15) .=> orders(U(15))

# Example 2.2

Z(10) .=> orders(Z(10))

# Example 3.9

U(10) == group(Mod{10}(3))

# Example 3.10

group(Mod{10}(2),+)

# Example 3.13

group(Mod{20}.([8,14]),+)
group([1,im])

# Exercise 3.1

Z(12) .=> orders(Z(12))
U(10) .=> orders(U(10))
U(12) .=> orders(U(12))
U(20) .=> orders(U(20))

# Exercise 3.22

U(14) == group(Mod{14}(3)) == group(Mod{14}(5))
iscyclic(U(14))

# Exercise 3.23

iscyclic(U(20))

# Exercise 3.62

order.(U.([3,4,12]))
order.(U.([5,7,35]))
order.(U.([4,5,20]))
order.(U.([3,5,15]))

# Exercise 3.64

order.(U.([4,10,40]))





