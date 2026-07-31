# https://github.com/chakravala
# math: Judson

Z(n,f=(+),g=groupinverse(f)) = Semimagma(Mod{n}.(0:n-1),f,g)
U(n,Zn=Z(n,*)) = Semimagma(Zn[findall(is_invertible,Zn)])

# Chapter 1

Naturals
Integers
Rationals

# Example 3.2

Mods.value.(cayley(Z(8,*)))

# Table 3.7

cayley(SymmetricGroup(3))

# Table 3.10

Mods.value.(cayley(Z(5)))

# Example 3.11

Mods.value.(cayeley(U(8)))

# Example 3.15

Q4 = Semimagma(TensorOperator.([[1 0; 0 1],[0 1; -1 0],[0 im; im 0],[im 0; 0 -im]]))
Q8 = Semimagma([Q4.v...,-Q4.v...])

# Section 3.3

Semimagma(2Integers(10))

# Example 3.24

Semimagma(NonzeroRationals(10))

# Example 3.25

H = Semimagma([1,-1,im,-im])

# Exercise 3.4.6

Mods.value.(vayley(U(12)))

# Exercise 3.4.37

2^Naturals(20)

# Chapter 4

# Example 4.1

Semimagma(3Integers(20))

# Example 4.5

Zp(6)
group(Mod{6}(1),+)
group(Mod{6}(5),+)
group(Mod{6}(2),+)

# Example 4.6

iscyclic(U(9))

# Example 4.15

Zp(16)
group(Mod{16}(1),+)
group(Mod{16}(3),+)
group(Mod{16}(5),+)
group(Mod{16}(7),+)
group(Mod{16}(9),+)
group(Mod{16}(11),+)
group(Mod{16}(13),+)
group(Mod{16}(15),+)

# Theorem 4.25

unityroots(3)
unityroots(4)

# Exercise 4.4.1

iscyclic(U(8))

# Exercise 4.4.2

group(Mod{12}(5),+)
group(Mod{240}(72),+)
group(Mod{471}(312),+)

# Exercise 4.4.3

7Integers
group(Mod{24}(15),+)
U(20)
group(Mod{20}(3))
U(18)
group(Mod{18}(5),+)

group(1im)
group(2im)
group((1+im)/sqrt(2))
group((1+sqrt(3)*im)/2)

# Chapter 6

# Example 6.1

Z(6)
H = Semimagma(Mod{6}.([0,3]),+)
Z(6)/H

# Example 6.2

S3 = SymmetricGroup(3)

H = Semimagma(Permutation.(Cycle{3}.([Values(1),Values(1,2,3),Values(1,3,2)])))
leftcosets(H,S3)
rightcosets(H,S3)

K = Semimagma(Permutation.(Cycle{3}.([Values(1),Values(1,2)])))
leftcosets(K,S3)
rightcosets(K,S3)

# Example 6.6

H = Semimagma(Mod{6}.([0,3]),+,-)
Z(6)/H

# Example 10.2

S3 = SymmetricGroup(3)
H = Semimagma(Permutation.(Cycle{3}.([Values(1),Values(1,2)])))
N = Semimagma(Permutation.(Cycle{3}.([Values(1),Values(1,2,3),Values(1,3,2)])))

S3/N

