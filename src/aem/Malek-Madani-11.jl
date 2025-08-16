# github.com/chakravala
# AEM: Malek-Madani
# 10 Introduction to COmplex Variables

## Curves in the Complex Domain

z = complexify(TensorField(ProductSpace(-3:0.03:3,-3:0.03:3)))

# 11.3.1

contour(z-1;levels[2])

# 11.3.2

contour((z-1)/(z+2);levels=[3])
contour((z^2-1)/(z+2);levels=[1])
contour((z^3-2z+4)/(z^4-4z+8);levels=[1/2])
contour((z^3-2z+4)/(z^4-4z+8);levels=[3/2])

# Complex-Valued Functions

# Example 11.4.5

f(z) = z^2-3z+4

