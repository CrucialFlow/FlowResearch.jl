# https://github.com/chakravala

# Example 1.4

V = Submanifold(2)
D = TensorField(ProductSpace{V}(LinRange(-π/2,π/2,100),LinRange(-π/2,π/2,100)));
scherk(x) = (z = cos(x[2])/cos(x[1]); Chain(x[1],x[2],bound(z≤0 ? 0.0 : log(z),π/2)))
mat = TensorOperator{V,Submanifold(3)}([0 0; π 0; 0 π])
k1 = TensorField(ProductSpace{Submanifold(2)}(-3:2:3,-3:2:3))
k2 = TensorField(ProductSpace{Submanifold(2)}(-2:2:2,-2:2:2))
mesh(Ref(scherk.(D)).+vec(fiber(mat⋅k1)),norm)
mesh!(Ref(scherk.(D)).+vec(fiber(mat⋅k2)),norm)

# Example 1.6

wireframe(revolve(unitcircle()+Chain(sqrt(2),0)))

# Example 1.7

t = TensorField(-2:0.1:2)
wireframe(revolve(Chain.(cosh(t),t+0)))

# Sierpinski

function sierpinski(n)
    pts = LinRange(0,1,4)
    out = lines(TensorField(ProductSpace(pts,pts)),color=:black)
    display(out)
    for i ∈ 1:n
        sierpinski!(i,pts)
        display(out)
    end
end
function sierpinski!(i,pts=0:0.25:1)
    pts2 = (pts/3^i).+inv(3^i)
    for j ∈ 1:3^(i-1)
        x = ((j-1)/(3^(i-1)))
        for k ∈ 1:3^(i-1)
            y = ((k-1)/(3^(i-1)))
            mesh!(TensorField(ProductSpace(pts2.+x,pts2.+y)),color=:black)
        end
    end
end

sierpinski(3)

