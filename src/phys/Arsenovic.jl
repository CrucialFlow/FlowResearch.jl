# Arsenovic notebook
# https://github.com/arsenovic/notebook

V = Submanifold(S"-+++")
D = Λ(V)

K =  rand(Chain{V,3,Float64})
x = rand(Chain{V,1,Float64})
F0 = rand(Chain{V,2,Float64})

k = !K

# Plane Waves

bm(x::LocalTensor) = bm(fiber(x))
bm(x::Chain{V,G}) where {V,G} = Chain{Submanifold(D"1,1,1,0"),G}(value(x))
bm(x::Spinor) = Spinor{Submanifold(D"1,1,1,0")}(value(x))
bm(x::CoSpinor) = CoSpinor{Submanifold(D"1,1,1,0")}(value(x))

up(x) =  !(x+Λ(Manifold(x)).v1)

function make_x(v23_N, v23_bounds, v1_N, v1_bounds, v4_N, v4_bounds;V=S"-+++")
    D = Λ(V)
    x1_range = range(-v1_bounds, v1_bounds, v1_N)
    x2_range = range(-v23_bounds, v23_bounds, v23_N)
    x3_range = range(-v23_bounds, v23_bounds, v23_N)
    x4_range = range(-v4_bounds, v4_bounds, v4_N)
    TensorField(ProductSpace{V}(x1_range,x2_range,x3_range,x4_range))
end

function render_frames(Fx;vector_scale=1)
    D = Λ(Manifold(Fx))
    E = (Fx>D.v1)*D.v1
    H = !(((!D.v1)>Fx)*!D.v1)
    vmag = vector_scale/size(Fx)[2]
    render_frame(E,vmag),render_frame(H,vmag)
end

render_frame(EH,vmag) = TensorField(base(EH),render_frame.(points(EH),fiber(EH),vmag))
function render_frame(x,EH,vmag)
    D = Λ(Manifold(EH))
    x_ = up(bm((D.v234>x)*D.v234))
    EH_ = bm(EH)
    exp(vmag*EH_)>>>x_
end

function F(x,V=S"-+++")
    D = Λ(V)
    F0 = D.v12 + D.v24
    K = (D.v234 + 0.01D.v123)#*exp(-0.1D.v1234)
    xK = x*K
    return exp(xK) >>> F0
end

function A(x,V=S"-+++")
    D = Λ(V)
    A0 = D.v4
    # K = x -> exp((0.1D.v1|x)*D.v23)>>(D.v234 + 0.5D.v124)#*exp(-1D.v1234)
    K = x -> (D.v23/D.v4 + 0.5D.v24/D.v1)#*exp(-1D.v1234)
    xK = x*K(x)
    return exp(xK)>>>A0
end

x = make_x(21,6,100,4,1,0)

#Fx = F.(x) # make F-field directly
Fx = d(A.(x))
E,H = render_frames(Fx,vector_scale=2)

