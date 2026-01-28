# Makie.jl

f = Figure(size = (800, 800))
Axis(f[1, 1])

xs = LinRange(0, 2pi, 20)
ys = LinRange(0, 3pi, 20)
us = [sin(x) * cos(y) for x in xs, y in ys]
vs = [-cos(x) * sin(y) for x in xs, y in ys]
xy = TensorField(OpenParameter(xs,ys),Chain.(us,vs))
strength = vec(fiber(norm(xy)))

arrows2d!(xy, lengthscale = 0.2, color = strength, colormap=:grays, rasterize=3)

save("arrows2d.pdf",f)

ps = OpenParameter(-5:2:5,-5:2:5,-5:2:5)
ns = map(p -> 0.1 * Chain(p[2], p[3], p[1]), ps)
save("arrows3d1.pdf",arrows(
    TensorField(ps, ns),
    shaftcolor = :gray, tipcolor = :black,
    align = :center, axis=(type=Axis3,), rasterize=3
   ))

lengths = vec(norm.(ns))
save("arrows3d2.pdf",arrows(
    TensorField(ps, ns), color = lengths, lengthscale = 1.5,
    align = :center, axis=(type=Axis3,), colormap=:grays, rasterize=3
))

f = Figure(); Axis(f[1, 1])
xs = LinRange(0, 10, 100)
ys = LinRange(0, 15, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]
xyz = TensorField(OpenParameter(xs,ys),zs)
contour!(xyz,colormap=:grays,rasterize=3,linewidth=3)
contour!(xyz,levels=-1:0.3:1,colormap=:grays,rasterize=3,linewidth=3)
f
save("contour1.pdf",f)

himmelblau(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
x = y = range(-6, 6; length=100)
z = himmelblau.(x, y')

levels = 10.0.^range(0.3, 3.5; length=10)
colorscale = ReversibleScale(x -> x^(1 / 10), x -> x^10)
xyz = TensorField(OpenParameter(x,y),z)
f, ax, ct = contour(xyz; labels=true, levels, colormap=:grays, colorscale,rasterize=3,linewidth=3)
save("contour2.pdf",f)

x = -10:10
y = -10:10
# The curvilinear grid:
xs = [x + 0.01y^3 for x in x, y in y]
ys = [y + 10cos(x/40) for x in x, y in y]

# Now, for simplicity, we calculate the `zs` values to be
# the radius from the center of the grid (0, 10).
zs = sqrt.(xs .^ 2 .+ (ys .- 10) .^ 2)
# We can use Makie's tick finders to get some nice looking contour levels:
levels = Makie.get_tickvalues(Makie.LinearTicks(7), extrema(zs)...)
xy = TensorField(ProductSpace(x,y),Chain.(xs,ys))
xyz = TensorField(ProductSpace(x,y),zs)
# and now, we plot!
fig, ax, srf = mesh(xy,xyz; shading = NoShading, colormap = :grays, rasterize=3,
    axis = (; type = Axis, aspect = DataAspect()))
ctr = contour!(xy,xyz; levels = levels, color=:black, rasterize=3,
    labels = true, labelfont = :bold, labelsize = 12,linewidth=3)
save("contour3.pdf",fig)

f = Figure()
ax1 = Axis(f[1, 1])
ctrf1 = contourf!(TensorField(ProductSpace(x,y)), xyz; levels = levels,colormap=:grays,rasterize=3)
save("contourf1.pdf",f)

f = Figure()
ax2 = Axis(f[1, 1])
ctrf2 = contourf!(xy,xyz; levels = levels,colormap=:grays,rasterize=3)
save("contourf2.pdf",f)

r = range(-pi, pi, length = 21)
data2d = [cos(x) + cos(y) for x in r, y in r]
data3d = [cos(x) + cos(y) + cos(z) for x in r, y in r, z in r]

f = Figure(size = (700, 400))
a1 = Axis3(f[1, 1], title = "3D contour()")
rrr = OpenParameter(r,r,r)
contour!(TensorField(rrr, data3d),colormap=:grays,rasterize=3)
save("contour4.pdf",f)

f = Figure(size = (700, 400))
a2 = Axis3(f[1, 1], title = "contour3d()")
rr = OpenParameter(r,r)
contour3d!(TensorField(rr,data2d), linewidth = 3, levels = 10,colormap=:grays)
save("contour5.pdf",f)

f = Figure(size = (700, 300))
a1 = Axis3(f[1, 1])
contour!(TensorField(rrr,data3d), isorange = 0.04,colormap=:grays)
save("contour6.pdf",f)

f = Figure(size = (700, 300))
# small alpha can be used to see into the contour plot
a2 = Axis3(f[1, 1])
contour!(TensorField(rrr,data3d), data3d, alpha = 0.05,colormap=:grays)
save("contour7.pdf",f)

f = Figure()
Axis3(f[1, 1], aspect=(0.5,0.5,1), perspectiveness=0.75)

xs = ys = LinRange(-0.5, 0.5, 100)
zs = [sqrt(x^2+y^2) for x in xs, y in ys]

xyz = TensorField(OpenParameter(xs,ys),zs)

contour3d!(-xyz, linewidth=2, color=:black,rasterize=2)
contour3d!(+xyz, linewidth=2, color=:gray,rasterize=2)

save("contour3d1.pdf",f)

f = Figure()
Axis3(f[1, 1], aspect=(0.5,0.5,1), perspectiveness=0.75)
contour3d!(-xyz, levels=-(.025:0.05:.475), linewidth=2, color=:black)
contour3d!(+xyz, levels=  .025:0.05:.475,  linewidth=2, color=:gray)

save("contour3d2.pdf",f)

f = Figure()
ax = Axis(f[1, 1])

centers_x = [1, 2, 4, 7, 11]
centers_y = [6, 7, 9, 12, 16]
data = reshape(1:25, 5, 5)
xy = ProductSpace(centers_x,centers_y)

heatmap!(TensorField(xy,data),colormap=:grays)
scatter!(TensorField(xy,collect(xy)), color=:white, strokecolor=:black, strokewidth=1)
save("heatmap0.pdf",f)

xs = range(0, 2pi, length=100)
ys = range(0, 2pi, length=100)
zs = [sin(x*y) for x in xs, y in ys]
xyz = TensorField(OpenParameter(xs,ys),zs)

fig, ax, hm = heatmap(xyz,colormap=:grays)
Colorbar(fig[:, end+1], hm)
save("heatmap1.pdf",fig)

x = 10.0.^(1:0.1:4)
y = 1.0:0.1:5.0
z = broadcast((x, y) -> x - 10, x, y')
xyz = TensorField(OpenParameter(x,collect(y)),z)

scale = ReversibleScale(x -> asinh(x/2) / log(10), x -> 2sinh(log(10) * x))
fig, ax, hm = heatmap(xyz; colorscale = scale, axis = (; xscale = scale),colormap=:grays)
Colorbar(fig[1, 2], hm)
save("heatmap2.pdf",fig)

f = Figure()
Axis(f[1, 1])

xs = TensorField(1:0.2:10)
ys = sin(xs)

linesegments!(ys,rasterize=2,colormap=:grays)
linesegments!(ys - 1, linewidth = 5,rasterize=2,colormap=:grays)
linesegments!(ys - 2, linewidth = 5, color = LinRange(1, 5, length(xs)),rasterize=2,colormap=:grays)
save("linesegments1.pdf",f)

xs = range(0, 10, length = 30)
ys = 0.5 .* sin.(xs)
xy = TensorField(xs,ys)
save("scatter1.pdf",scatter(xy,color=:black))

xs = range(0, 10, length = 30)
ys = 0.5 .* sin.(xs)
pts = TensorField(xs,Chain.(xs, ys))
save("scatter2.pdf",scatter(pts, color = 1:30, markersize = range(5, 30, length = 30),
                            colormap = :grays))

struct FitzhughNagumo{T}
    e::T
    s::T
    y::T
    b::T
end

P = FitzhughNagumo(0.1, 0.0, 1.5, 0.8)

fun(x) = fun(x, P)
fun(x, P::FitzhughNagumo) = Chain(
    (x[1]-x[2]-x[1]^3+P.s)/P.e,
    P.y*x[1]-x[2] + P.b
)

xy = OpenParameter(-1.5:0.1:1.5,-1.5:0.1:1.5)

save("streamplot1.pdf",streamplot(fun.(xy), colormap = :grays,rasterize=3))

xs = LinRange(0, 10, 100)
ys = LinRange(0, 15, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]
xyz = TensorField(OpenParameter(xs,ys),zs)

save("surface0.pdf",surface(xyz, axis=(type=Axis3,), colormap=:grays,rasterize=3))

rs = 1:10
thetas = 0:10:360

xs = rs .* cosd.(thetas')
ys = rs .* sind.(thetas')
zs = sin.(rs) .* cosd.(thetas')
xyz = TensorField(ProductSpace(rs,thetas),Chain.(xs,ys,zs))

save("mesh1.pdf",mesh(xyz,TensorField(xyz,zs),colormap=:grays,rasterize=3))

xy = TensorField(ProductSpace(rs,thetas),Chain.(xs,ys))

save("mesh2.pdf",mesh(xy,TensorField(xy,zs),colormap=:grays,rasterize=3))

r = LinRange(-1, 1, 100)
rrr = OpenParameter(r,r,r)
cube = Real(abs2(rrr))
save("volume1.pdf",contour(abs2(rrr),alpha=0.5))

colormap = [:gray33, :transparent, :transparent,
           :gray66, :transparent, :black]

# Same as volume example
r = LinRange(-1, 1, 100)
rrr = OpenParameter(r,r,r)
cube = Real(abs2(rrr))
cube_with_holes = cube*(cube .> 1.4)

# To match the volume example with isovalue=1.7 and
# isorange=0.05 we map all values outside the
# range (1.65..1.75) to invisible air blocks with is_air
save("voxels1.pdf",voxels(cube_with_holes,
                          is_air = x -> !(1.65 <= x <= 1.75), colormap=:grays))

chunk = TensorField(OpenParameter(8,8,8),
    reshape(collect(1:512), 8, 8, 8))

f, a, p = voxels(chunk,
    colorrange = (65, 448), colorscale = log10,
    lowclip = :gray, highclip = :gray,
    colormap = :grays)

x, y = -8:0.5:8, -8:0.5:8
z = [sinc(sqrt(X^2 + Y^2) / pi) for X in x, Y in y]
xyz = TensorField(ProductSpace(x,y),z)
save("wireframe1.pdf",wireframe(graph(xyz), axis=(type=Axis3,), color=:black,rasterize=3))

