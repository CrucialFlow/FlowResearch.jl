
# Γ(z,x)

Z = complexify(TensorField(ProductSpace(-5:0.01:5,-5:0.01:5)))

save("gamma1.pdf",mesh(graph(bound(norm(gamma(-Z)),10.0)),angle(gamma(-Z)),colormap=:grays,rasterize=10))

t = TensorField(-5:0.01:4)
save("gamma2.pdf",graylines(bound(real(gamma(t+0im)),5.0),3,rasterize=2))
save("gamma3.pdf",graylines(bound(real(loggamma(t+0im)),5.0),3,rasterize=2))

save("digamma1.pdf",mesh(graph(bound(norm(digamma(-Z)),10.0)),angle(digamma(-Z)),colormap=:grays,rasterize=8))
save("trigamma1.pdf",mesh(graph(bound(norm(trigamma(-Z)),10.0)),angle(trigamma(-Z)),colormap=:grays,rasterize=6))

t = TensorField(-5:0.0001:4);
save("digamma2.pdf",graylines(bound(real(digamma(t+0im)),5.0),3,rasterize=2))
save("trigamma2.pdf",graylines(bound(real(trigamma(t+0im)),5.0),3,rasterize=2))

# B(x,y)

XY = TensorField(ProductSpace(-3:0.01:3,-3:0.01:3));
fun(x) = beta(x[1],x[2])
save("beta1.pdf",contourf(log10(norm(fun.(XY))),levels=50,colormap=:grays))

t = TensorField(0:0.01:1)
f(t,x,y) = t^(x-1)*(1-t)^(y-1)/beta(x,y)
out = graylines(f(t,0.5,0.5),3,rasterize=2)
graylines!(f(t,5,1),3,rasterize=2)
graylines!(f(t,1,3),3,rasterize=2)
graylines!(f(t,2,2),3,rasterize=2)
graylines!(f(t,2,5),3,rasterize=2)
save("beta2.pdf",out)

t = TensorField(0:0.01:1)
out = graylines(ncbeta(0.5,0.5,0,t),3,rasterize=2)
graylines!(ncbeta(5,1,0,t),3,rasterize=2)
graylines!(ncbeta(1,3,0,t),3,rasterize=2)
graylines!(ncbeta(2,2,0,t),3,rasterize=2)
graylines!(ncbeta(2,5,0,t),3,rasterize=2)
save("ncbeta.pdf",out)

# expint

z = complexify(TensorField(ProductSpace(-2:0.01:2,-2:0.01:2)));
save("expint.pdf",mesh(graph(bound(norm(expint(z)),7)),angle(expint(z)),colormap=:grays,rasterize=5))

t = TensorField(-6pi:0.01:6pi)
save("sinint.pdf",graylines(sinint(t),3,rasterize=2))

t = TensorField(0:0.01:6pi)
save("cosint.pdf",graylines(bound(cosint(t),1),3,rasterize=2))

# erf

z = complexify(TensorField(ProductSpace(-3:0.01:3,-3:0.01:3)));
save("erf.pdf",mesh(graph(bound(norm(erf(z)),3)),angle(erf(z)),colormap=:grays,rasterize=5))
save("erfc.pdf",mesh(graph(bound(norm(erfc(z)),3)),angle(erfc(z)),colormap=:grays,rasterize=5))
save("erfi.pdf",mesh(graph(bound(norm(erfi(z)),3)),angle(erfi(z)),colormap=:grays,rasterize=5))
save("erfcx.pdf",mesh(graph(bound(norm(erfcx(z)),3)),angle(erfcx(z)),colormap=:grays,rasterize=5))
save("faddeeva.pdf",mesh(graph(bound(norm(faddeeva(z)),3)),angle(faddeeva(z)),colormap=:grays,rasterize=5))
save("dawson.pdf",mesh(graph(bound(norm(dawson(z)),3)),angle(dawson(z)),colormap=:grays,rasterize=5))

t = TensorField(-2π:0.01:2π)
save("logerfc.pdf",graylines(logerfc(t),10))
save("logerfcx.pdf",graylines(logerfcx(t),10))

t = TensorField(0:0.01:1)
save("erfinv.pdf",graylines(bound(erfinv(t),2),10))
save("erfcinv.pdf",graylines(bound(erfcinv(t),2),10))

# besselj

t = TensorField(0:0.01:10)
out = graylines(besselj(0,t),3,rasterize=2)
for i in 1:5
    graylines!(besselj(0.2i,t),3,rasterize=2)
end
save("besselj.pdf",out)

t = TensorField(0:0.01:10)
out = graylines(besselj0(t+0im),3,rasterize=2)
for i in 1:5
    graylines!(besselj0(t+0.2im*i),3,rasterize=2)
end
save("besselj0.pdf",out)

t = TensorField(0:0.01:10)
out = graylines(besselj1(t+0im),3,rasterize=2)
for i in 1:5
    graylines!(besselj1(t+0.2im*i),3,rasterize=2)
end
save("besselj1.pdf",out)

t = TensorField(0:0.01:10)
out = graylines(besseljx(0,t+0im),3,rasterize=2)
for i in 1:5
    graylines!(besseljx(0,t+0.2im*i),3,rasterize=2)
end
save("besseljx.pdf",out)

t = TensorField(0:0.01:10)
out = graylines(sphericalbesselj(1,t+0im),3,rasterize=2)
for i in 1:5
    graylines!(sphericalbesselj(1,t+0.2im*i),3,rasterize=2)
end
save("sphericalbesselj.pdf",out)

t = TensorField(0:0.01:10)
out = graylines(jinc(t+0im),3,rasterize=2)
for i in 1:5
    graylines!(jinc(t+0.2im*i),3,rasterize=2)
end
save("jinc.pdf",out)

# bessely

t = TensorField(0:0.01:10)
out = graylines(bound(bessely(0,t),3),3,rasterize=2)
for i in 1:5
    graylines!(bound(bessely(i,t),3),3,rasterize=2)
end
save("bessely.pdf",out)

t = TensorField(0:0.01:10)
out = graylines(bessely0(t+0.1im),3,rasterize=2)
for i in 1:5
    graylines!(bessely0(t+0.2im*i),3,rasterize=2)
end
save("bessely0.pdf",out)

t = TensorField(0:0.01:10)
out = graylines(bessely1(t+0.1im),3,rasterize=2)
for i in 1:5
    graylines!(bessely1(t+0.2im*i),3,rasterize=2)
end
save("bessely1.pdf",out)

t = TensorField(0:0.01:10)
out = graylines(besselyx(0,t+0.1im),3,rasterize=2)
for i in 1:5
    graylines!(besselyx(0,t+0.2im*i),3,rasterize=2)
end
save("besselyx.pdf",out)

t = TensorField(0:0.001:10)
out = graylines(sphericalbessely(1,t+0.1im),3,rasterize=2)
for i in 1:5
    graylines!(sphericalbessely(1,t+0.2im*i),3,rasterize=2)
end
save("sphericalbessely.pdf",out)

# besselh

t = TensorField(1:0.01:10)
out = graylines(hankelh1(0,t),3,rasterize=2)
for i in 1:5
    graylines!(hankelh1(i,t+i),3,rasterize=2)
end
save("hankelh1.pdf",out)

t = TensorField(1:0.01:10)
out = graylines(hankelh1x(0,t),3,rasterize=2)
for i in 1:5
    graylines!(hankelh1x(i,t+i),3,rasterize=2)
end
save("hankelh1x.pdf",out)

t = TensorField(1:0.01:10)
out = graylines(hankelh2(0,t),3,rasterize=2)
for i in 1:5
    graylines!(hankelh2(i,t+i),3,rasterize=2)
end
save("hankelh2.pdf",out)

t = TensorField(1:0.01:10)
out = graylines(hankelh2x(0,t),3,rasterize=2)
for i in 1:5
    graylines!(hankelh2x(i,t+i),3,rasterize=2)
end
save("hankelh2x.pdf",out)

t = TensorField(0:0.01:5)
out = graylines(besseli(0,t),3,rasterize=2)
for i in 1:5
    graylines!(besseli(0.2i,t),3,rasterize=2)
end
save("besseli.pdf",out)

t = TensorField(0:0.01:5)
out = graylines(besselix(0,t),3,rasterize=2)
for i in 1:5
    graylines!(besselix(0.2i,t),3,rasterize=2)
end
save("besselix.pdf",out)

t = TensorField(0:0.01:5)
out = graylines(besselk(0,t),3,rasterize=2)
for i in 1:5
    graylines!(besselk(0.2i,t),3,rasterize=2)
end
save("besselk.pdf",out)

t = TensorField(0:0.01:5)
out = graylines(besselkx(0,t),3,rasterize=2)
for i in 1:5
    graylines!(besselkx(0.2i,t),3,rasterize=2)
end
save("besselkx.pdf",out)

# Airy

z = complexify(TensorField(ProductSpace(-6:0.01:6,-6:0.01:6)));
for fun ∈ (:airyai,:airyaiprime,:airybi,:airybiprime)
    @eval save($("$fun.pdf"),mesh(graph(bound(norm($fun(z)),3)),angle($fun(z)),colormap=:grays,rasterize=5))
end
for fun ∈ (:airyaix,:airyaiprimex,:airybix,:airybiprimex)
    @eval save($("$fun.pdf"),mesh(graph(norm($fun(z))),angle($fun(z)),colormap=:grays,rasterize=5))
end

# elliptic

t = TensorField(0:0.01:1)
out = graylines(bound(F(π/2,t),3),3,rasterize=2)
for i in 3:7
    graylines!(F(π/i,t),3,rasterize=2)
end
save("EllipticF.pdf",out)

t = TensorField(0:0.01:1)
out = graylines(bound(E(π/2,t),3),3,rasterize=2)
for i in 3:7
    graylines!(E(π/i,t),3,rasterize=2)
end
save("EllipticE.pdf",out)

t = TensorField(0:0.01:1)
out = graylines(bound(Pi(0.5,π/2,t),5),3,rasterize=2)
for i in 0:-1:-3
    graylines!(bound(Pi(i,π/2,t),5),3,rasterize=2)
end
save("EllipticPi.pdf",out)

t = TensorField(0:0.01:1)
out = graylines(bound(real(ellipticZ(π/2,t)),0.5),3,rasterize=2)
for i in 3:7
    graylines!(bound(real(ellipticZ(π/i,t)),0.5),3,rasterize=2)
end
save("EllipticZ.pdf",out)

# jacobi theta

z = complexify(TensorField(ProductSpace(-6:0.01:6,-3:0.01:3)));
save("jtheta1.pdf",mesh(graph(bound(norm(jtheta1(z,0.33)),3)),angle(jtheta1(z,0.33)),colormap=:grays,rasterize=5))
save("jtheta2.pdf",mesh(graph(bound(norm(jtheta2(z,0.33)),3)),angle(jtheta2(z,0.33)),colormap=:grays,rasterize=5))
save("jtheta3.pdf",mesh(graph(bound(norm(jtheta3(z,0.33)),3)),angle(jtheta3(z,0.33)),colormap=:grays,rasterize=5))
save("jtheta4.pdf",mesh(graph(bound(norm(jtheta4(z,0.33)),3)),angle(jtheta4(z,0.33)),colormap=:grays,rasterize=5))

save("ljtheta1.pdf",mesh(graph(bound(norm(ljtheta1(z,0.33)),3)),angle(ljtheta1(z,0.33)),colormap=:grays,rasterize=5))
save("ljtheta2.pdf",mesh(graph(bound(norm(ljtheta2(z,0.33)),3)),angle(ljtheta2(z,0.33)),colormap=:grays,rasterize=5))
save("ljtheta3.pdf",mesh(graph(bound(norm(ljtheta3(z,0.33)),3)),angle(ljtheta3(z,0.33)),colormap=:grays,rasterize=5))
save("ljtheta4.pdf",mesh(graph(bound(norm(ljtheta4(z,0.33)),3)),angle(ljtheta4(z,0.33)),colormap=:grays,rasterize=5))

save("jtheta1dash.pdf",mesh(graph(bound(norm(jtheta1dash(z,0.33)),3)),angle(jtheta1dash(z,0.33)),colormap=:grays,rasterize=5))

# neville theta

z = complexify(TensorField(ProductSpace(-6:0.01:6,-6:0.01:6)));
save("thetaC.pdf",mesh(graph(bound(norm(thetaC(z,m=1.7)),3)),angle(thetaC(z,m=1.7)),colormap=:grays,rasterize=5))
save("thetaD.pdf",mesh(graph(bound(norm(thetaD(z,m=1.7)),3)),angle(thetaD(z,m=1.7)),colormap=:grays,rasterize=5))
save("thetaN.pdf",mesh(graph(bound(norm(thetaN(z,m=1.7)),3)),angle(thetaN(z,m=1.7)),colormap=:grays,rasterize=5))
save("thetaS.pdf",mesh(graph(bound(norm(thetaS(z,m=1.7)),3)),angle(thetaS(z,m=1.7)),colormap=:grays,rasterize=5))

for fun ∈ ("sn","cn","dn","sd","cd","nd","dc","nc","sc","ns","ds","cs")
    save("$fun.pdf",mesh(graph(bound(norm(jellip(fun,z,m=1.7)),3)),angle(jellip(fun,z,m=1.7)),colormap=:grays,rasterize=5))
end

t = TensorField(-10:0.01:10)
out = graylines(Jacobi.am(t,0),3,rasterize=2)
for k in [0.8,0.9,0.95,1]
    graylines!(Jacobi.am(t,k^2),3,rasterize=2)
end
save("am.pdf",out)

t = TensorField(-10:0.01:10)
out = graylines(bound(Jacobi.cn(t,0),3),3,rasterize=2)
for k in [0.8,0.9,0.95,1]
    graylines!(Jacobi.cn(t,k^2),3,rasterize=2)
end
save("cn.pdf",out)

for fun ∈ (:cn,:sn,:dn,:sd,:cd,:nd,:dc,:nc,:sc,:ns,:ds,:cs)
    @eval begin
        t = TensorField(-10:0.01:10)
        out = graylines(bound(Jacobi.$fun(t,0),3),3,rasterize=2)
        for k in [0.8,0.9,0.95,1]
            graylines!(bound(Jacobi.$fun(t,k^2),3),3,rasterize=2)
        end
        save("$($fun).pdf",out)
    end
end

# zeta

z = complexify(TensorField(ProductSpace(-20:0.03:3,-20:0.03:20)));
save("eta.pdf",mesh(graph(bound(norm(eta(z)),10)),angle(eta(z)),colormap=:grays,rasterize=5))
save("zeta.pdf",mesh(graph(bound(norm(zeta(z)),10)),angle(zeta(z)),colormap=:grays,rasterize=5))

z = complexify(TensorField(ProductSpace(-2:0.001:2,0.01:0.01:2)));
save("etaDedekind.pdf",mesh(graph(bound(norm(etaDedekind(z)),10)),angle(etaDedekind(z)),colormap=:grays,rasterize=7))
save("kleinj.pdf",mesh(graph(bound(norm(kleinj(z)),3)),angle(kleinj(z)),colormap=:grays,rasterize=7))
save("lambda.pdf",mesh(graph(bound(norm(lambda(z)),pi/2)),angle(lambda(z)),colormap=:grays,rasterize=7))

# FewSpecialFuntions

t = TensorField(0:0.01:5)
out = graylines(FresnelC(t),3,rasterize=2)
graylines!(FresnelS(t),3,rasterize=2)
save("FresnelCS.pdf",out)
save("FresnelE.pdf",graylines(FresnelE(t),3,rasterize=2))

t = TensorField(LinRange(0,2pi,100))
out = graylines(bound(Clausen(1,t),3),3,rasterize=2)
for i in 2:4
    graylines!(Clausen(i,t),3,rasterize=2)
end
save("Clausen.pdf",out)

t = TensorField(0:0.01:100)
out = graylines(FermiDiracIntegralNorm(-1/2,t),3,rasterize=2)
for i in 1:2:5
    graylines!(FermiDiracIntegralNorm(i/2,t),3,rasterize=2)
end
save("FermiDirac.pdf",out)

t = TensorField(0:0.01:25)
out = graylines(debye_function(1,1,t),3,rasterize=2)
for i in 2:3
    graylines!(debye_function(i,1,t),3,rasterize=2)
end
save("Debye.pdf",out)

t = TensorField(-2.5:0.05:2.5)
out = graylines(bound(V(0.5,t),3),3,rasterize=2)
for i in 1:3
    graylines!(bound(V(0.5+1.5i,t),3),3,rasterize=2)
end
save("V.pdf",out)

t = TensorField(-2.5:0.05:2.5)
out = graylines(U(0.5,t),3,rasterize=2)
for i in 1:3
    graylines!(U(0.5+1.5i,t),3,rasterize=2)
end
save("U.pdf",out)

t = TensorField(0:0.01:10)
out = graylines(MarcumQ(1,0.2,t),3,rasterize=2)
for i in [1.3,2.5,4.7]
    graylines!(MarcumQ(1,i,t),3,rasterize=2)
end
save("MarcumQ.pdf",out)

t = TensorField(0:0.01:10)
out = graylines(dQdb(1,0.2,t),3,rasterize=2)
for i in [1.3,2.5,4.7]
    graylines!(dQdb(1,i,t),3,rasterize=2)
end
save("dQdb.pdf",out)

t = TensorField(0:0.01:25)
out = graylines(FewSpecialFunctions.F(0,0.3,t),3,rasterize=2)
graylines!(FewSpecialFunctions.F(0,-0.3,t),3,rasterize=2)
save("CoulombF.pdf",out)

out = graylines(FewSpecialFunctions.F(1e-5,5,t),3,rasterize=2)
for i in 1:3
    graylines!(FewSpecialFunctions.F(i,5,t),3,rasterize=2)
end
save("CoulombFl.pdf",out)

z = complexify(TensorField(ProductSpace(-3:0.01:2,0:0.01:5)))
save("CoulombFcomplex.pdf",mesh(graph(bound(norm(F(0,z,z)),3)),angle(F(0,z,z)),colormap=:grays,rasterize=7))

t = TensorField(0:0.01:25)
out = graylines(G(0,0.3,t),3,rasterize=2)
graylines!(G(0,-0.3,t),3,rasterize=2)
save("CoulombG.pdf",out)

out = graylines(bound(G(1e-5,5,t),3),3,rasterize=2)
for i in 1:3
    graylines!(bound(G(i,5,t),3),3,rasterize=2)
end
save("CoulombGl.pdf",out)

t = TensorField(0:0.01:25)
out = graylines(H⁺(0,0.3,t),3,rasterize=2)
graylines!(H⁺(0,-0.3,t),3,rasterize=2)
save("CoulombHp.pdf",out)

out = graylines(bound(H⁺(1e-5,5,t),3),3,rasterize=2)
for i in 1:3
    graylines!(bound(H⁺(i,5,t),3),3,rasterize=2)
end
save("CoulombHpl.pdf",out)

z = complexify(TensorField(ProductSpace(-3:0.01:2,0:0.01:5)))
save("CoulombHpcomplex.pdf",mesh(graph(bound(norm(H⁺(0,z,z)),3)),angle(H⁺(0,z,z)),colormap=:grays,rasterize=7))

t = TensorField(0:0.01:25)
out = graylines(H⁻(0,0.3,t),3,rasterize=2)
graylines!(H⁻(0,-0.3,t),3,rasterize=2)
save("CoulombHm.pdf",out)

out = graylines(bound(H⁻(1e-5,5,t),3),3,rasterize=2)
for i in 1:3
    graylines!(bound(H⁻(i,5,t),3),3,rasterize=2)
end
save("CoulombHml.pdf",out)

