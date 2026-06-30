# https://github.com/chakravala
# AEM: Lopez, Chapter 8

# Table 8.1

function plotseq(fun,x,n=1:3)
    out = lines(fun(x,n[1]))
    [lines!(fun(x,k)) for k ∈ n[2:end]]
    return out
end

seq1(x,k) = x^k
plotseq(seq1,TensorField(range(0,1,201)),0:2)

seq2(x,k) = sin(k*x)
plotseq(seq2,TensorField(range(0,2pi,201)))

seq3(x,k) = k*x/(1+k^2*x^2)
plotseq(seq3,TensorField(range(0,2,201)))

seq4(x,k) = x/k
plotseq(seq4,TensorField(range(-10,10,201)))

seq5(x,k) = x*exp(-k*x)
plotseq(seq5,TensorField(range(0,1,201)))

seq6(x,k) = x^k/k
plotseq(seq6,TensorField(range(0,1,201)))

seq8(x,k) = k*x*exp(-k*x)
plotseq(seq8,TensorField(range(0,1,201)))

seq11(x,k) = x/(1+k*x^2)
plotseq(seq11,TensorField(range(-1,1,201)))

seq12(x,k) = (1-k*x^2)/(1+k*x^2)^2
plotseq(seq12,TensorField(range(-1,1,201)))





