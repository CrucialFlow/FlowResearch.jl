

t = TensorField(-1:0.01:1)

rect(t) = abs(fiber(t))≤1 ? 0.5 : 0.0
rect(t::TensorField) = rect.(t)

d1(t,T) = rect(t/T)/T

d2(t::TensorField,T) = d2.(t,T)
function d2(t,T)
    at = abs(fiber(t))
    at ≤ T ? (T-at)/T^2 : 0.0
end

d3(t,T) = exp(-(t/T)^2)/(sqrt(pi)*T)
d3(t::TensorField,T) = d3.(t,T)

lines(d1(t,0.5))
lines!(d2(t,0.5))
lines!(d3(t,0.5))

dd2(t,T) = abs(fiber(t))≤T ? -sign(fiber(t))/T^2 : 0.0
dd2(t::TensorField,T) = dd2.(t,T)

lines(gradient(d2(t,0.5)))
lines!(dd2(t,0.5))



