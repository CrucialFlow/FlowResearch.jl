# https://github.com/chakravala
# dsp: Antoniou


# Eq. 7.24

wk(t::TensorField,τ,α) = TensorField(base(t),wk.(fiber(t),τ,α))
wk(t,τ,α) = abs(t) ≤ τ ? besseli(0,α*sqrt(1-(t/τ)^2))/besseli(0,α) : 0.0

