# https://github.com/chakravala
# Math: Differential Equations and Dynamical Systems, Perko

# 2.4 Example 3

f(x) = Chain(-x[2]/x[3]^2,x[1]/x[3]^2,1)
ic = InitialCondition(f,Chain(0,-1,1/8π),1pi)
stepper = MultistepIntegrator{4}(2^-15)
sol = odesolve(ic,stepper)


