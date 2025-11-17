using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as der
# using OrdinaryDiffEq
using DifferentialEquations
# using Plots

# R = 287.0

pars = @parameters begin
    # p = 101325.0
    # m = 0.0046
    # A = pi/4 * (3/39.37)^2
    S = 0.2 / 39.37
    V = 2 / 39.37
    L = 4 / 39.37
end

tmid = (L - 2*S)/V

vars = @variables begin
    x(t) = S
    # T(t)
end

# vol = A*x

eqs = [
    # p * vol ~ m * R * T
    der(x) ~ V
    ]

@mtkcompile sys = System(eqs, t, vars, pars)

prob_working = ODEProblem(sys, [], (0.0, 10.0)) # does work
prob_nowork = ODEProblem(sys, [], (0.0, tmid)) # doesn't work

sol_working = solve(prob_working)
sol_err = solve(prob_nowork)

# plot(sol, idxs = [x])
