using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as der
using DifferentialEquations
using Printf

patm = 101325.0 # Pa
k = 1.4 # gamma
R = 287.0 # J/K
Tatm = 25 + 273.15 # K
prCrit = ((k+1)/2)^(k/(k-1)) # min pressure ratio for choked flow
cp = 1005.0 # J/kg-K
cv = cp / k

pars = @parameters begin
    p0 = 100 * 6894.76 + patm # Pa
    T0 = 25 + 273.15 # K
    D = 3 / 39.37 # m
    L = 4 / 39.37 # m
    de = 10 / 1000 # m
    S = 0.2 / 39.37 # m
    V = 2 / 39.37 # m/s
    hint = 100.0 # W/m2-K
    hext = 5.0 # W/m2-K
    Cpwall = 2 * 500 # J/K, non-specific heat capacity of wall
end

Ae = pi/4 * de^2 # m2, area of inlet/exit
A = pi/4 * D^2 # m2, area of piston
Tstar = T0 * (2/(k+1)) # critical temperature w.r.t. reservoir
tmid = (L - 2*S) / V # time at which spool switches
nCycle = floor(t/tmid) + 1

function isitodd(symb)
    symb % 2 == 1
end

vars = @variables begin
    p(t)[1:2] = [patm + 7; patm + 7]
    T(t)[1:2] = [Tatm + 1e-3; Tatm + 1e-3]
    vol(t)[1:2]
    m(t)[1:2] = [p[1] * vol[1] / (R * T[1]); p[2] * vol[2] / (R * T[2])]
    Twall(t) = Tatm
    Me(t)[1:2] # mach number at inlet/ exit
    Te(t)[1:2] # temperature at inlet/ exit
    rhoe(t)[1:2] # density at inlet/ exit
    Q(t)[1:2] # rate of heat transfer FROM wall into chambers
    x(t)[1:2]
    Ve(t)[1:2] # velocity at inlet/ exit
    pr(t)[1:2] # ratio of pressures across inlet/ exit
    tr(t)[1:2] # temperature ratio across inlet/ exit
    pe(t)[1:2] # pressure at inlet/ exit
end

# define dictionary of initial guesses to replace symbols in initialization system
initvals = Dict([
    Te => [Tstar; Tatm]
    tr => [T0/Tstar; 1]
    ])

eqs = [
    # derivative of ideal gas law (same as equation of mass balance)
    der.(p) ./ p + der.(vol) ./ vol ~ der.(m) ./ m + der.(T) ./ T

    # motion of cylinder (defines volume change)
    x[1] ~ ifelse(isitodd(nCycle),
        S + V*(t - (nCycle-1)*tmid)
        ,
        L - S - V*(t-(nCycle-1)*tmid)
    )
    x[2] ~ L - x[1]
    vol ~ A .* x

    # mach number at inlet/ exit (choked at 1 for pressure ratio more than critical)
    Me ~ ifelse.(pr .> prCrit, 
        [1.0; 1.0]
        ,
        sqrt.(2/(k-1) .* ((pr).^((k-1)/k) .- 1))
    )
    
    # definition of pressure ratio
    pr[1] ~ ifelse(isitodd(nCycle),
        p0/p[1]
        ,
        p[1]/patm
    )
    pr[2] ~ ifelse(isitodd(nCycle),
        p[2]/patm
        ,
        p0/p[2]
    )

    # definition of temperature ratio
    tr[1] ~ ifelse(isitodd(nCycle),
        T0/Te[1]
        ,
        T[1]/Te[1]
    )
    tr[2] ~ ifelse(isitodd(nCycle),
        T[2]/Te[2]
        ,
        T0/Te[2]
    )

    # isentropic flow at the inlet/ exit
    tr[1] ^ (k/(k-1)) ~ ifelse(isitodd(nCycle),
        p0/pe[1]
        ,
        p[1]/pe[1]
        )
    tr[2] ^ (k/(k-1)) ~ ifelse(isitodd(nCycle),
        p[2]/pe[2]
        ,
        p0/pe[2]
        )

    # ideal gas law at the inlet/ exit
    R .* rhoe .* Te ~ pe

    pe[1] ~ 
    ifelse(isitodd(nCycle),
        # e denotes inlet
        ifelse(pr[1] > prCrit,
            p0/prCrit, # critical pressure w.r.t. reservoir
            p[1]) # same as inside
        , # e denotes exit
        ifelse(pr[1] > prCrit,
            p[1]/prCrit, # critical pressure w.r.t. atmosphere
            patm) # same as outside
    )
    pe[2] ~ 
    ifelse(isitodd(nCycle),
        # e denotes inlet
        ifelse(pr[2] > prCrit,
            p[2]/prCrit, # critical pressure w.r.t. atmosphere
            patm) # same as outside
        , # e denotes exit
        ifelse(pr[2] > prCrit,
            p0/prCrit, # critical pressure w.r.t. reservoir
            p[2]) # same as inside
        )
    
    # velocity at inlet/ exit
    Ve ~ Me .* sqrt.(k .* R .* Te)

    # mass flow rate at inlet/ exit
    der(m)[1] ~ ifelse(isitodd(nCycle),
        rhoe[1] * Ae * Ve[1]
        ,
        -rhoe[1] * Ae * Ve[1]
    )
    der(m)[2] ~ ifelse(isitodd(nCycle),
        -rhoe[2] * Ae * Ve[2]
        ,
        rhoe[2] * Ae * Ve[2]
    )

    # rate of heat "addition" from part of wall in contact with top chamber
    Q ~ -pi .* D .* x .* hint .* (T .- Twall)

    # rate of heat "loss" from entire wall (including both chambers, atmosphere, and thermal mass) is 0 (because wall doesn't store heat)
    Q[1] + Q[2] + Cpwall * der(Twall) + hext * pi * D * L * (Twall - Tatm) ~ 0
    # Q ~ -hext * pi * D * L .* (Twall .- Tatm)

    # energy balance (dU/dT = m(cp - h/T) = 0)
    m.*(cv.*T.*der.(p)./p .+ cp.*T.*der.(vol)./vol) ~ der.(m).*(cp.*Te + 1/2 .* Ve.^2) .+ Q
    
]

@mtkcompile sys = System(eqs, t, vars, pars)

prob = ODEProblem(sys, [], (0.0, 10.8); fully_determined = true, guesses = initvals)

sol = solve(prob, Rodas5(), force_dtmin = true, maxiters = Inf)

diameters = (4:1:6) ./ 39.37
all_sols = []
for dia in diameters
    prob_new = remake(prob; p=[D => dia])
    sol_new = solve(prob_new, Rodas5(), force_dtmin = true, maxiters = Inf)
    push!(all_sols, sol_new)
    # Plots.plot!(pfirst, sol_new[t], sol_new[Q[1]+Q[2]], label = "D = $(@sprintf("%.1f",dia*39.37)) in")
end

diameters = pushfirst!(collect(diameters), prob.ps[sys.D])
pushfirst!(all_sols, sol)

import PlotlyJS

# spl = [PlotlyJS.scatter(x=sol[t], y=sol[m[1]*1000], name="m1 (D = 3.0 in)")]
spl = []
for (dia, sol) in zip(diameters, all_sols)
    push!(spl, PlotlyJS.scatter(x=sol[t], y=sol[Ve[1]], name="Ve (D = $(@sprintf("%.1f",dia*39.37)) in)"))
end
spl = [promote(spl...)...]

pl = PlotlyJS.plot(spl, PlotlyJS.Layout(
    title_text="Velocity at inlet/exit"
    , xaxis_title_text="sec"
    , yaxis_title_text="m/s"
    # ,
    # font = PlotlyJS.attr(
    #     family="\"Open Sans\", verdana, arial, sans-serif",
    #     size=24
    # )
))

# spr = [PlotlyJS.scatter(x=sol[t], y=sol[T[1]-273.15], name="T1 (D = 3.0 in)", yaxis="y2")]
# for (dia, sol) in zip(diameters, all_sols)
#     push!(spr, PlotlyJS.scatter(x=sol[t], y=sol[T[1]-273.15], name="T1 (D = $(@sprintf("%.1f",dia*39.37)) in)", yaxis="y2"))
# end
# spr = [PlotlyJS.scatter(x=sol[t], y=sol[hext * pi * D * L * (Twall - Tatm)], name="Qatm (D = 3.0 in)", yaxis="y2")]
# for (dia, sol) in zip(diameters, all_sols)
#     push!(spr, PlotlyJS.scatter(x=sol[t], y=sol[hext * pi * D * L * (Twall - Tatm)], name="Qatm (D = $(@sprintf("%.1f",dia*39.37)) in)", yaxis="y2"))
# end

# ply2 = PlotlyJS.plot(
#     vcat(spl, spr)
#     , PlotlyJS.Layout(
#         title_text="Q"
#         , xaxis_title_text="sec"
#         , yaxis_title_text="Q[1]+Q[2] (W)"
#         , yaxis2=PlotlyJS.attr(
#             title="Qatm (W)"
#             , overlaying="y"
#             , side="right"
#         )
#     )
# )

# PlotlyJS.add_hline!(pl, Tatm, name="Tatm")

PlotlyJS.savefig(pl, "plots/twall_dias.png")

# import Plots; Plots.pythonplot()

# derm2 = getfield(equations(sys)[8], :rhs)
# derp2 = getfield(equations(sys)[9], :rhs)
# derT2 = getfield(equations(sys)[11], :rhs)
# derm1 = getfield(equations(sys)[13], :rhs)
# derp1 = getfield(equations(sys)[14], :rhs)
# derT1 = getfield(equations(sys)[16], :rhs)
# derv1 = getfield(observed(sys)[6], :rhs)
# derv2 = getfield(observed(sys)[15], :rhs)

# pfirst = Plots.plot(sol[t], sol[m[1]*cv*T[1]*derp1/p[1] - derm1*(cp*Te[1] + 1/2 * Ve[1]^2)], label="f(derp/p)")
# Plots.plot!(pfirst, sol[t], sol[m[1]*cp*T[1]*derv1/vol[1]], label="f(derv/v)")
# Plots.plot!(pfirst, sol[t], sol[derm1*(cp*Te[1] + 1/2 * Ve[1]^2)], label="derm*cp*Tt")

# pleft = Plots.plot(sol[t], sol[Q[1]+Q[2]], label="Q1+Q2", legend=:topleft, ylabel="W", color=Plots.theme_palette(:default)[1])
# Plots.plot!(pleft
#     , y_guidefontcolor=Plots.theme_palette(:default)[1]
#     , y_foreground_color_axis=Plots.theme_palette(:default)[1]
#     , y_foreground_color_text=Plots.theme_palette(:default)[1]
#     , y_foreground_color_border=Plots.theme_palette(:default)[1]
#     )
# pright = Plots.twinx();
# Plots.plot!(pright, sol[t], sol[Me[1]], label="Me[1]", xticks=:none, color=Plots.theme_palette(:default)[2])
# Plots.plot!(pright, sol[t], sol[Me[2]], label="Me[2]", xticks=:none, color=Plots.theme_palette(:default)[3])
# Plots.plot!(pright
#     , y_guidefontcolor=Plots.theme_palette(:default)[2]
#     , y_foreground_color_axis=Plots.theme_palette(:default)[2]
#     , y_foreground_color_text=Plots.theme_palette(:default)[2]
#     , y_foreground_color_border=Plots.theme_palette(:default)[2]
#     )

# Plots.plot!(pright, xlabel="\nsec", thickness_scaling=2)
# plot!(pleft, sol_new[t], sol_new[Me[1]], label = "Me (D = 18 in)", color=3)
# plot!(pright, sol_new[t], sol_new[pr[1]], label = "pr (D = 18 in)", xticks=:none, color=4)
