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

# print.(unknowns(sys),"\n")

prob = ODEProblem(sys, [], (0.0, 1200.0); fully_determined = true, guesses = initvals)

sol = solve(prob, Rodas5(), force_dtmin = true, maxiters = Inf)

import Plots; Plots.pythonplot()

pleft = Plots.plot(sol[t/60], sol[Q[1]+Q[2]], label="from wall into both chambers (combined)", legend=:topleft, ylabel="W", color=Plots.theme_palette(:default)[1])
Plots.plot!(pleft
    , y_guidefontcolor=Plots.theme_palette(:default)[1]
    , y_foreground_color_axis=Plots.theme_palette(:default)[1]
    , y_foreground_color_text=Plots.theme_palette(:default)[1]
    , y_foreground_color_border=Plots.theme_palette(:default)[1]
    )
# pleft = Plots.plot(sol[t], sol[Q[1]+Q[2]], label = "Me (D = 3 in)", legend=:topleft, color=1)
pright = Plots.twinx();
Plots.plot!(pright, sol[t/60], sol[hext*pi*D*L*(Tatm-Twall)], label="from atmosphere into wall", xticks=:none, ylabel="W", color=Plots.theme_palette(:default)[2])
Plots.plot!(pright
    , y_guidefontcolor=Plots.theme_palette(:default)[2]
    , y_foreground_color_axis=Plots.theme_palette(:default)[2]
    , y_foreground_color_text=Plots.theme_palette(:default)[2]
    , y_foreground_color_border=Plots.theme_palette(:default)[2]
    )

Plots.plot!(pright, xlabel="\nmin", title="Rate of heat transfer", thickness_scaling=2)
# plot!(pleft, sol_new[t], sol_new[Me[1]], label = "Me (D = 18 in)", color=3)
# plot!(pright, sol_new[t], sol_new[pr[1]], label = "pr (D = 18 in)", xticks=:none, color=4)


# import PlotlyJS

# pl = plot(
#     [
#         scatter(x=sol[t], y=sol[Q[1]+Q[2]], name="Q[1]+Q[2]")
#     ]
#     , Layout(
#         title_text="rate of heat transfer from wall into chambers",
#         xaxis_title_text="sec",
#         yaxis_title_text="W"
#         # ,
#         # yaxis2=attr(
#         #     title="pr",
#         #     overlaying="y",
#         #     side="right"
#         # )
#     )
# )

# ptemp = plot(
#     [
#         scatter(x=sol[t], y=sol[getfield(equations(sys)[7], :rhs)], name="der(Twall)")
#     ]
#     , Layout(
#         title_text="rate of wall temperature",
#         xaxis_title_text="sec",
#         yaxis_title_text="K/s"
#     )
# )

# plotatm = plot(
#     scatter(x=sol[t], y=sol[hext * pi * D * L * (Twall - Tatm)], name="qatm")
#     , Layout(
#         title_text="rate of heat transfer from wall to atm",
#         xaxis_title_text="sec",
#         yaxis_title_text="K/s"
#     )
# )

# prob_new = remake(prob; p=[D => 18/39.37])
# sol_new = solve(prob_new, Rodas5(), force_dtmin = true)#, maxiters = Inf)

# pl = plot(
#     [
#         scatter(x=sol[t], y=sol[Me[1]], name="Me (D=3)"),
#         scatter(x=sol_new[t], y=sol_new[Me[1]], name="Me (D=18)"),
#         scatter(x=sol[t], y=sol[pr[1]], name="pr (D=3)", yaxis="y2"),
#         scatter(x=sol_new[t], y=sol_new[pr[1]], name="pr (D=18)", yaxis="y2")
#     ],
#     Layout(
#         title_text="Me and pr",
#         xaxis_title_text="sec",
#         yaxis_title_text="Me",
#         yaxis2=attr(
#             title="pr",
#             overlaying="y",
#             side="right"
#         )
#     )
# )

pfirst = Plots.plot(sol[t/60], sol[Twall-273.15], label="D = 3.0 in")

diameters = [8 9 13 18] ./ 39.37

for dia in diameters
    prob_new = remake(prob; p=[D => dia])
    sol_new = solve(prob_new, Rodas5(), force_dtmin = true, maxiters = Inf)
    Plots.plot!(pfirst, sol_new[t/60], sol_new[Twall-273.15], label = "D = $(@sprintf("%.1f",dia*39.37)) in")
end

Plots.plot!(pfirst, xlabel="min", ylabel="degC", title="wall temperature", thickness_scaling=2)
# savefig(plt, "twall_dias_no_freezing.png")
