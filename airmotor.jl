using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as der
using DifferentialEquations
using Plots

patm = 101325.0
k = 1.4
R = 287.0
Tatm = 25 + 273.15 # normal temperature
prCrit = ((k+1)/2)^(k/(k-1))
cp = 1005.0
cv = cp / k

pars = @parameters begin
    p0 = 100 * 6894.76 + patm
    T0 = 25 + 273.15
    D = 3 / 39.37
    L = 4 / 39.37
    de = 10 / 1000
    S = 0.2 / 39.37
    V = 2 / 39.37
    hint = 100.0
    hext = 5.0
end

Ae = pi/4 * de^2
A = pi/4 * D^2
Tstar = T0 * (2/(k+1))
rhostar = p0/(R*T0) * (2/(k+1))^(1/(k-1))
tmid = (L - 2*S) / V

vars = @variables begin
    p(t) = patm + 10
    T(t) = Tatm + 1e-3
    vol(t)
    m(t) = p * vol / (R * T)
    Twall(t)# = Tatm
    Me(t)
    Te(t)
    rhoe(t)
    Q(t)# = 0
    x(t)
    Ve(t)
    pr(t)
    tr(t)
    pe(t)
end

function smooth_transition(in, mid, lo, hi; steepness = 100)
    # if in < mid, lo, else, hi, end
    (hi-lo)/2 * tanh(steepness*(in-mid)) + (hi+lo)/2
end

eqs = [
    # derivative of ideal gas law (same as equation of mass balance)
    der(p) / p + der(vol) / vol ~ der(m) / m + der(T) / T

    # motion of cylinder (defines volume change)
    x ~ ifelse(t < tmid,
        S + V*t
        ,
        L - S - V*(t-tmid)
        )
    
    vol ~ A * x

    # mach number at inlet/ exit (choked at 1 for pressure ratio more than critical)
    Me ~ ifelse(pr > prCrit, 
        1.0
        ,
        sqrt(2/(k-1) * ((pr)^((k-1)/k) - 1))
        )
    
    # definition of pressure ratio
    pr ~ ifelse(t < tmid,
        p0/p
        ,
        p/patm
        )

    # definition of temperature ratio
    tr ~ ifelse(t < tmid,
        T0/Te
        ,
        T/Te
        )

    # isentropic flow at the inlet/ exit
    tr ^ (k/(k-1)) ~ ifelse(t < tmid,
        p0/pe
        ,
        p/pe
        )

    # ideal gas law at the inlet/ exit
    R * rhoe * Te ~ pe

    pe ~ smooth_transition(t, tmid,
            ifelse(pr > prCrit,
                p0/prCrit, # critical pressure w.r.t. reservoir
                p) # same as inside
            , # e denotes exit
            ifelse(pr > prCrit,
                p/prCrit, # critical pressure w.r.t. atmosphere
                patm) # same as outside
        )
    
    # velocity at inlet/ exit
    Ve ~ Me * sqrt(k * R * Te)

    # mass flow rate at inlet/ exit
    der(m) ~ ifelse(t < tmid,
        rhoe * Ae * Ve
        ,
        -rhoe * Ae * Ve
        )

    # rate of heat "addition" from part of wall in contact with top chamber
    Q ~ -pi * D * x * hint * (T - Twall)

    # rate of heat "loss" from entire wall to atmosphere (same Q because wall doesn't store heat)
    # -Q ~ hext * pi * D * L * (Twall - Tatm)

    # energy balance (dU/dT = m(cp - h/T) = 0)
    m*(cv*T*der(p)/p + cp*T*der(vol)/vol) ~ der(m)*(cp*Te + 1/2 * Ve^2) + Q
    
]

@named top = System(eqs, t, vars, pars)

btmeqs = [
    # derivative of ideal gas law (same as equation of mass balance)
    der(p) / p + der(vol) / vol ~ der(m) / m + der(T) / T

    # motion of cylinder (defines volume change)
    # x ~ ifelse(t < tmid,
    #     S + V*t
    #     ,
    #     L - S - V*(t-tmid)
    #     )
    
    vol ~ A * x

    # mach number at inlet/ exit (choked at 1 for pressure ratio more than critical)
    Me ~ ifelse(pr > prCrit, 
        1.0
        ,
        sqrt(2/(k-1) * ((pr)^((k-1)/k) - 1))
        )
    
    # definition of pressure ratio
    pr ~ ifelse(t > tmid,
        p0/p
        ,
        p/patm
        )

    # definition of temperature ratio
    tr ~ ifelse(t > tmid,
        T0/Te
        ,
        T/Te
        )

    # isentropic flow at the inlet/ exit
    tr ^ (k/(k-1)) ~ ifelse(t > tmid,
        p0/pe
        ,
        p/pe
        )

    # ideal gas law at the inlet/ exit
    R * rhoe * Te ~ pe

    pe ~ smooth_transition(t, tmid,
            patm,
            ifelse(pr > prCrit,
                p0/prCrit, # critical pressure w.r.t. reservoir
                p) # same as inside
            # , # e denotes exit
            # ifelse(pr > prCrit,
            #     p/prCrit, # critical pressure w.r.t. atmosphere
            #     patm) # same as outside
        )
    
    # velocity at inlet/ exit
    Ve ~ Me * sqrt(k * R * Te)

    # mass flow rate at inlet/ exit
    der(m) ~ ifelse(t > tmid,
        rhoe * Ae * Ve
        ,
        -rhoe * Ae * Ve
        )

    # rate of heat "addition" from part of wall in contact with top chamber
    Q ~ -pi * D * x * hint * (T - Twall)

    # rate of heat "loss" from entire wall to atmosphere (same Q because wall doesn't store heat)
    # -Q ~ hext * pi * D * L * (Twall - Tatm)

    # energy balance (dU/dT = m(cp - h/T) = 0)
    m*(cv*T*der(p)/p + cp*T*der(vol)/vol) ~ der(m)*(cp*Te + 1/2 * Ve^2) + Q
    
]

@named btm = System(btmeqs, t, vars, pars)

fulleqs = [
    top.x + btm.x ~ L
    top.Twall ~ btm.Twall
    top.Q + btm.Q + top.hext * pi * top.D * top.L * (top.Twall - Tatm) ~ 0
]
# fulleqs = [
#     btm.x ~ L -
#     ifelse(t < tmid,
#         S + V*t
#         ,
#         L - S - V*(t-tmid)
#         )
#     # top.Twall ~ btm.Twall
#     btm.Q + btm.hext * pi * btm.D * btm.L * (btm.Twall - Tatm) ~ 0
# ]

@named full = System(fulleqs, t)
@mtkcompile sys = compose(full, top, btm)
# @mtkcompile sys = compose(full, btm)

# @mtkcompile sys = System(eqs, t, vars, pars)
initvals = Dict([
    top.Te => substitute(Tstar, [T0 => top.T0])
    top.tr => T0/Tstar
    btm.Te => Tatm
    btm.tr => 1
    # top.p => patm
    # btm.p => patm
    # top.T => Tatm
    # btm.T => Tatm
    # btm.Twall => Tatm
    # top.m => substitute(patm * A * S / (R * Tatm), [D => top.D, S => top.S])
    # btm.m => substitute(patm * A * (L-S) / (R * Tatm), [D => btm.D, S => btm.S])
    # top.rhoe => substitute(rhostar, [p0 => top.p0, T0 => top.T0])
    # btm.rhoe => patm/ (R * Tatm)
    ])

prob = ODEProblem(sys, [], (0.0, 1.2); fully_determined = true, guesses = initvals)#, warn_cyclic_dependency = true, substitution_limit = 50)
# prob = ODEProblem(sys, [], (0.0, 1.2); fully_determined = false, guesses = initvals)#, warn_cyclic_dependency = true, substitution_limit = 50)
sol = solve(prob)
# alwaysFalse(a...) = false
# sol = solve(prob; adaptive = false, dt = 1e-10, abstol = 1e-4, reltol = 1e-3)
# sol = solve(prob; adaptive = false, dt = 1e-7, unstable_check = alwaysFalse, force_dtmin = true, abstol = 1e-2, reltol = 1e-1)
# sol = solve(prob; dtmin = 1e-3, unstable_check = alwaysFalse, force_dtmin = true, abstol = 1e-2, reltol = 1e-1)

display(plot(sol[t], sol[top.T-273.15], xlabel="sec", ylabel="degC", label="T"))
# display(plot(sol[t], sol[btm.p], xlabel="sec", ylabel="Pa", label="p"))
# display(plot(sol[t], sol[btm.m], xlabel="sec", ylabel="kg", label="m"))
# display(plot(sol[t], sol[btm.m/btm.vol], xlabel="sec", ylabel="kg/m^3", label="rho"))
# display(plot(sol[t], sol[btm.Twall-273.15], xlabel="sec", ylabel="degC", label="Twall"))
# display(plot(sol[t], sol[-btm.Q], xlabel="sec", ylabel="W", label="-Q", title="Rate of heat transfer from gas to wall"))