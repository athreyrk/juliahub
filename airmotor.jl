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
    p(t) = patm
    T(t)
    vol(t)
    m(t) = p * vol / (R * T)
    Twall(t) = Tatm
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

initvals = Dict([
    Te => Tstar
    tr => T0/Tstar
    ])

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
    -Q ~ hext * pi * D * L * (Twall - Tatm)

    # energy balance (dU/dT = m(cp - h/T) = 0)
    m*(cv*T*der(p)/p + cp*T*der(vol)/vol) ~ der(m)*(cp*Te + 1/2 * Ve^2) + Q
    
]

@mtkcompile sys = System(eqs, t, vars, pars)

prob = ODEProblem(sys, [], (0.0, 3.6); fully_determined = true, guesses = initvals)
sol = solve(prob)

display(plot(sol[t], sol[T-273.15], xlabel="sec", ylabel="degC", label="T"))
display(plot(sol[t], sol[p], xlabel="sec", ylabel="Pa", label="p"))
display(plot(sol[t], sol[m], xlabel="sec", ylabel="kg", label="m"))
display(plot(sol[t], sol[m/vol], xlabel="sec", ylabel="kg/m^3", label="rho"))
display(plot(sol[t], sol[Twall-273.15], xlabel="sec", ylabel="degC", label="Twall"))
display(plot(sol[t], sol[-Q], xlabel="sec", ylabel="W", label="-Q", title="Rate of heat transfer from gas to wall"))