using ModelingToolkit
using DifferentialEquations
# using Symbolics
# using Plots
using ModelingToolkit: t_nounits as t, D_nounits as der
# using IfElse

patm = 101325.0
k = 1.4
R = 287.0
Tatm = 20 + 273.15 # normal temperature
wrootT_pA = sqrt(k/R*(2/(k+1))^((k+1)/(k-1)))
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
    T(t)# = Tatm
    # rho(t) = patm / (R * T0)
    # mdot(t) = wrootT_pA * p0 / sqrt(T0) * pi/4 * de^2
    vol(t)# = A * x
    m(t) = p * vol / (R * T)
    Twall(t) = Tatm
    Me(t)# = 1
    Te(t)# = Tstar
    rhoe(t)# = rhostar
    # phi(t)
    Q(t)
    x(t)# = S
    # h(t)
    Ve(t)# = sqrt(k*R*Tstar)
    pr(t)# = p0/p
    tr(t)# = T0/Te
    pe(t)# = p0/prCrit
end
# Me = 1.0
# pr = p0/p
# tr = T0/Te
# vol = A*(L-S)
initvals = Dict([
    # p => patm
    # T => T0
    # vol => A*S
    # m => patm * A*S/(R*T0)
    # Twall => Tatm
    # Me => 1.0
    Te => Tstar
    # rhoe => rhostar
    # x => S
    # Ve => sqrt(k*R*Tstar)
    # pr => p0/patm
    tr => T0/Tstar
])
# @register_symbolic(t < tmid)

eqs = [
    # derivative of ideal gas law (same as equation of mass balance)
    der(p) / p + der(vol) / vol ~ der(m) / m + der(T) / T
    # der(p) / p ~ der(m) / m + der(T) / T
    # vol * rho / p * der(p) - vol * rho / T * der(T) + rho * der(vol) ~ der(m)#mdot

    # p * vol ~ m * R * T
    # motion of cylinder (defines volume change)
    x ~ ifelse(t < tmid, 
        S + V*t
        ,
        L - S - V*t
        )
    # x ~ S + V*t
    vol ~ A * x

    # mach number at inlet/ exit (choked at 1 for pressure ratio more than critical)
    Me ~ ifelse(pr > prCrit, 
        1.0
        ,
        sqrt(2/(k-1) * ((pr)^((k-1)/k) - 1))
        )
    # Me ~ 1.0
    
    # definition of pressure ratio
    pr ~ ifelse(t < tmid,
        p0/p
        ,
        p/patm
        )
    # pr ~ p0/p

    # definition of temperature ratio
    tr ~ ifelse(t < tmid,
        T0/Te
        ,
        T/Te
        )
    # tr ~ T0/Te

    # isentropic flow at the inlet/ exit
    tr ^ (k/(k-1)) ~ ifelse(t < tmid,
        p0/pe
        ,
        p/pe
        )
    
    # ideal gas law at the inlet/ exit
    R * rhoe * Te ~ pe

    pe ~
        ifelse(t < tmid, # e denotes inlet
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
    der(m) ~ rhoe * Ae * Ve

    # energy flow at inlet
    # phi ~ cp * Te + (1/2)*Ve^2

    # heat "addition" from part of wall in contact with top chamber
    Q ~ -pi * D * x * hint * (T - Twall)

    # heat "loss" from entire wall to atmosphere (same Q because wall doesn't store heat)
    -Q ~ hext * pi * D * L * (Twall - Tatm)

    # enthalpy of gas
    # h ~ cp * T

    # energy balance
    # vol * (h/(R*T) - 1) * der(p) + vol * rho * cv * der(T) + rho * h * der(vol) ~ cp*Te + (1/2)*Ve^2 + Q/der(m)
    m*(cv*T*der(p)/p + cv*der(T) + cp*T*der(vol)/vol) ~ der(m)*(cp*Te + 1/2 * Ve^2) + der(Q)
    
]

@mtkcompile sys = System(eqs, t, vars, pars)

# @named nlsys = ODESystem(eqs, t, vars, pars)
# sys = structural_simplify(nlsys)
prob = ODEProblem(sys, [], (0.0, 1.8); fully_determined = true, guesses = initvals)
sol = solve(prob)
# traced_sys = modelingtoolkitize(prob)
# compiled_sys = mtkcompile(dae_index_lowering(traced_sys))
# new_prob = ODEProblem(compiled_sys, Pair[], (0, tmid), [])
# sol = solve(new_prob)
# plot(sol.t, sol[T])