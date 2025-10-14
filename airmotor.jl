using ModelingToolkit
using DifferentialEquations
using Symbolics
# using Plots
using ModelingToolkit: t_nounits as t, D_nounits as der

patm = 101325
k = 1.4
R = 287
Tatm = 20 + 273.15 # normal temperature
wrootT_pA = sqrt(k/R*(2/(k+1))^((k+1)/(k-1)))
prCrit = ((k+1)/2)^(k/(k-1))
cp = 1005

pars = @parameters begin
    p0 = 100 * 6894.76 + patm
    T0 = 25 + 273.15
    D = 3 / 39.37
    L = 4 / 39.37
    de = 10 / 1000
    S = 0.2 / 39.37
    V = 2 / 39.37
    hint = 100
    hext = 5
end

vars = @variables begin
    p(t) = patm
    T(t) = T0
    rho(t) = patm / (R * T0)
    mdot(t) = wrootT_pA * p0 / sqrt(T0) * pi/4 * de^2
    vol(t) = pi/4 * D^2 * S
    m(t) = patm / (R * T0) * S * pi/4 * D^2
    Twall(t) = Tatm
    M(t) = 1
    Tin(t)
    rhoin(t)
    phi(t)
    Q(t)
    x(t) = S
    h(t)
end

eqs = [
    # mass balance
    vol * rho / p * der(p) - vol * rho / T * der(T) + rho * der(vol) ~ mdot

    # definition of derivatives
    der(m) ~ mdot
    der(x) ~ V

    vol ~ pi/4 * D^2 * x

    p ~ rho * R * T # ideal gas

    rho ~ m / vol

    # mach number
    ifelse(p0/p > prCrit, 
        M ~ 1
        ,
        M ~ sqrt(2/(k-1) * ((p0/p)^((k-1)/k) - 1))
        )

    # inlet temp
    Tin ~ T0 / (1 + (k-1)/2 * M^2)

    # isentropic flow between reservoir and inlet
    p/p0 ~ (rhoin/(p0/(R*T0)))^k

    # mass flow rate
    mdot ~ rhoin * pi/4 * de^2 * M * sqrt(k * R * Tin)

    # energy flow at inlet
    phi ~ mdot * cp * Tin

    # heat transfer between part of wall in contact with top chamber
    Q ~ pi * D * x * hint * (der(T) - der(Twall))

    # heat transfer between entire wall and atmosphere (same Q because wall doesn't store heat)
    Q ~ hext * pi * D * L * (der(Twall) - Tatm)

    # enthalpy of gas
    h ~ m * cp * T

    # energy balance
    vol * (h/(R*T) - 1) * der(p) + vol * rho * (cp - h / T) * der(T) + rho * h * der(vol) ~ phi + Q
    
]

@named sys = System(eqs, vars, pars)
