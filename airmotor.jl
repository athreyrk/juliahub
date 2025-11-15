using ModelingToolkit
using DifferentialEquations
# using Symbolics
# using Plots
using ModelingToolkit: t_nounits as t, D_nounits as der

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
rhostar = p0/(R*T0) * ((2/k+1)^(1/(k-1)))

vars = @variables begin
    p(t) = patm
    T(t) = Tatm
    rho(t) = patm / (R * Tatm)
    # mdot(t) = wrootT_pA * p0 / sqrt(T0) * pi/4 * de^2
    vol(t)# = A * S
    m(t) = rho * vol
    Twall(t) = Tatm
    M(t)# = 1
    Tin(t) = Tstar
    rhoin(t) = rhostar
    phi(t)
    Q(t)
    x(t) = S
    h(t)
    Vin(t)
end

eqs = [
    # mass balance
    vol * rho / p * der(p) - vol * rho / T * der(T) + rho * der(vol) ~ der(m)#mdot

    # definition of derivatives
    # der(m) ~ mdot
    # der(x) ~ V
    # V ~ der(x)
    x ~ S + V*t

    vol ~ pi/4 * D^2 * x
    # der(vol) ~ A * V

    # p ~ rho * R * T # ideal gas

    rho ~ m / vol

    # mach number
    M ~ ifelse(p0/p > prCrit, 
        1
        ,
        sqrt(2/(k-1) * ((p0/p)^((k-1)/k) - 1))
        )


    # isentropic flow between reservoir and inlet
    Tin ~ T0 / (1 + (k-1)/2 * M^2)

    # ideal gas equation at inlet, same pressure as inside
    rhoin*Tin ~ rho*T

    # p/p0 ~ (rhoin/(p0/(R*T0)))^k

    # mass flow rate
    # mdot ~ rhoin * pi/4 * de^2 * M * sqrt(k * R * Tin)
    Vin ~ M * sqrt(k * R * Tin)
    der(m) ~ rhoin * Ae * Vin

    # energy flow at inlet
    phi ~ cp * Tin + (1/2)*Vin^2

    # heat "addition" from part of wall in contact with top chamber
    Q ~ -pi * D * x * hint * (der(T) - der(Twall))

    # heat "loss" from entire wall to atmosphere (same Q because wall doesn't store heat)
    -Q ~ hext * pi * D * L * (der(Twall) - Tatm)

    # enthalpy of gas
    h ~ cp * T

    # energy balance
    vol * (h/(R*T) - 1) * der(p) + vol * rho * cv * der(T) + rho * h * der(vol) ~ phi + Q
    
]

@named nlsys = ODESystem(eqs, t, vars, pars)
sys = structural_simplify(nlsys)
prob = NonlinearProblem(sys, [], [])
sol = solve(prob)