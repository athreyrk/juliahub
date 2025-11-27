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
    p(t)[1:2] = [patm + 7; patm + 7]
    T(t)[1:2] = [Tatm + 1e-3; Tatm + 1e-3]
    vol(t)[1:2]
    m(t)[1:2] = [p[1] * vol[1] / (R * T[1]); p[2] * vol[2] / (R * T[2])]
    Twall(t)#[1:2]# = Tatm
    Me(t)[1:2]
    Te(t)[1:2]
    rhoe(t)[1:2]
    Q(t)[1:2]# = 0
    x(t)[1:2]
    Ve(t)[1:2]
    pr(t)[1:2]
    tr(t)[1:2]
    pe(t)[1:2]
end

initvals = Dict([
    Te => [Tstar;Tatm]
    tr => [T0/Tstar;1]
    # T[2] => Tatm
    ])

function smooth_transition(in, mid, lo, hi; steepness = 100)
    # if in < mid, lo, else, hi, end
    (hi-lo)/2 * tanh(steepness*(in-mid)) + (hi+lo)/2
end

eqs = [
    # derivative of ideal gas law (same as equation of mass balance)
    der.(p) ./ p + der.(vol) ./ vol ~ der.(m) ./ m + der.(T) ./ T
    # der(p)[1] / p[1] + der(vol)[1] / vol[1] ~ der(m)[1] / m[1] + der(T)[1] / T[1]
    # if t < tmid
    #     der(p)[2] ~ 0
    #     der(T)[2] ~ 0
    # else
    #     der(p)[2] / p[2] + der(vol)[2] / vol[2] ~ der(m)[2] / m[2] + der(T)[2] / T[2]
    # end
    # 0 ~ ifelse(t < tmid,
    #     der(T)[2],
    #     der(p)[2] / p[2] + der(vol)[2] / vol[2] - der(m)[2] / m[2] - der(T)[2] / T[2]
    # )

    # motion of cylinder (defines volume change)
    x[1] ~ ifelse(t < tmid,
        S + V*t
        ,
        L - S - V*(t-tmid)
    )
    x[2] ~ L - x[1]
    vol ~ A .* x

    # mach number at inlet/ exit (choked at 1 for pressure ratio more than critical)
    Me ~ ifelse.(pr .> prCrit, 
        [1.0;1.0]
        ,
        sqrt.(2/(k-1) .* ((pr).^((k-1)/k) .- 1))
    )
    
    # definition of pressure ratio
    pr[1] ~ ifelse(t < tmid,
        p0/p[1]
        ,
        p[1]/patm
    )
    pr[2] ~ ifelse(t < tmid,
        p[2]/patm
        ,
        p0/p[2]
    )

    # definition of temperature ratio
    tr[1] ~ ifelse(t < tmid,
        T0/Te[1]
        ,
        T[1]/Te[1]
    )
    tr[2] ~ ifelse(t < tmid,
        T[2]/Te[2]
        ,
        T0/Te[2]
    )

    # isentropic flow at the inlet/ exit
    tr[1] ^ (k/(k-1)) ~ ifelse(t < tmid,
        p0/pe[1]
        ,
        p[1]/pe[1]
        )
    tr[2] ^ (k/(k-1)) ~ ifelse(t < tmid,
        p[2]/pe[2]
        ,
        p0/pe[2]
        )

    # ideal gas law at the inlet/ exit
    R .* rhoe .* Te ~ pe

    pe[1] ~ 
    ifelse(t < tmid,
    # smooth_transition(t, tmid,
        ifelse(pr[1] > prCrit,
            p0/prCrit, # critical pressure w.r.t. reservoir
            p[1]) # same as inside
        , # e denotes exit
        ifelse(pr[1] > prCrit,
            p[1]/prCrit, # critical pressure w.r.t. atmosphere
            patm) # same as outside
    )
    pe[2] ~ 
    ifelse(t < tmid,
    # smooth_transition(t, tmid,
            patm # same as outside
        , # e denotes exit
        ifelse(pr[2] > prCrit,
            p0/prCrit, # critical pressure w.r.t. reservoir
            p[2]) # same as inside
        )
    
    # velocity at inlet/ exit
    Ve ~ Me .* sqrt.(k .* R .* Te)

    # mass flow rate at inlet/ exit
    der(m)[1] ~ ifelse(t < tmid,
        rhoe[1] * Ae * Ve[1]
        ,
        -rhoe[1] * Ae * Ve[1]
    )
    der(m)[2] ~ ifelse(t < tmid,
        -rhoe[2] * Ae * Ve[2]
        ,
        rhoe[2] * Ae * Ve[2]
    )

    # rate of heat "addition" from part of wall in contact with top chamber
    Q ~ -pi .* D .* x .* hint .* (T .- Twall)

    # rate of heat "loss" from entire wall (including both chambers and atmosphere) is 0 (because wall doesn't store heat)
    Q[1] + Q[2] + hext * pi * D * L * (Twall - Tatm) ~ 0
    # Q ~ -hext * pi * D * L .* (Twall .- Tatm)

    # energy balance (dU/dT = m(cp - h/T) = 0)
    m.*(cv.*T.*der.(p)./p + cp.*T.*der.(vol)./vol) ~ der.(m).*(cp.*Te + 1/2 .* Ve.^2) .+ Q
    # m[1]*(cv*T[1]*der(p)[1]/p[1] + cp*T[1]*der(vol)[1]/vol[1]) ~ der(m)[1]*(cp*Te[1] + 1/2 * Ve[1]^2) + Q[1]
    
]

@mtkcompile sys = System(eqs, t, vars, pars)

prob = ODEProblem(sys, [], (0.0, 3.6); fully_determined = true, guesses = initvals)#, substitution_limit = 5, warn_cyclic_dependency = true)
# alwaysFalse(a...) = false
sol = solve(prob, Rosenbrock23())#, unstable_check = alwaysFalse, force_dtmin = true)

# display(plot(sol[t], sol[T[1]-273.15], xlabel="sec", ylabel="degC", label="T"))
# display(plot(sol[t], sol[p], xlabel="sec", ylabel="Pa", label="p"))
# display(plot(sol[t], sol[m], xlabel="sec", ylabel="kg", label="m"))
# display(plot(sol[t], sol[m/vol], xlabel="sec", ylabel="kg/m^3", label="rho"))
# display(plot(sol[t], sol[Twall-273.15], xlabel="sec", ylabel="degC", label="Twall"))
# display(plot(sol[t], sol[-Q], xlabel="sec", ylabel="W", label="-Q", title="Rate of heat transfer from gas to wall"))