## Libraries ----

using CairoMakie
using DataFrames
using DifferentialEquations
using LinearAlgebra
using StatsBase
using Statistics
using SymPy

## Function ----

function model(u, p, time)
    # Extract conditions
    S1, E1, I1, R1, S2, E2, I2, R2, Sv, Ev, Iv = u
    β, α, γ, ζ, prob, ρ, r, μ_A, k, Sp, Tp, k_h, k_v, center_h ,center_v, σ_h, σ_v = p

    # Populations
    I    = I1 + I2                               # Total infected
    Nact = S1 + E1 + R1 + S2 + E2 + R2           # Total active population
    Nh   = S1 + E1 + I1 + R1 + S2 + E2 + I2 + R2 # Total population
    Nv   = Sv + Ev + Iv                          # Total vectors
    py1  = S1 + E1 + R1                          # Total cooperators
    py2  = S2 + E2 + R2                          # Total defectors
    x1   = py1 > 0 ? py1 / Nact : 0.0            # Proportion of cooperators
    x2   = 1 - x1                                # Proportion of defectors

    # FOIs
    λ_h = prob * β * Iv / Nh
    λ_v = prob * β * I / Nh

    # Awareness
    aw_h = σ_h * (1 / (1 + exp(-k_h * ((I/ 37.25) - center_h))))
    aw_v = σ_v * (1 / (1 + exp(-k_v * ((Nv / k) - center_v))))
    σ_t  = aw_h + aw_v

    # Game theory
    A     = Float64[1.0 Sp; (Tp - σ_t) (- σ_t)] #* 1/5# Payoff's matrix
    v     = [x1, x2]                                  # Strategy frequencies
    f_i   = A * v                                     # Fitness
    f1    = f_i[1]                                    # Fitness cooperators
    f2    = f_i[2]                                    # Fitness defectors
    f_bar = v[1] * f1 + v[2] * f2                     # Average fitness
    rate1 = f1 - f_bar                                # Coopertors' rate
    rate2 = f2 - f_bar                                # Defectors' rate

    # Equations
    ## Coopertors
    dS1 = -λ_h * S1 + ζ * R1 + rate1 * S1
    dE1 =  λ_h * S1 - α * E1 + rate1 * E1
    dI1 =    α * E1 - γ * I1
    dR1 =    γ * I1 - ζ * R1 + rate1 * R1

    ## Defectors
    dS2 = -λ_h * S2 + ζ * R2 + rate2 * S2
    dE2 =  λ_h * S2 - α * E2 + rate2 * E2
    dI2 =    α * E2 - γ * I2
    dR2 =    γ * I2 - ζ * R2 + rate2 * R2

    ## Vectors
    dSv =  r * Nv * (1 - (Nv / k))  * exp(-x1) - λ_v * Sv - μ_A * Sv
    dEv =  λ_v * Sv - ρ * Ev - μ_A * Ev
    dIv =  ρ * Ev - μ_A * Iv

    return [dS1, dE1, dI1, dR1, dS2, dE2, dI2, dR2, dSv, dEv, dIv]
end

## Initial conditions and Parameters ----

Nhi = 100

state = [
  Nhi*0.5,   # Susceptible cooperators
  0,         # Exposed cooperators
  0,         # Infected cooperators
  0,         # Recovered cooperators
  Nhi*0.5,   # Susceptible defectors
  0,         # Exposed defectors
  0,         # Infected defectors
  0,         # Recovered defectors
  199,       # Susceptible vectors
  0,         # Exposed vectors
  1          # Infected vectors
]

params = [
    4.5 / 10, # Bititng rate
    1 / 7,    # Human progression rate
    1 / 14,   # Human recovery rate
    0,        # Loss of immunity
    4 / 10,   # Transmission probability
    1 / 7,    # Vector progression rate
    1 / 7,    # Recruitment rate
    1 / 14,   # Adult mosquito mortality rate
    2000,     # Carrying capacity
    -0.5,     # Prisoner's dilemma payoff
    1.5,      # Prisoner's dilemma payoff
    15,       # Slope of human awareness
    15,       # Slope of vector awareness
    0.3,      # Infected human threshold
    0.3,      # Mosquitoes threshold
    0.0,      # Sigma humans
    0.0,      # Sigma vectors
]

tspan = (0, 365*2)

aw_range_b = 0:0.1:1

## Plots ----

f = Figure(size = (900, 350), fontsize = 12)

### A) Amplitude

Z1 = zeros(length(aw_range_b), length(aw_range_b))

for i in eachindex(aw_range_b)
  for j in eachindex(aw_range_b)
    p = copy(params)
    p[end-1] = aw_range_b[i]
    p[end]   = aw_range_b[j]

    res = DataFrame(solve(ODEProblem(model, state, tspan, p); saveat = 0.01))
    rename!(res, [:time, :S1, :E1, :I1, :R1, :S2, :E2, :I2, :R2, :Sv, :Ev, :Iv])

    res.Nact = res.S1 .+ res.E1 .+ res.R1 .+ res.S2 .+ res.E2 .+ res.R2
    res.x1 = (res.S1 .+ res.E1 .+ res.R1) ./ res.Nact
    res.x2 = 1 .- res.x1

    n = length(res.time)
    idx = Int(floor(n*0.5)):n
    amp = maximum(res.x1[idx]) - minimum(res.x1[idx])
    Z1[i, j] = amp
  end
end

ax1 = Axis(f[1, 1],
           title = "A)",
           titlealign = :left,
           xlabel = "Disease awareness (σₕ)",
           ylabel = "Vector awareness (σᵥ)",
           xticks = 0:0.2:1,
           yticks = 0:0.2:1
)

hm1 = heatmap!(ax1, aw_range_b, aw_range_b, Z1, colormap = :matter)

Colorbar(f[1, 2], hm1, label = "Amplitude (x₁)")

### B) Eigen oscilations

@syms S1 E1 I1 R1 S2 E2 I2 R2 Sv Ev Iv
@syms prob_h prob_v β ζ α γ r k μ_A ρ Sp Tp σ_h σ_v k_h k_v center_h center_v

I    = I1 + I2
Nact = S1 + E1 + R1 + S2 + E2 + R2
Nh   = S1 + E1 + I1 + R1 + S2 + E2 + I2 + R2
Nv   = Sv + Ev + Iv
py1  = S1 + E1 + R1
py2  = S2 + E2 + R2
x1   = py1 / Nact
x2   = 1 - x1

λ_h = prob_h * β * Iv / Nh
λ_v = prob_v * β * I / Nh

# I_ref = 37.25
aw_h = σ_h * (1 / (1 + exp(-k_h * ((I / 37.25) - center_h))))
aw_v = σ_v * (1 / (1 + exp(-k_v * ((Nv / k) - center_v))))
σ_t  = aw_h + aw_v

A     = [1.0 Sp; (Tp - σ_t) (-σ_t)]
v_g   = [x1, x2]
f_i   = A * v_g
f1    = f_i[1]
f2    = f_i[2]
f_bar = v_g[1] * f1 + v_g[2] * f2
rate1 = f1 - f_bar
rate2 = f2 - f_bar

dS1 = -λ_h * S1 + ζ * R1 + rate1 * S1
dE1 =  λ_h * S1 - α * E1 + rate1 * E1
dI1 =    α * E1 - γ * I1
dR1 =    γ * I1 - ζ * R1 + rate1 * R1

dS2 = -λ_h * S2 + ζ * R2 + rate2 * S2
dE2 =  λ_h * S2 - α * E2 + rate2 * E2
dI2 =    α * E2 - γ * I2
dR2 =    γ * I2 - ζ * R2 + rate2 * R2

dSv =  r * Nv * (1 - (Nv / k)) * exp(-x1) - λ_v * Sv - μ_A * Sv
dEv =  λ_v * Sv - ρ * Ev - μ_A * Ev
dIv =  ρ * Ev - μ_A * Iv

reduced_eqs = [dS1, dE1, dI1, dR1, dS2, dE2, dI2, dR2, dSv, dEv, dIv]
reduced_vars = [S1, E1, I1, R1, S2, E2, I2, R2, Sv, Ev, Iv]

J_red = simplify(reduced_eqs.jacobian(reduced_vars))

#### Simulation

aw_range = 0.0:0.05:1.0
tspan = (0.0, 365.0 * 5.0)

Z1 = zeros(length(aw_range), length(aw_range)) # Real part (stability)
Z2 = zeros(length(aw_range), length(aw_range)) # Imag part (oscilation frequency)

for i in eachindex(aw_range)
    for j in eachindex(aw_range)
        p = copy(params)
        p[end-1] = aw_range[i] # σ_h
        p[end]   = aw_range[j] # σ_v

        # Solve ODE towards equilibria
        prob_ode = ODEProblem(model, state, tspan, p)
        sol = solve(prob_ode, Vern9(), saveat = 1.0)
        res_jac = sol.u[end] # Final state

        dfe_subs = Dict(
            S1 => res_jac[1], E1 => res_jac[2], I1 => res_jac[3], R1 => res_jac[4],
            S2 => res_jac[5], E2 => res_jac[6], I2 => res_jac[7], R2 => res_jac[8],
            Sv => res_jac[9], Ev => res_jac[10], Iv => res_jac[11],

            β => 0.45, α => 1/7, γ => 1/14, ζ => 0.0,
            prob_h => 0.4, prob_v => 0.4,
            ρ => 1/7, r => 1/7, μ_A => 1/14, k => 2000.0,
            Sp => -0.5, Tp => 1.5,
            k_h => 15.0, k_v => 15.0,
            center_h => 0.3, center_v => 0.3,
            σ_h => aw_range[i],
            σ_v => aw_range[j]
        )

        J_num = Float64.(N(J_red.subs(dfe_subs)))
        eigen_vals = eigvals(J_num)

        # Filter structural eigenvalues near 0 
        non_zero_ev = filter(e -> abs(e) > 1e-8, eigen_vals)

        if !isempty(non_zero_ev)
            # Dominant = higher real part
            dom_idx = argmax(real.(non_zero_ev))
            dom_ev = non_zero_ev[dom_idx]

            Z1[i, j] = real(dom_ev)
            Z2[i, j] = abs(imag(dom_ev))
        else
            Z1[i, j] = 0.0
            Z2[i, j] = 0.0
        end
    end
end

ax2 = Axis(f[1,3], 
           title = "B)", 
           titlealign = :left, 
           xlabel = "Disease awareness (σₕ)",
           ylabel = "Vector awareness (σᵥ)",
           xticks = 0:0.2:1,
           yticks = 0:0.2:1
)

hm2 = heatmap!(ax2, aw_range, aw_range, Z2, colormap = :dense)
Colorbar(f[1,4], hm2, label = "Im(λ)")

rowgap!(f.layout, 20)

f
