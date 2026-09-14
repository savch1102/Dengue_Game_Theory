## Libraries ----

using CairoMakie
using DataFrames
using DifferentialEquations

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

## Plot ----

aw_range = 0:0.5:1

f = Figure(size = (800, 1000), fontsize = 12)

gl_A = f[1, 1:2] = GridLayout()
gl_B = f[2, 1]   = GridLayout()
gl_C = f[2, 2]   = GridLayout()

### A) Cooperators and Defectors (3x3 Grid)

idx = CartesianIndices((3,3))
ax1 = nothing

for id in eachindex(idx) 
    r = id[1]
    c = id[2]

    p = copy(params)
    p[end-1] = aw_range[c]
    p[end]   = aw_range[r]

    res = DataFrame(solve(ODEProblem(model, state, tspan, p), Vern9(), saveat = 0.01))
    rename!(res, [:time, :S1, :E1, :I1, :R1, :S2, :E2, :I2, :R2, :Sv, :Ev, :Iv])

    res.Nact = res.S1 .+ res.E1 .+ res.R1 .+ res.S2 .+ res.E2 .+ res.R2
    res.x1 = (res.S1 .+ res.E1 .+ res.R1) ./ res.Nact
    res.x2 = 1 .- res.x1
    
    ax1 = Axis(gl_A[r+1, c],
               title  = "σₕ = $(aw_range[c]) | σᵥ = $(aw_range[r])",
               xlabel = "Time (days)",
               ylabel = "Proportion",
               #xticks = 0:100:730,
               yticks = 0:0.25:1
    )

    if c != 1
        hideydecorations!(ax1, grid = false)
    end

    if r != 3
        hidexdecorations!(ax1, grid = false)
    end

    lines!(ax1, res.time, res.x1, color = :black, label = "Cooperators")
    lines!(ax1, res.time, res.x2, color = :red, label = "Defectors")
end

Label(gl_A[0, 1], "A)", font = :bold, halign = :left, tellwidth = false)
Legend(gl_A[1, 1:3], ax1, orientation = :horizontal, tellwidth = false)

### B) Peak infection (Percent change)

aw_range_b = 0:0.1:1

Z1 = zeros(length(aw_range_b), length(aw_range_b))
base = 0.0

for i in eachindex(aw_range_b)
    for j in eachindex(aw_range_b)
        p = copy(params)
        p[end-1] = aw_range_b[i]
        p[end]   = aw_range_b[j]

        res = DataFrame(solve(ODEProblem(model, state, tspan, p), Vern9(), saveat = 0.01))
        rename!(res, [:time, :S1, :E1, :I1, :R1, :S2, :E2, :I2, :R2, :Sv, :Ev, :Iv])

        res.Nh = res.S1 .+ res.E1 .+ res.I1 .+ res.R1 .+ res.S2 .+ res.I2 .+ res.E2 .+ res.R2
        res.pI = (res.I1 .+ res.I2) ./ res.Nh
        I_max = maximum(res.pI)

        Z1[i, j] = I_max
        if (aw_range_b[i] == 0.0) && (aw_range_b[j] == 0.0)
            base = I_max
        end
    end
end

Z1 = ((Z1 .- base) ./ base) .* 100

Label(gl_B[1, 1], "B)", font = :bold, halign = :left, tellwidth = false)

ax2 = Axis(gl_B[2, 1],
           xlabel = "Disease awareness (σₕ)",
           ylabel = "Vector awareness (σᵥ)",
           xticks = 0:0.2:1,
           yticks = 0:0.2:1
)

hm1 = heatmap!(ax2, aw_range_b, aw_range_b, Z1, colormap = Reverse(:viridis))

Colorbar(gl_B[2, 2], hm1, label = "Difference to Control (%)")

### C) Cooperators Proportion

Z2 = zeros(length(aw_range_b), length(aw_range_b))

for i in eachindex(aw_range_b)
  for j in eachindex(aw_range_b)
    p = copy(params)
    p[end-1] = aw_range_b[i]
    p[end]   = aw_range_b[j]

    res = DataFrame(solve(ODEProblem(model, state, tspan, p), Vern9(), saveat = 0.01))
    rename!(res, [:time, :S1, :E1, :I1, :R1, :S2, :E2, :I2, :R2, :Sv, :Ev, :Iv])

    res.Nact = res.S1 .+ res.E1 .+ res.R1 .+ res.S2 .+ res.E2 .+ res.R2
    res.x1 = (res.S1 .+ res.E1 .+ res.R1) ./ res.Nact
    res.x2 = 1 .- res.x1

    pcoop = res.x1[end]
    Z2[i, j] = pcoop
  end
end

Label(gl_C[1, 1], "C)", font = :bold, halign = :left, tellwidth = false)

ax3 = Axis(gl_C[2, 1],
           xlabel = "Disease awareness (σₕ)",
           ylabel = "Vector awareness (σᵥ)",
           xticks = 0:0.2:1,
           yticks = 0:0.2:1
)

hm2 = heatmap!(ax3, aw_range_b, aw_range_b, Z2, colormap = :inferno)

Colorbar(gl_C[2, 2], hm2, label = "Proportion of cooperators (x₁)")

rowgap!(f.layout, 20)
rowsize!(f.layout, 1, Relative(2.1/3))

f
