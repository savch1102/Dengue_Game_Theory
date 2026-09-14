## Libraries ----

using CairoMakie
using DataFrames
using DifferentialEquations
using LinearAlgebra
using StatsBase
using Statistics
using SymPy

## System ----

@syms S1 E1 I1 S2 E2 I2 Sv Ev Iv
@syms prob β α γ ρ μ_A k N_h r

infected = [E1, I1, E2, I2, Ev, Iv]

I = I1 + I2

λ_h = prob * β * Iv / N_h
λ_v = prob * β * I / N_h

f_vec = [
    λ_h * S1,  # E1
    0,         # I1
    λ_h * S2,  # E2
    0,         # I2
    λ_v * Sv,  # Ev
    0          # Iv
]

v_vec = [
    α * E1,               # E1
    γ * I1 - α * E1,      # I1
    α * E2,               # E2
    γ * I2 - α * E2,      # I2
    (ρ + μ_A) * Ev,       # Ev
    μ_A * Iv - ρ * Ev     # Iv
]

## Jacobians ----

F = f_vec.jacobian(infected)
V = v_vec.jacobian(infected)

# Substitution = 0
dfe_infected = Dict(E1 => 0, I1 => 0, E2 => 0, I2 => 0, Ev => 0, Iv => 0)
F_dfe = F.subs(dfe_infected)
V_dfe = V.subs(dfe_infected)

K = F_dfe * V_dfe.inv()

K_eigen = K.eigenvals()

R0 = collect(keys(K_eigen))[3]

## Numeric evaluation ----

params = Dict(
    β => 4.5 / 10,    # Biting rate
    α => 1 / 7,       # Human progression rate
    γ => 1 / 14,      # Human recovery rate
    prob => 4 / 10,   # Transmission probability
    ρ => 1 / 7,       # Vector progression rate
    r => 1 / 7,       # Recruitment rate
    μ_A => 1 / 14,    # Adult mosquito mortality rate
    k => 2000,        # Carrying capacity
    N_h => 100        # Total proportion
)

x1_vec = 0:0.001:1
R0_vec = zeros(length(x1_vec))

dSv = r * Sv * (1 - (Sv / k)) * exp(-S1 / N_h) - μ_A * Sv
Sv_eq = solve(dSv, Sv)[2]

for i in eachindex(x1_vec) 
    r_eff = params[r] * exp(-x1_vec[i])
    
    p_i = copy(params)
    
    p_i[S1] = x1_vec[i] * p_i[N_h]
    p_i[S2] = p_i[N_h] - p_i[S1]

    if r_eff > p_i[μ_A]
        Sval = Float64(N(Sv_eq.subs(p_i)))
    else
        Sval = 0.0
    end

    p_i[Sv] = Sval

    K_num = K.subs(p_i)
    K_mat = Float64.(N(K_num))

    autovalores = eigvals(K_mat)
    R0_vec[i] = maximum(abs.(autovalores))
end

## Plot ----

blue_vd = cgrad(:dense, categorical = true, 9)[8]

f = Figure(size = (800, 350), fontsize = 12)
ax1 = Axis(f[1,1],
          title = "A)",
          titlealign = :left,
          xlabel = "Proportion of cooperators (x₁)",
          ylabel = "Basic reproductive number (R₀)",
          xticks = 0:0.2:1,
          yticks = 0:1:7
)

lines!(ax1, x1_vec, R0_vec, color = blue_vd)
band!(ax1, x1_vec, min.(1, R0_vec), R0_vec, color = (blue_vd, 0.1))
hlines!(1, color = :black, linestyle = :dash)

## R0 in the DFE by awareness ----

@syms S1 E1 I1 R1 S2 E2 I2 R2 Sv Ev Iv N_h x1 x2
@syms prob β ζ α γ r k μ_A ρ σ_h σ_v k_h k_v center_h center_v

dSv = r * Sv * (1 - (Sv / k)) * exp(-x1) - μ_A * Sv
Sv_eq = solve(dSv, Sv)[2]

aw_h = σ_h * (1 / (1 + exp(-k_h * ((I / N_h) - center_h))))
aw_v = σ_v * (1 / (1 + exp(-k_v * ((Sv_eq / k) - center_v))))

aw_h = aw_h.subs(Dict(I1 => 0, I2 => 0))

σ_t = aw_h + aw_v

A = [1.0 -0.5; (1.5 - σ_t) (- σ_t)]

a = A[1, 1]
b = A[1, 2]
c = A[2, 1]
d = A[2, 2]

f1 = a * x1 + b * x2
f2 = c * x1 + d * x2

phi = x1 * f1 + x2 * f2

x1_dot = x1 * (f1 - phi)
x2_dot = x2 * (f2 - phi)

x1_dot_mono = x1_dot.subs(Dict(x2 => 1 - x1))
x1_eq = solve(x1_dot_mono, x1)[3]

Sv_eq = Sv_eq.subs(Dict(x1 => x1_eq))

## Solution ----

params = Dict(
    β => 4.5 / 10,    # Biting rate
    α => 1 / 7,       # Human progression rate
    γ => 1 / 14,      # Human recovery rate
    prob => 4 / 10,   # Transmission probability
    ρ => 1 / 7,       # Vector progression rate
    r => 1 / 7,       # Recruitment rate
    μ_A => 1 / 14,    # Adult mosquito mortality rate
    k => 2000,        # Carrying capacity
    N_h => 100,       # Total proportion
    k_h => 15,        # Slope of human awareness
    k_v => 15,        # Slope of vector awareness
    center_h => 0.3,  # Infected human threshold
    center_v => 0.3   # Mosquitoes threshold
)

params[x1] = 1
aw_vec = 0:0.1:1

c_h = 1/(1+exp(params[k_h]*params[center_h]))
D   = 1 - params[center_v]
M   = params[μ_A]/params[r]

function x1_star(σh, σv, p; c_h=c_h, D=D, M=M)
    Y = 0.5 - σh*c_h
    inner_arg = σv/Y - 1
    if inner_arg <= 0
        return σh*c_h + 0.0 > 0.5 ? 1.0 : 0.0
    end
    outer_arg = D + log(inner_arg)/p[k_v]
    if outer_arg <= 0
        return σh*c_h > 0.5 ? 1.0 : 0.0
    end
    x1 = log(outer_arg / M)
    return clamp(x1, 0.0, 1.0)
end

Z_R0 = zeros(length(aw_vec), length(aw_vec))

for i in eachindex(aw_vec), j in eachindex(aw_vec)
    σh_val, σv_val = aw_vec[i], aw_vec[j]
    x1v = x1_star(σh_val, σv_val, params)

    p_i = copy(params)
    p_i[x1] = x1v
    p_i[S1] = x1v * p_i[N_h]
    p_i[S2] = p_i[N_h] - p_i[S1]

    r_eff = params[r]*exp(-x1v)
    Sval = r_eff > params[μ_A] ? params[k]*(1-params[μ_A]*exp(x1v)/params[r]) : 0.0
    p_i[Sv] = Sval

    K_mat = [Float64(N(celda.subs(p_i))) for celda in K]
    Z_R0[i,j] = maximum(abs.(eigvals(K_mat)))
end

## Plot ----

ax2 = Axis(f[1,2],
          title = "B)",
          titlealign = :left,
          xlabel = "Disease awareness (σₕ)",
          ylabel = "Vector awareness (σᵥ)",
          xticks = 0:0.2:1,
          yticks = 0:0.2:1
)


hm = heatmap!(ax2, aw_vec, aw_vec, Z_R0, colormap = :matter)
Colorbar(f[1,3], hm, label = "Basic reproductive number (R₀)")

colgap!(f.layout, 20)
f
