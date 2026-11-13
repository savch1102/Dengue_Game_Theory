## Libraries ----

using CairoMakie
using DataFrames
using DifferentialEquations
using LinearAlgebra
using Statistics

## Functions ----

function model(u, p, t)
    S1, E1, I1, R1, S2, E2, I2, R2, Sv, Ev, Iv = u
    β, α, γ, ζ, prob, ρ, r, μ_A, k, Sp, Tp, k_h, k_v, center_h, center_v, σ_h, σ_v = p

    I = I1 + I2
    Nact = S1 + E1 + R1 + S2 + E2 + R2
    Nh = S1 + E1 + I1 + R1 + S2 + E2 + I2 + R2
    Nv = Sv + Ev + Iv
    py1 = S1 + E1 + R1
    py2 = S2 + E2 + R2
    x1 = py1 > 0 ? py1 / Nact : 0.0
    x2 = 1.0 - x1

    λ_h = prob * β * Iv / Nh
    λ_v = prob * β * I / Nh

    aw_h = σ_h * (1.0 / (1.0 + exp(-k_h * ((I / 37.25) - center_h))))
    aw_v = σ_v * (1.0 / (1.0 + exp(-k_v * ((Nv / k) - center_v))))
    σ_t = aw_h + aw_v

    A = Float64[1.0 Sp; (Tp - σ_t) (-σ_t)]
    v = [x1, x2]
    f_i = A * v
    f1 = f_i[1]
    f2 = f_i[2]
    f_bar = v[1] * f1 + v[2] * f2
    rate1 = f1 - f_bar
    rate2 = f2 - f_bar

    dS1 = -λ_h * S1 + ζ * R1 + rate1 * S1
    dE1 = λ_h * S1 - α * E1 + rate1 * E1
    dI1 = α * E1 - γ * I1
    dR1 = γ * I1 - ζ * R1 + rate1 * R1

    dS2 = -λ_h * S2 + ζ * R2 + rate2 * S2
    dE2 = λ_h * S2 - α * E2 + rate2 * E2
    dI2 = α * E2 - γ * I2
    dR2 = γ * I2 - ζ * R2 + rate2 * R2

    dSv = r * Nv * (1.0 - (Nv / k)) * exp(-x1) - λ_v * Sv - μ_A * Sv
    dEv = λ_v * Sv - (ρ + μ_A) * Ev
    dIv = ρ * Ev - μ_A * Iv

    return [dS1, dE1, dI1, dR1, dS2, dE2, dI2, dR2, dSv, dEv, dIv]
end

function compute_jacobian_fd(u_ss, p; eps=1e-7)
    n = length(u_ss)
    J = zeros(Float64, n, n)
    for j in 1:n
        u_plus = copy(u_ss)
        u_minus = copy(u_ss)
        u_plus[j] += eps
        u_minus[j] -= eps
        f_plus = model(u_plus, p, 0.0)
        f_minus = model(u_minus, p, 0.0)
        J[:, j] = (f_plus .- f_minus) ./ (2.0 * eps)
    end
    return J
end

## Parameters and Setup ----

Nhi = 100.0

state = [
    Nhi * 0.5, 0.0, 0.0, 0.0,
    Nhi * 0.5, 0.0, 0.0, 0.0,
    199.0, 0.0, 1.0
]

base_params = [
    4.5 / 10, 1 / 7, 1 / 14, 0.0, 4 / 10,
    1 / 7, 1 / 7, 1 / 14, 2000.0,
    -0.5, 1.5, 15.0, 15.0, 0.3, 0.3, 0.0, 0.0
]

aw_range = 0.0:0.02:1.0
tspan = (0.0, 365.0 * 5.0)

re_lambda_matrix = zeros(length(aw_range), length(aw_range))
tau_damp_matrix = zeros(length(aw_range), length(aw_range))

## Simulation Sweep ----

for (i, σ_h_val) in enumerate(aw_range)
    for (j, σ_v_val) in enumerate(aw_range)
        p = copy(base_params)
        p[end-1] = σ_h_val
        p[end] = σ_v_val

        prob_ode = ODEProblem((u, p_in, t) -> model(u, p_in, t), state, tspan, p)
        sol = solve(prob_ode, Vern9(); reltol=1e-8, abstol=1e-10)
        u_ss = sol.u[end]

        J_num = compute_jacobian_fd(u_ss, p)
        evals = eigvals(J_num)

        # Filter structural conservation zero eigenvalues
        non_zero_evals = filter(e -> abs(e) > 1e-5, evals)

        if !isempty(non_zero_evals)
            dom_idx = argmax(real.(non_zero_evals))
            dom_ev = non_zero_evals[dom_idx]
            
            re_val = real(dom_ev)
            re_lambda_matrix[i, j] = re_val
            tau_damp_matrix[i, j] = re_val < 0.0 ? -1.0 / re_val : NaN
        else
            re_lambda_matrix[i, j] = NaN
            tau_damp_matrix[i, j] = NaN
        end
    end
end

## Figure Rendering ----

fig_s1 = Figure(size=(900, 380), fontsize=13)

ax_a = Axis(
    fig_s1[1, 1],
    title="A) Dominant Real Eigenvalue Re(λ_dom)",
    titlealign=:left,
    xlabel="Disease awareness (σₕ)",
    ylabel="Vector awareness (σᵥ)",
    xticks=0:0.2:1,
    yticks=0:0.2:1
)

hm_a = heatmap!(ax_a, aw_range, aw_range, re_lambda_matrix, colormap=:viridis)
Colorbar(fig_s1[1, 2], hm_a, label="Re(λ_dom)")

ax_b = Axis(
    fig_s1[1, 3],
    title="B) Damping Time τ_damp (days)",
    titlealign=:left,
    xlabel="Disease awareness (σₕ)",
    ylabel="Vector awareness (σᵥ)",
    xticks=0:0.2:1,
    yticks=0:0.2:1
)

hm_b = heatmap!(ax_b, aw_range, aw_range, tau_damp_matrix, colormap=:inferno)
Colorbar(fig_s1[1, 4], hm_b, label="τ_damp (days)")

colgap!(fig_s1.layout, 15)

display(fig_s1)
