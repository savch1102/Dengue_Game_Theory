## Libraries ----

using CairoMakie
using DifferentialEquations
using LinearAlgebra

## Model Definition ----

function model_s2(u, p, t)
    S1, E1, I1, R1, S2, E2, I2, R2, Sv, Ev, Iv = u
    β, α, γ, ζ, prob_h, prob_v, ρ, r, μ_A, k, Sp, Tp, k_h, k_v, center_h, center_v, σ_h, σ_v = p

    I = I1 + I2
    Nact = S1 + E1 + R1 + S2 + E2 + R2
    Nh = S1 + E1 + I1 + R1 + S2 + E2 + I2 + R2
    Nv = Sv + Ev + Iv
    py1 = S1 + E1 + R1
    py2 = S2 + E2 + R2
    x1 = py1 > 0 ? py1 / Nact : 0.0
    x2 = 1.0 - x1

    λ_h = prob_h * β * Iv / Nh
    λ_v = prob_v * β * I / Nh

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

## Parameters and Theoretical Functions ----

r_val = 1 / 7
μ_val = 1 / 14
γ_val = 1 / 14
ρ_val = 1 / 7
β_val = 0.45
ph_val = 0.4
pv_val = 0.4
K_val = 2000.0
Nh_val = 100.0

σ_h_val = 0.0
σ_v_val = 1.0
κ_h_val = 15.0
Tr_h_val = 0.3

c_h = 1.0 / (1.0 + exp(κ_h_val * Tr_h_val))
a_val = 0.5 - σ_h_val * c_h

# Baseline R0 without control (x1 = 0)
Sv0 = K_val * (1.0 - μ_val / r_val)
R0_base = β_val * sqrt((ph_val * pv_val * ρ_val * Sv0) / (Nh_val * γ_val * μ_val * (μ_val + ρ_val)))

function analytical_sigma_v_star(Tr_v, κ_v)
    s_v = 1.0 / (1.0 + exp(-κ_v * ((1.0 - μ_val / r_val) - Tr_v)))
    return (0.5 - σ_h_val * c_h) / s_v
end

function analytical_x1_star(Tr_v, κ_v)
    σ_v_star = analytical_sigma_v_star(Tr_v, κ_v)
    if σ_v_val <= σ_v_star
        return 0.0
    end
    arg_inner = κ_v * (1.0 - Tr_v) + log((σ_v_val - a_val) / a_val)
    arg_outer = (r_val / (κ_v * μ_val)) * arg_inner
    return arg_outer > 1.0 ? log(arg_outer) : 0.0
end

function analytical_Sv_ratio(Tr_v, κ_v)
    x1_s = analytical_x1_star(Tr_v, κ_v)
    return max(0.0, 1.0 - (μ_val / r_val) * exp(x1_s))
end

function analytical_R0_star(Tr_v, κ_v)
    sv_ratio = analytical_Sv_ratio(Tr_v, κ_v)
    return R0_base * sqrt(sv_ratio / (1.0 - μ_val / r_val))
end

## Simulation and Data Collection ----

Tr_v_dense = 0.1:0.005:0.6
Tr_v_points = 0.1:0.05:0.6
kappa_v_list = [5.0, 15.0, 50.0]
colors = [:navy, :crimson, :forestgreen]

state_init = [
    Nh_val * 0.5, 0.0, 0.0, 0.0,
    Nh_val * 0.5, 0.0, 0.0, 0.0,
    Sv0 - 1.0, 0.0, 1.0
]

## Figure Rendering ----

fig_s2 = Figure(size=(850, 700), fontsize=13)

ax1 = Axis(fig_s2[1, 1], title="A) Cooperation Threshold (σ_v*)", xlabel="Perception threshold (Tr_v)", ylabel="σ_v*")
ax2 = Axis(fig_s2[1, 2], title="B) Equilibrium Cooperators (x₁*)", xlabel="Perception threshold (Tr_v)", ylabel="x₁*")
ax3 = Axis(fig_s2[2, 1], title="C) Vector Abundance (S_v*/K)", xlabel="Perception threshold (Tr_v)", ylabel="S_v*/K")
ax4 = Axis(fig_s2[2, 2], title="D) Basic Reproduction Number (R₀*)", xlabel="Perception threshold (Tr_v)", ylabel="R₀*")

for (idx, κ_v_val) in enumerate(kappa_v_list)
    col = colors[idx]
    label_str = "κ_v = $(Int(κ_v_val))"

    # Lines: Analytical
    sigma_v_star_lines = [analytical_sigma_v_star(tr, κ_v_val) for tr in Tr_v_dense]
    x1_star_lines = [analytical_x1_star(tr, κ_v_val) for tr in Tr_v_dense]
    sv_ratio_lines = [analytical_Sv_ratio(tr, κ_v_val) for tr in Tr_v_dense]
    r0_star_lines = [analytical_R0_star(tr, κ_v_val) for tr in Tr_v_dense]

    lines!(ax1, Tr_v_dense, sigma_v_star_lines, color=col, linewidth=2, label=label_str)
    lines!(ax2, Tr_v_dense, x1_star_lines, color=col, linewidth=2, label=label_str)
    lines!(ax3, Tr_v_dense, sv_ratio_lines, color=col, linewidth=2, label=label_str)
    lines!(ax4, Tr_v_dense, r0_star_lines, color=col, linewidth=2, label=label_str)

    # Points: Numerical ODE Simulation
    x1_pts = Float64[]
    sv_pts = Float64[]
    r0_pts = Float64[]

    for tr in Tr_v_points
        p = [
            β_val, 1/7, γ_val, 0.0, ph_val, pv_val, ρ_val, r_val, μ_val, K_val,
            -0.5, 1.5, κ_h_val, κ_v_val, Tr_h_val, tr, σ_h_val, σ_v_val
        ]
        prob_ode = ODEProblem((u, p_in, t) -> model_s2(u, p_in, t), state_init, (0.0, 365.0 * 5.0), p)
        sol = solve(prob_ode, Vern9(); reltol=1e-8, abstol=1e-10)
        u_end = sol.u[end]

        S1_e, E1_e, R1_e = u_end[1], u_end[2], u_end[4]
        S2_e, E2_e, R2_e = u_end[5], u_end[6], u_end[8]
        Sv_e, Ev_e, Iv_e = u_end[9], u_end[10], u_end[11]

        Nact_e = S1_e + E1_e + R1_e + S2_e + E2_e + R2_e
        x1_sim = (S1_e + E1_e + R1_e) / Nact_e
        Nv_sim = Sv_e + Ev_e + Iv_e
        sv_ratio_sim = Nv_sim / K_val
        r0_sim = R0_base * sqrt(max(0.0, sv_ratio_sim) / (1.0 - μ_val / r_val))

        push!(x1_pts, x1_sim)
        push!(sv_pts, sv_ratio_sim)
        push!(r0_pts, r0_sim)
    end

    scatter!(ax2, Tr_v_points, x1_pts, color=col, markersize=8)
    scatter!(ax3, Tr_v_points, sv_pts, color=col, markersize=8)
    scatter!(ax4, Tr_v_points, r0_pts, color=col, markersize=8)
end

hlines!(ax4, [1.0], color=:gray, linestyle=:dash)

axislegend(ax1, position=:rb)

display(fig_s2)
