## Libraries ----

using CSV
using DataFrames
using DifferentialEquations

## Function Anapoima ----

function model_an(u, p, time)
    # Extract conditions
    S1, E1, I1, R1, S2, E2, I2, R2, Sv, Ev, Iv = u
    α, γ, ζ, r, k, Sp, Tp, k_h, k_v, center_h ,center_v, σ_h, σ_v = p

    temp = 6.648e-01 * sin(3.199e-02 * (time - 3.832e+02)) + 2.554e+01
    β    = 0.2635 * exp(-((temp - 27.8951)^2) / (2 * (5.5970^2)))      # Bititng rate
    ρ    = 0.1779 * exp(-((temp - 35.2827)^2) / (2 * (8.6530^2)))      # Vector progression rate
    p_v  = 0.7914 * exp(-((temp - 30.8379)^2) / (2 * (5.3989^2)))      # Probability of infection vectors
    p_h  = 0.7462 * exp(-((temp - 30.3610)^2) / (2 * (4.8879^2)))      # Probability of infection humans
    μ_A  = 1 / (29.732 * exp(-((temp - 22.875)^2) / (2 * (10.119^2)))) # Mortality vectors
    

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
    λ_h = p_h * β * Iv / Nh
    λ_v = p_v * β * I / Nh

    # Awareness
    aw_h = σ_h * (1 / (1 + exp(-k_h * ((I / 37) - center_h))))
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

## Function Narathiwat ----

function model_na(u, p, time)
    # Extract conditions
    S1, E1, I1, R1, S2, E2, I2, R2, Sv, Ev, Iv = u
    α, γ, ζ, r, k, Sp, Tp, k_h, k_v, center_h ,center_v, σ_h, σ_v = p

    t_ciclo = mod(time, 365)
    temp = 2.433e+01 + 3.101e-02 * t_ciclo - 8.630e-05 * t_ciclo^2
    β    = 0.2635 * exp(-((temp - 27.8951)^2) / (2 * (5.5970^2)))      # Bititng rate
    ρ    = 0.1779 * exp(-((temp - 35.2827)^2) / (2 * (8.6530^2)))      # Vector progression rate
    p_v  = 0.7914 * exp(-((temp - 30.8379)^2) / (2 * (5.3989^2)))      # Probability of infection vectors
    p_h  = 0.7462 * exp(-((temp - 30.3610)^2) / (2 * (4.8879^2)))      # Probability of infection humans
    μ_A  = 1 / (29.732 * exp(-((temp - 22.875)^2) / (2 * (10.119^2)))) # Mortality vectors
    

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
    λ_h = p_h * β * Iv / Nh
    λ_v = p_v * β * I / Nh

    # Awareness
    aw_h = σ_h * (1 / (1 + exp(-k_h * ((I / 37.25) - center_h))))
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
    1 / 7,    # Human progression rate
    1 / 14,   # Human recovery rate
    1 / 180,  # Loss of immunity
    1 / 7,    # Recruitment rate
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

tspan = (0, 365*5)
## Signal simulation ----

aw_vector = 0:0.1:1.0

df_an = DataFrame(
    sigma_h = Float64[],
    sigma_v = Float64[],
    time = Float64[],
    temp = Float64[],
    x1 = Float64[]
)

df_na = DataFrame(
    sigma_h = Float64[],
    sigma_v = Float64[],
    time = Float64[],
    temp = Float64[],
    x1 = Float64[]
)

for i in eachindex(aw_vector) 
    for j in eachindex(aw_vector)
        p = copy(params)
        p[end-1] = aw_vector[i]
        p[end] = aw_vector[j]

        res_an = DataFrame(solve(ODEProblem(model_an, state, tspan, p), Vern9(), saveat = 1))
        rename!(res_an, [:time, :S1, :E1, :I1, :R1, :S2, :E2, :I2, :R2, :Sv, :Ev, :Iv])
        t_vec_an  = res_an.time
        temp_v_an = @. 6.648e-01 * sin(3.199e-02 * (t_vec_an - 3.832e+02)) + 2.554e+01
        Nact_an   = res_an.S1 .+ res_an.E1 .+ res_an.R1 .+ res_an.S2 .+ res_an.E2 .+ res_an.R2
        x1_an     = (res_an.S1 .+ res_an.E1 .+ res_an.R1) ./ Nact_an
        x1_an     = @. max(x1_an, 0)

        res_na = DataFrame(solve(ODEProblem(model_na, state, tspan, p), Vern9(), saveat = 1))
        rename!(res_na, [:time, :S1, :E1, :I1, :R1, :S2, :E2, :I2, :R2, :Sv, :Ev, :Iv])
        t_vec_na  = mod.(res_na.time, 365) 
        temp_v_na = @. 2.433e+01 + 3.101e-02 * t_vec_na - 8.630e-05 * t_vec_na^2
        Nact_na   = res_na.S1 .+ res_na.E1 .+ res_na.R1 .+ res_na.S2 .+ res_na.E2 .+ res_na.R2
        x1_na     = (res_na.S1 .+ res_na.E1 .+ res_na.R1) ./ Nact_na
        x1_na     = @. max(x1_na, 0)

        sigma_h_vec = fill(aw_vector[i], length(t_vec_an))
        sigma_v_vec = fill(aw_vector[j], length(t_vec_an))

        new_an = DataFrame(
            sigma_h = sigma_h_vec,
            sigma_v = sigma_v_vec,
            time = t_vec_an,
            temp = temp_v_an,
            x1 = x1_an
        )

        df_an = vcat(df_an, new_an)

        new_na = DataFrame(
            sigma_h = sigma_h_vec,
            sigma_v = sigma_v_vec,
            time = t_vec_na,
            temp = temp_v_na,
            x1 = x1_na
        )

        df_na = vcat(df_na, new_na)
    end
end

df_an.time = df_an.time ./ 365
df_na.time = df_na.time ./ 365

CSV.write("anap_df.csv", df_an)
CSV.write("nara_df.csv", df_na)
