## Libraries ----

using CSV
using CairoMakie
using DataFrames
using DataInterpolations
using StatsBase
using Statistics

## Data ----

time0 = 0:1:365

temp_v_an = @. 6.648e-01 * sin(3.199e-02 * (time0 - 3.832e+02)) + 2.554e+01
temp_v_na = @. 2.433e+01 + 3.101e-02 * time0 - 8.630e-05 * time0^2

df_an = dropmissing(CSV.read("anap_mean_coh.csv", DataFrame))
df_na = dropmissing(CSV.read("nara_mean_coh.csv", DataFrame))

df_an = combine(groupby(df_an, :period), 
    :mean_coh => mean => :mean,
    :lower_coh => mean => :lower,
    :upper_coh => mean => :upper)

df_na = combine(groupby(df_na, :period),
    :mean_coh => mean => :mean,
    :lower_coh => mean => :lower,
    :upper_coh => mean => :upper)

sort!(df_an, :period)
sort!(df_na, :period)

signals_an = CSV.read("anap_df.csv", DataFrame)
signals_na = CSV.read("nara_df.csv", DataFrame)

## Interpolation ----

period_seq = 1.0:0.1:12.0

func_mean_an  = CubicSpline(df_an.mean, df_an.period)
func_lower_an = CubicSpline(df_an.lower, df_an.period)
func_upper_an = CubicSpline(df_an.upper, df_an.period)

func_mean_na  = CubicSpline(df_na.mean, df_na.period)
func_lower_na = CubicSpline(df_na.lower, df_na.period)
func_upper_na = CubicSpline(df_na.upper, df_na.period)

int_mean_an  = func_mean_an(period_seq)
int_lower_an = clamp.(func_lower_an(period_seq), 0.0, 1.0)
int_upper_an = clamp.(func_upper_an(period_seq), 0.0, 1.0)

int_mean_na  = func_mean_na(period_seq)
int_lower_na = clamp.(func_lower_na(period_seq), 0.0, 1.0)
int_upper_na = clamp.(func_upper_na(period_seq), 0.0, 1.0)

## Helper function para calculo de fase ----

function compute_phase_matrix(df_signals, ranges)
    lags = -365:365
    Z = zeros(length(ranges), length(ranges))
    
    for i in eachindex(ranges)
        for j in eachindex(ranges)
            ex = filter(r -> r.sigma_h == ranges[i] && r.sigma_v == ranges[j], df_signals)
            if nrow(ex) > 0
                ccf = crosscor(ex.x1, ex.temp, lags)
                max_idx = findmax(ccf)[2]
                lag_days = lags[max_idx]
                Z[i, j] = lag_days # (lag_days / 365.0) * 360.0 to degrees
            end
        end
    end
    return Z
end

ranges = 0.0:0.1:1.0
Z_an = compute_phase_matrix(signals_an, ranges)
Z_na = compute_phase_matrix(signals_na, ranges)

## Plot ----

f = Figure(size = (900, 850), fontsize = 12)

colors = cgrad(:inferno, categorical = true, 10)
color_an = colors[2]
color_na = colors[7]

### A) Temperature profiles
ax1 = Axis(f[1,1],
           title = "A)",
           titlealign = :left,
           xlabel = "Time (days)",
           ylabel = "Temperature (°C)"
)

lines!(ax1, time0, temp_v_an, color = color_an, label = "Anapoima")
lines!(ax1, time0, temp_v_na, color = color_na, label = "Narathiwat")
# axislegend(ax1, position = :rt)

### B) Mean coherence
ax2 = Axis(f[1,2],
           title = "B)",
           titlealign = :left,
           xlabel = "Mean coherence",
           ylabel = "Period (months)",
           xticks = 0:0.2:1,
           yticks = 0:2:12,
           limits = ((-0.05, nothing), (nothing, 12.5))
)


band!(ax2, period_seq, int_lower_an, int_upper_an, direction = :y, color = (color_an, 0.08))
lines!(ax2, int_mean_an, period_seq, color = color_an, label = "Anapoima")

band!(ax2, period_seq, int_lower_na, int_upper_na, direction = :y, color = (color_na, 0.08))
lines!(ax2, int_mean_na, period_seq, color = color_na, label = "Narathiwat")

axislegend(ax2, position = :rt)

### C) | D) Phase displacement

ax3 = Axis(f[2,1],
          title = "C)",
          titlealign = :left,
          xlabel = "Disease awareness (σₕ)",
          ylabel = "Vector awareness (σᵥ)",
          xticks = 0:0.2:1,
          yticks = 0:0.2:1
)

hm = heatmap!(ax3, ranges, ranges, Z_an, colormap = :inferno)

ax4 = Axis(f[2,2],
          title = "D)",
          titlealign = :left,
          xlabel = "Disease awareness (σₕ)",
          ylabel = "Vector awareness (σᵥ)",
          xticks = 0:0.2:1,
          yticks = 0:0.2:1
)

hm = heatmap!(ax4, ranges, ranges, Z_na, colormap = :inferno)
Colorbar(f[3,1:2], hm, label = "Phase (days)", vertical = false, flipaxis = false)

colgap!(f.layout, 30)

f
