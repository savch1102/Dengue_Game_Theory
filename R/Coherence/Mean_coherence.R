library(psd)
library(rio)
library(tidyverse)

data_an = import("anap_df.csv")
data_na = import("nara_df.csv")

calc_coherence_stats = function(df, conf_level = 0.95) {
  sigma_h_unique = unique(df$sigma_h)
  sigma_v_unique = unique(df$sigma_v)
  
  results_list = list()
  
  for (i in seq_along(sigma_h_unique)) {
    for (j in seq_along(sigma_v_unique)) {
      current = df %>% 
        filter(sigma_h == sigma_h_unique[i], sigma_v == sigma_v_unique[j])
      
      if (nrow(current) < 365) next
      
      my_ts = ts(current[, c("temp", "x1")], frequency = 365)
      spec_res = spec.pgram(my_ts, spans = c(7, 7), plot = FALSE)
      
      freqs = spec_res$freq 
      coh_vals = spec_res$coh
      
      period_months = 12 / freqs
      
      res_months = tibble(
        sigma_h = sigma_h_unique[i],
        sigma_v = sigma_v_unique[j],
        period = 1:12,
        mean_coh = NA_real_,
        lower_coh = NA_real_,
        upper_coh = NA_real_
      )
      
      for (k in 1:12) {
        idx = which(period_months >= (k - 0.5) & period_months < (k + 0.5))
        
        if (length(idx) > 0) {
          vals = coh_vals[idx]
          n = length(vals)
          m = mean(vals, na.rm = TRUE)
          
          if (n > 1 && sd(vals, na.rm = TRUE) > 0) {
            se = sd(vals, na.rm = TRUE) / sqrt(n)
            margin = qt(1 - (1 - conf_level) / 2, df = n - 1) * se
            low = max(0, m - margin)
            high = min(1, m + margin)
          } else {
            low = m
            high = m
          }
          
          res_months$mean_coh[k] = m
          res_months$lower_coh[k] = low
          res_months$upper_coh[k] = high
        }
      }
      
      results_list[[length(results_list) + 1]] = res_months
    }
  }
  
  bind_rows(results_list)
}

df_an_coh = calc_coherence_stats(data_an)
df_na_coh = calc_coherence_stats(data_na)

export(df_an_coh, "anap_mean_coh.csv")
export(df_na_coh, "nara_mean_coh.csv")
