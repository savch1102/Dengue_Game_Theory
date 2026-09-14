## Libraries ----

library(minpack.lm)
library(rio)
library(tidyverse)

## Parameter data ----

br = rio::import("biting_rate.csv")
eip = rio::import("eip_rate.csv")
infection = rio::import("infection_prob.csv")
transm = rio::import("transm_prob.csv")
lifespan = rio::import("lifespan.csv")

## Fit BR ----

gauss_model = function(x, a, b, c) {
  a * exp(-((x - b)^2) / (2 * c^2))
}

fit_br = nlsLM(
  y ~ gauss_model(x, a, b, c),
  data = br,
  start = list(a = max(br$y),    # amplitud inicial
               b = mean(br$x),           # posición del pico
               c = sd(br$x)/2),          # anchura
  control = nls.lm.control(maxiter = 200)
)

summary(fit_br)

x = seq(10, 41)
y =  0.2635 * exp(-((x - 27.8951)^2) / (2 * (5.5970^2)))
fit_df_br = tibble(x, y)

ggplot() +
  geom_point(data = br, aes(x = x, y = y), alpha = 0.3, size = 2)+
  geom_line(data = fit_df_br, aes(x = x, y = y), size = 0.7, linetype = "dashed", color = "blue") +
  labs(title = "Biting rate")+
  theme_bw()

## Fit EIP ----

fit_eip = nlsLM(
  y ~ gauss_model(x, a, b, c),
  data = eip,
  start = list(a = max(eip$y),    # amplitud inicial
               b = mean(eip$x),           # posición del pico
               c = sd(eip$x)/2),          # anchura
  control = nls.lm.control(maxiter = 200)
)

summary(fit_eip)

x = seq(10, 46)
y =  0.1779 * exp(-((x - 35.2827)^2) / (2 * (8.6530^2)))
fit_df_eip = tibble(x, y)

ggplot() +
  geom_point(data = eip, aes(x = x, y = y), alpha = 0.3, size = 2)+
  geom_line(data = fit_df_eip, aes(x = x, y = y), size = 0.7, linetype = "dashed", color = "blue") +
  labs(title = "Extrinsic incubation period")+
  theme_bw()

## Fit Infection ----

fit_infection = nlsLM(
  y ~ gauss_model(x, a, b, c),
  data = infection,
  start = list(a = max(infection$y),
               b = mean(infection$x),
               c = sd(infection$x)/2),
  control = nls.lm.control(maxiter = 200)
)

summary(fit_infection)

x = seq(10, 46)
y = 0.79143 * exp(-((x - 30.83790)^2) / (2 * (5.39891^2)))
fit_df_infection = tibble(x, y)

ggplot() +
  geom_point(data = infection, aes(x = x, y = y), alpha = 0.3, size = 2) +
  geom_line(data = fit_df_infection, aes(x = x, y = y),
            size = 0.7, linetype = "dashed", color = "blue") +
  labs(title = "Infection probability") +
  theme_bw()


## Fit Transmission ----

fit_transm = nlsLM(
  y ~ gauss_model(x, a, b, c),
  data = transm,
  start = list(a = max(transm$y),
               b = mean(transm$x),
               c = sd(transm$x)/2),
  control = nls.lm.control(maxiter = 200)
)

summary(fit_transm)

x = seq(10, 46)
y = 0.7462 * exp(-((x - 30.3610)^2) / (2 * (4.8879^2)))
fit_df_transm = tibble(x, y)

ggplot() +
  geom_point(data = transm, aes(x = x, y = y), alpha = 0.3, size = 2) +
  geom_line(data = fit_df_transm, aes(x = x, y = y),
            size = 0.7, linetype = "dashed", color = "blue") +
  labs(title = "Transmission probability") +
  theme_bw()


## Fit Lifespan ----

fit_lifespan = nlsLM(
  y ~ gauss_model(x, a, b, c),
  data = lifespan,
  start = list(a = max(lifespan$y),
               b = mean(lifespan$x),
               c = sd(lifespan$x)/2),
  control = nls.lm.control(maxiter = 200)
)

summary(fit_lifespan)


x = seq(10, 46)
y = 29.732 * exp(-((x - 22.875)^2) / (2 * (10.119^2)))
fit_df_lifespan = tibble(x, y)

ggplot() +
  geom_point(data = lifespan, aes(x = x, y = y), alpha = 0.3, size = 2) +
  geom_line(data = fit_df_lifespan, aes(x = x, y = y),
            size = 0.7, linetype = "dashed", color = "blue") +
  labs(title = "Adult lifespan") +
  theme_bw()

