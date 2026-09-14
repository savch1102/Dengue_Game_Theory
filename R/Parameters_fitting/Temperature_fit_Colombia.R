## Libraries ----

library(minpack.lm)
library(rio)
library(tidyverse)

## Municipalities data ----

jerusalen = rio::import("Jerusalen.xlsx") %>% 
  separate(Fecha, into = c("Fecha", "Hour"), sep = " ") %>% 
  select(Fecha, Temperatura) %>% 
  mutate(Temperatura = as.numeric(Temperatura), Fecha = ymd(Fecha),
         day = yday(Fecha)) %>% 
  group_by(day) %>% 
  summarise(mean_temp = mean(Temperatura, na.rm = TRUE))

plot(jerusalen$day, jerusalen$mean_temp)

anapoima = rio::import("Anapoima.xlsx") %>% 
  separate(Fecha, into = c("Fecha", "Hour"), sep = " ") %>% 
  select(Fecha, Temperatura) %>% 
  mutate(Temperatura = as.numeric(Temperatura), Fecha = ymd(Fecha),
         day = yday(Fecha)) %>% 
  group_by(day) %>% 
  summarise(mean_temp = mean(Temperatura, na.rm = TRUE))

plot(anapoima$day, anapoima$mean_temp)

silvania = rio::import("Silvania.xlsx") %>% 
  separate(Fecha, into = c("Fecha", "Hour"), sep = " ") %>% 
  select(Fecha, Temperatura) %>% 
  mutate(Temperatura = as.numeric(Temperatura), Fecha = ymd(Fecha),
         day = yday(Fecha)) %>% 
  group_by(day) %>% 
  summarise(mean_temp = mean(Temperatura, na.rm = TRUE))

plot(silvania$day, silvania$mean_temp)

## Function fit ----

sin_model = function(x, a, b, c, d) {
  a * sin(b * (x - c)) + d
}

## ---- Fit Jerusalen ----

fit_jes = nlsLM(
  mean_temp ~ sin_model(day, a, b, c, d),
  data = jerusalen,
  start = list(
    a = sd(jerusalen$mean_temp),
    b = 2*pi/365,
    c = 100,
    d = 25
  )
)

summary(fit_jes)

x = seq(1, 366)
y1 = 8.037e-01 * sin(3.273e-02 * (x - 1.827e+02)) + 2.866e+01
jerusalen$estimate = y1

ggplot(jerusalen, aes(x = day)) +
  geom_point(aes(y = mean_temp), alpha = 0.3, size = 2) +
  geom_line(aes(y = estimate), color = "blue", size = 0.7, linetype = "dashed") +
  labs(title = "Jerusalen – Mean Temperature (Sinusoidal fit)",
       x = "Day of year", y = "Mean temperature (°C)") +
  theme_bw()

## ---- Fit Anapoima ----

fit_ana = nlsLM(
  mean_temp ~ sin_model(day, a, b, c, d),
  data = anapoima,
  start = list(
    a = sd(anapoima$mean_temp),
    b = 2*pi/365,
    c = 100,
    d = 25
  )
)

summary(fit_ana)

y2 = 6.648e-01 * sin(3.199e-02 * (x - 3.832e+02)) + 2.554e+01
anapoima$estimate = y2

ggplot(anapoima, aes(x = day)) +
  geom_point(aes(y = mean_temp), alpha = 0.3, size = 2) +
  geom_line(aes(y = estimate), color = "blue", size = 0.7, linetype = "dashed") +
  labs(title = "Anapoima – Mean Temperature (Sinusoidal fit)",
       x = "Day of year", y = "Mean temperature (°C)") +
  theme_bw()


## ---- Fit Silvania ----

fit_sil = nlsLM(
  mean_temp ~ sin_model(day, a, b, c, d),
  data = silvania,
  start = list(
    a = (max(silvania$mean_temp) - min(silvania$mean_temp)) / 2,
    b = 2*pi/365,
    c = 180,     
    d = mean(silvania$mean_temp)
  ),
  control = nls.lm.control(maxiter = 500)
)

summary(fit_sil)

y3 = 2.653e-01 * sin(3.151e-02 * (x - 2.062e+02)) + 2.072e+01
silvania$estimate = y3

ggplot(silvania, aes(x = day)) +
  geom_point(aes(y = mean_temp), alpha = 0.3, size = 2) +
  geom_line(aes(y = estimate), color = "blue", size = 0.7, linetype = "dashed") +
  labs(title = "Silvania – Mean Temperature (Sinusoidal fit)",
       x = "Day of year", y = "Mean temperature (°C)") +
  theme_bw()
