## Libraries ----

library(minpack.lm)
library(rio)
library(tidyverse)

## Thailand data ----

tai = rio::import("daily_temp_ERA5_provinces.csv")

lampang = tai %>% 
  filter(province == "Lampang") %>% 
  mutate(day = yday(date)) %>% 
  group_by(day) %>% 
  summarise(m_temp = mean(mean_temp))

chumphon = tai %>% 
  filter(province == "Chumphon")%>% 
  mutate(day = yday(date)) %>% 
  group_by(day) %>% 
  summarise(m_temp = mean(mean_temp))

narathiwat = tai %>% 
  filter(province == "Narathiwat")%>% 
  mutate(day = yday(date)) %>% 
  group_by(day) %>% 
  summarise(m_temp = mean(mean_temp))

modelo_lam = lm(m_temp ~ poly(day, 2, raw = TRUE), data = lampang)
summary(modelo_lam)

modelo_chu = lm(m_temp ~ poly(day, 2, raw = TRUE), data = chumphon)
summary(modelo_chu)

modelo_nara = lm(m_temp ~ poly(day, 2, raw = TRUE), data = narathiwat)
summary(modelo_nara)

x = seq(1, 366)
y1 = 22.15 + 0.07619*x - 0.0002185*x^2
lampang$estimate = y1

ggplot(data = lampang, aes(x = day)) +
  geom_point(aes(y = m_temp), alpha = 0.3, size = 2)+
  geom_line(aes(y = estimate), size = 0.7, linetype = "dashed", color = "blue") +
  labs(title = "Lampang - Mean Temperature")+
  theme_bw()

x = seq(1, 366)
y2 = 2.521e+01 + 2.206e-02*x - 6.710e-05*x^2
chumphon$estimate = y2

ggplot(data = chumphon, aes(x = day)) +
  geom_point(aes(y = m_temp), alpha = 0.3, size = 2)+
  geom_line(aes(y = estimate), size = 0.7, linetype = "dashed", color = "blue") +
  labs(title = "Chumphon - Mean Temperature")+
  theme_bw()

x = seq(1, 366)
y3 = 2.433e+01 + 3.101e-02*x - 8.630e-05*x^2
narathiwat$estimate = y3

ggplot(data = narathiwat, aes(x = day)) +
  geom_point(aes(y = m_temp), alpha = 0.3, size = 2)+
  geom_line(aes(y = estimate), size = 0.7, linetype = "dashed", color = "blue") +
  labs(title = "Narathiwat - Mean Temperature")+
  theme_bw()


# Gaussian -------------------------------------------------------------------

## Gaussian model ----
gauss_model = function(x, a, b, c) {
  a * exp(-((x - b)^2) / (2 * c^2))
}

## ---- Fit Lampang ----
fit_lam = nlsLM(
  m_temp ~ gauss_model(day, a, b, c),
  data = lampang,
  start = list(a = max(lampang$m_temp),
               b = 180,
               c = 60)
)

summary(fit_lam)


a_lam = 28.85525
b_lam = 172.86028
c_lam = 243.38327

x = seq(1, 366)
y1 = a_lam * exp(-((x - b_lam)^2) / (2 * c_lam^2))
lampang$estimate = y1

ggplot(lampang, aes(x = day)) +
  geom_point(aes(y = m_temp), alpha = 0.3, size = 2) +
  geom_line(aes(y = estimate), color = "blue", size = 0.7, linetype = "dashed") +
  labs(title = "Lampang – Mean Temperature (Gaussian fit)",
       x = "Day of year", y = "Mean temperature (°C)") +
  theme_bw()

## ---- Fit Chumphon ----
fit_chu = nlsLM(
  m_temp ~ gauss_model(day, a, b, c),
  data = chumphon,
  start = list(a = max(chumphon$m_temp),
               b = 180,
               c = 60)
)

summary(fit_chu)

a_chu = 27.03677
b_chu = 163.56506
c_chu = 439.75645

x = seq(1, 366)
y2 = a_chu * exp(-((x - b_chu)^2) / (2 * c_chu^2))
chumphon$estimate = y2

ggplot(chumphon, aes(x = day)) +
  geom_point(aes(y = m_temp), alpha = 0.3, size = 2) +
  geom_line(aes(y = estimate), color = "blue", size = 0.7, linetype = "dashed") +
  labs(title = "Chumphon – Mean Temperature (Gaussian fit)",
       x = "Day of year", y = "Mean temperature (°C)") +
  theme_bw()


## ---- Fit Narathiwat ----
fit_nara = nlsLM(
  m_temp ~ gauss_model(day, a, b, c),
  data = narathiwat,
  start = list(a = max(narathiwat$m_temp),
               b = 180,
               c = 60)
)

summary(fit_nara)

a_nara = 27.13282
b_nara = 179.26958
c_nara = 386.69663

x = seq(1, 366)
y3 = a_nara * exp(-((x - b_nara)^2) / (2 * c_nara^2))
narathiwat$estimate = y3

ggplot(narathiwat, aes(x = day)) +
  geom_point(aes(y = m_temp), alpha = 0.3, size = 2) +
  geom_line(aes(y = estimate), color = "blue", size = 0.7, linetype = "dashed") +
  labs(title = "Narathiwat – Mean Temperature (Gaussian fit)",
       x = "Day of year", y = "Mean temperature (°C)") +
  theme_bw()




