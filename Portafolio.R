# R script: cap4_pipeline.R
# Requeridos:
# install.packages(c("readxl","data.table","xts","PerformanceAnalytics","quadprog","MASS","ggplot2"))

library(readxl)
library(data.table)
library(xts)
library(PerformanceAnalytics)
library(quadprog)
library(MASS)
library(ggplot2)

# 1. Cargar ambos Excel

f1 <- read_xlsx("C:/Users/Steven/Downloads/data_parte1.xlsx",sheet = "Sheet 1", guess_max = 4000)
f2 <- read_xlsx("C:/Users/Steven/Downloads/data_parte2.xlsx",sheet = "Sheet 1", guess_max = 4000)

set.seed(123)

df1 <- as.data.table(f1)
df2 <- as.data.table(f2)

# 2. Unificar (suponiendo mismas columnas) ---

# Estandarizar nombres de columnas clave

setDT(df1)
setDT(df2)

# Normalizar nombres
setnames(df1, "BANCO", "INST_ABRV")

# Asegurar formato de fecha
df1[, FECHA := as.IDate(FECHA)]
df2[, FECHA := as.IDate(FECHA)]

# Verificar claves duplicadas

df1[, .N, by = .(FECHA, INST_ABRV)][N > 1]
df2[, .N, by = .(FECHA, INST_ABRV)][N > 1]

# Unir datos
# LEFT JOIN

df <- merge(
  df1,
  df2,
  by = c("FECHA", "INST_ABRV"),
  all.x = TRUE
)

# INNER JOIN
df <- merge(
  df1,
  df2,
  by = c("FECHA", "INST_ABRV"),
  all = FALSE
)
 
# Resultado

dim(df)
names(df)

# Crear el data frame

setDT(df)

# Crear estructura compatible

df_mk <- df[, .(
  date   = as.Date(FECHA),
  ticker = INST_ABRV,
  price  = ROA.x
)]

# Para verificar
head(df_mk)
str(df_mk)


# 2. Limpieza básica

# Ordenar y eliminar duplicados

# Eliminar NA en columnas clave
df_mk <- df_mk[!is.na(price) & !is.na(date) & !is.na(ticker)]

# Ordenar
setorder(df_mk, ticker, date)

# Verificación
nrow(df_mk)
head(df_mk)

# Crear tabla wide (fecha × banco)

price_wide <- dcast(df_mk, date ~ ticker, value.var = "price")
price_wide

# Convertir a series temporales y retornos

names(price_wide)[1]
names(df_mk)
price_wide <- dcast(df_mk, date ~ ticker, value.var = "price")


dates <- price_wide$date
price_wide[, date := NULL]

price_xts <- xts(price_wide, order.by = dates)

# Retornos

setDT(price_wide)

# Eliminar columna date para cálculos
dates <- price_wide$date
R <- as.matrix(price_wide[, -1, with = FALSE])

# Eliminar activos con demasiados NA
keep <- colMeans(is.na(R)) < 0.2
R <- R[, keep]

# Omitir filas incompletas
R <- na.omit(R)


# 3. Estadística descriptiva

# Media y covarianza (mensual → anual)
mu <- colMeans(R) * 12
Sigma <- cov(R) * 12


#  4. Optimización: frontera eficiente
# Función para obtener pesos de mínima varianza dado retorno objetivo

efficient_portfolio <- function(mu, Sigma, target_return) {
  n <- length(mu)
  Dmat <- 2 * Sigma
  dvec <- rep(0, n)
  # Constraints: A^T w >= b  (quadprog uses: min 1/2 w^T D w - d^T w s.t. A^T w >= b)
  # We need equality constraints: sum(w)=1 and mu^T w = target_return
  Amat <- cbind(rep(1, n), mu)   # columns correspond to constraints
  bvec <- c(1, target_return)
  # quadprog wants A^T w >= b -> we will solve using equality by providing both >= and <= (trick)
  # Use solve.QP with meq = 2 for equality constraints
  res <- solve.QP(Dmat, dvec, Amat, bvec, meq = 2)
  return(res$solution)
}

# Generate frontier
n_points <- 50
rets <- seq(min(mu)*0.5, max(mu)*1.1, length.out = n_points)
weights_mat <- matrix(NA, nrow = length(mu), ncol = n_points)
vars <- numeric(n_points)
for (i in seq_along(rets)) {
  w <- efficient_portfolio(mu, Sigma, rets[i])
  weights_mat[, i] <- w
  vars[i] <- t(w) %*% Sigma %*% w
}
sd_vals <- sqrt(vars)

# 5. Portafolio tangente (max Sharpe)
rf <- 0.02 # ejemplo tasa libre de riesgo anual; ajustar
excess_mu <- mu - rf
w_tangent <- solve.QP(2*Sigma, rep(0,length(mu)), cbind(rep(1,length(mu)), excess_mu),
                      c(1,0), meq = 1)$solution
# Normalizar si fuera necesario
w_tangent <- w_tangent / sum(w_tangent)

# 6. Simulaciones Monte Carlo (multivar normal como ejemplo) ---
n_sim <- 5000
n_steps <- 1  # un periodo para distribución de retornos anuales
sim_returns <- mvrnorm(n = n_sim, mu = mu, Sigma = Sigma)
# Rendimiento del portafolio tangente en cada simulación
port_ret_sim <- as.vector(sim_returns %*% w_tangent)
# Métricas
VaR_95 <- quantile(port_ret_sim, probs = 0.05)
CVaR_95 <- mean(port_ret_sim[port_ret_sim <= VaR_95])

# 7. Visualizaciones ---
# Frontera eficiente
df_front <- data.frame(sd = sd_vals, mean = rets)
ggplot(df_front, aes(x = sd, y = mean)) +
  geom_line() + xlab("Volatilidad (anual)") + ylab("Rendimiento esperado (anual)") +
  ggtitle("Frontera eficiente")

# Histograma del portafolio tangente
qplot(port_ret_sim, bins = 50, xlab = "Rendimiento del portafolio (simulado)", ylab = "Frecuencia")

# 8. Guardar resultados
write.csv(weights_mat, "weights_frontier.csv", row.names = TRUE)
write.csv(data.frame(port_ret_sim = port_ret_sim), "sim_port_tangent.csv")
saveRDS(list(mu = mu, Sigma = Sigma, w_tangent = w_tangent), "params_processed.rds")


# ---- Resultados ----
cat("\n--- RESULTADOS PRINCIPALES ---\n")
print(mu)
print(VaR_95)
print(CVaR_95)

cat("\nPesos portafolio tangente:\n")
print(w_tangent)

# ---- Gráficos ----
p1 <- ggplot(df_front, aes(x = sd, y = mean)) +
  geom_line() +
  xlab("Volatilidad (anual)") +
  ylab("Rendimiento esperado (anual)") +
  ggtitle("Frontera eficiente")

print(p1)

p2 <- qplot(
  port_ret_sim,
  bins = 50,
  xlab = "Rendimiento del portafolio (simulado)",
  ylab = "Frecuencia"
)

print(p2)

# ==============================
# Tabla de rendimientos esperados
# ==============================

# Volatilidad anual por activo
sigma <- sqrt(diag(Sigma))

# Crear tabla resumen
tabla_rendimientos <- data.frame(
  Activo = names(mu),
  Rendimiento_Esperado_Anual = mu,
  Volatilidad_Anual = sigma
)

# Redondear para presentación
tabla_rendimientos$Rendimiento_Esperado_Anual <- round(tabla_rendimientos$Rendimiento_Esperado_Anual, 4)
tabla_rendimientos$Volatilidad_Anual <- round(tabla_rendimientos$Volatilidad_Anual, 4)

# Ver tabla
print(tabla_rendimientos)

write.csv(
  tabla_rendimientos,
  "tabla_rendimientos_esperados.csv",
  row.names = FALSE
)


