############################################################
# cap5_pipeline_clean.R
# Optimización de Portafolios y Simulación de Monte Carlo
############################################################

# Paquetes requeridos
# install.packages(c("readxl","data.table","quadprog","MASS","ggplot2"))

library(readxl)
library(data.table)
library(quadprog)
library(MASS)
library(ggplot2)

set.seed(123)

############################################################
# 1. Carga y preprocesamiento de los datos
############################################################

# Cargar bases de datos
df1 <- as.data.table(
  read_xlsx("C:/Users/Steven/Downloads/data_parte1.xlsx", guess_max = 4000)
)
df2 <- as.data.table(
  read_xlsx("C:/Users/Steven/Downloads/data_parte2.xlsx", guess_max = 4000)
)

# Estandarizar variables clave
setnames(df1, "BANCO", "INST_ABRV")
df1[, FECHA := as.IDate(FECHA)]
df2[, FECHA := as.IDate(FECHA)]

# Unir bases de datos (inner join para consistencia)
df <- merge(df1, df2, by = c("FECHA", "INST_ABRV"), all = FALSE)

# Construcción del conjunto de trabajo (proxy de rendimientos)
df_mk <- df[, .(
  date   = as.Date(FECHA),
  ticker = INST_ABRV,
  value  = ROA.x
)]

# Eliminación de valores faltantes y ordenamiento
df_mk <- df_mk[!is.na(value)]
setorder(df_mk, ticker, date)

############################################################
# 2. Construcción de la matriz de rendimientos
############################################################

# Formato ancho: fechas × activos
R_wide <- dcast(df_mk, date ~ ticker, value.var = "value")

# Eliminar columna de fechas y convertir a matriz
R <- as.matrix(R_wide[, -1, with = FALSE])

# Eliminar activos con exceso de valores faltantes
keep <- colMeans(is.na(R)) < 0.2
R <- R[, keep]

# Eliminar observaciones incompletas
R <- na.omit(R)

############################################################
# 3. Estadísticos descriptivos (anualizados)
############################################################

# Rendimientos esperados y matriz de covarianzas
mu    <- colMeans(R) * 12
Sigma <- cov(R) * 12

############################################################
# 4. Optimización de portafolios
############################################################

n_assets <- length(mu)
ones <- rep(1, n_assets)

# ---- 4.1 Portafolio de mínima varianza (benchmark) ----
Sigma_inv <- solve(Sigma)
w_mvp <- as.vector(
  Sigma_inv %*% ones / as.numeric(t(ones) %*% Sigma_inv %*% ones)
)

# ---- 4.2 Portafolio tangente (máximo ratio de Sharpe) ----
rf <- 0.02  # tasa libre de riesgo anual

excess_mu <- mu - rf
w_tangent <- solve.QP(
  Dmat = 2 * Sigma,
  dvec = rep(0, n_assets),
  Amat = cbind(rep(1, n_assets), excess_mu),
  bvec = c(1, 0),
  meq  = 1
)$solution

# Normalización de los pesos
w_tangent <- w_tangent / sum(w_tangent)

############################################################
# 5. Frontera eficiente (análisis descriptivo)
############################################################

efficient_portfolio <- function(mu, Sigma, target_return) {
  n <- length(mu)
  solve.QP(
    Dmat = 2 * Sigma,
    dvec = rep(0, n),
    Amat = cbind(rep(1, n), mu),
    bvec = c(1, target_return),
    meq  = 2
  )$solution
}

targets <- seq(min(mu) * 0.5, max(mu) * 1.1, length.out = 50)

frontier_sd <- sapply(targets, function(tr) {
  w <- efficient_portfolio(mu, Sigma, tr)
  sqrt(t(w) %*% Sigma %*% w)
})

############################################################
# 6. Simulación de Monte Carlo (portafolio tangente)
############################################################

n_sim <- 5000

sim_returns <- mvrnorm(n = n_sim, mu = mu, Sigma = Sigma)
port_sim <- as.vector(sim_returns %*% w_tangent)

# Medidas de riesgo
VaR_95  <- quantile(port_sim, 0.05)
CVaR_95 <- mean(port_sim[port_sim <= VaR_95])

MC_metrics <- data.frame(
  Mean_Return = mean(port_sim),
  Volatility  = sd(port_sim),
  VaR_95      = VaR_95,
  CVaR_95     = CVaR_95
)

############################################################
# 7. Análisis de escenarios (portafolio tangente)
############################################################

# Choque determinístico sobre los rendimientos esperados
delta <- 0.20

mu_base        <- mu
mu_optimistic  <- mu * (1 + delta)
mu_pessimistic <- mu * (1 - delta)

portfolio_metrics <- function(mu_vec, Sigma, w, rf) {
  ret <- sum(w * mu_vec)
  vol <- sqrt(t(w) %*% Sigma %*% w)
  shr <- (ret - rf) / vol
  c(Return = ret, Volatility = vol, Sharpe = shr)
}

escenarios <- data.frame(
  Scenario = c("Pesimista", "Base", "Optimista"),
  rbind(
    portfolio_metrics(mu_pessimistic, Sigma, w_tangent, rf),
    portfolio_metrics(mu_base,        Sigma, w_tangent, rf),
    portfolio_metrics(mu_optimistic,  Sigma, w_tangent, rf)
  )
)

escenarios[, -1] <- round(escenarios[, -1], 5)

############################################################
# 8. Salidas
############################################################

write.csv(escenarios, "scenario_comparison.csv", row.names = FALSE)
write.csv(MC_metrics, "MC_metrics.csv", row.names = FALSE)
