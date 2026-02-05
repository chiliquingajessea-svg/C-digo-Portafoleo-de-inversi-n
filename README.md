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

names(df1)
names(df2)

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


# 3. Limpieza básica

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


# 4. Estadística descriptiva

# Media y covarianza (mensual → anual)
mu <- colMeans(R) * 12
Sigma <- cov(R) * 12


#  5. Optimización: frontera eficiente
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

# 6. Portafolio tangente (max Sharpe)
rf <- 0.02 # tasa libre de riesgo anual;
excess_mu <- mu - rf
w_tangent <- solve.QP(2*Sigma, rep(0,length(mu)), cbind(rep(1,length(mu)), excess_mu),
                      c(1,0), meq = 1)$solution

mu_p_M <- sum(w_tangent * mu)
mu_p_M
sigma_p_M <- sqrt(t(w_tangent) %*% Sigma %*% w_tangent)
sigma_p_M

# Normalizar si fuera necesario
w_tangent <- w_tangent / sum(w_tangent)

# 7. Simulaciones Monte Carlo (multivar normal como ejemplo) ---
n_sim <- 5000
n_steps <- 1  # un periodo para distribución de retornos anuales
sim_returns <- mvrnorm(n = n_sim, mu = mu, Sigma = Sigma)
# Rendimiento del portafolio tangente en cada simulación
port_ret_sim <- as.vector(sim_returns %*% w_tangent)
# Métricas
VaR_95 <- quantile(port_ret_sim, probs = 0.05)
CVaR_95 <- mean(port_ret_sim[port_ret_sim <= VaR_95])

# 8. Visualizaciones ---
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


# ================================
# Métricas Monte Carlo
# ================================

mu_p_MC     <- mean(port_ret_sim)
sigma_p_MC <- sd(port_ret_sim)

VaR_95  <- quantile(port_ret_sim, 0.05)
CVaR_95 <- mean(port_ret_sim[port_ret_sim <= VaR_95])

# Resumen
MC_metrics <- data.frame(
  Media_MC = mu_p_MC,
  Volatilidad_MC = sigma_p_MC,
  VaR_95 = VaR_95,
  CVaR_95 = CVaR_95
)

print(MC_metrics)
write.csv(MC_metrics, "MC_metrics.csv", row.names = FALSE)



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
# Gráficos
# =============================

# Histograma de rendimientos simulados del portafolio
library(ggplot2)

df_hist <- data.frame(Rendimiento = port_ret_sim)

g_hist <- ggplot(df_hist, aes(x = Rendimiento)) +
  geom_histogram(
    bins = 50,
    fill = "gray70",
    color = "black"
  ) +
  labs(
    title = "Distribución de los rendimientos simulados del portafolio",
    x = "Rendimiento anual simulado",
    y = "Frecuencia"
  ) +
  theme_minimal()

print(g_hist)


# Curva de densidad + VaR y CVaR
ggplot(data.frame(R = port_ret_sim), aes(R)) +
  geom_density(fill = "gray70", alpha = 0.5) +
  geom_vline(xintercept = VaR_95, linetype = "dashed") +
  geom_vline(xintercept = CVaR_95, linetype = "dotted") +
  labs(x = "Rendimiento", y = "Densidad")

# Trayectoria media y bandas percentiles
g_traj <- ggplot(traj_df, aes(x = Time)) +
  geom_ribbon(
    aes(ymin = P5, ymax = P95),
    fill = "gray80",
    alpha = 0.5
  ) +
  geom_line(
    aes(y = Mean),
    linewidth = 1,
    color = "black"
  ) +
  labs(
    title = "Trayectoria media y bandas percentiles (5%–95%)",
    x = "Tiempo",
    y = "Valor del portafolio"
  ) +
  theme_minimal()

print(g_traj)


#  Rendimiento vs Riesgo

ggplot(comp_df, aes(x = Volatilidad, y = Rendimiento, label = Metodo)) +
  geom_point(size = 2.1) +
  geom_text(vjust = -1, size = 2.1) +
  labs(
    title = "Comparación de desempeño entre metodologías",
    x = "Volatilidad del portafolio",
    y = "Rendimiento esperado"
  ) +
  theme_minimal()






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


ROA_matrix <- as.matrix(tabla_rendimientos)

# R ya es tu matriz de rendimientos (filas = fechas, columnas = bancos)
ROA_matrix <- R

# Media mensual/esperada (anualizada)
colMeans(ROA_matrix) * 12

# Desviación estándar mensual/esperada (anualizada)
apply(ROA_matrix, 2, sd) * sqrt(12)

# Matriz de covarianza anualizada
cov(ROA_matrix) * 12


# Instala los paquetes si no los tienes
install.packages("gridExtra")
install.packages("grid")

########################################

# ================================
# Trayectorias Monte Carlo en R
# ================================

library(MASS)     # mvrnorm
library(ggplot2)
library(reshape2)

set.seed(123)

# Parámetros
T <- 250                 # horizonte temporal (ej. días)
n_sim <- 300             # número de trayectorias
V0 <- 1                  # valor inicial del portafolio


# mu_vec  -> vector de rendimientos esperados (21x1)
mu_vec <- colMeans(R, na.rm = TRUE)
# Sigma   -> matriz de covarianza (21x21)
Sigma <- cov(R, use = "complete.obs")
# w_opt   -> pesos del portafolio óptimo (21x1)
# Número de activos
n <- ncol(R)

# Vector de unos
ones <- rep(1, n)

# Inversa de la covarianza
Sigma_inv <- solve(Sigma)

# Pesos del portafolio de mínima varianza
w_opt <- Sigma_inv %*% ones / as.numeric(t(ones) %*% Sigma_inv %*% ones)

# Convertir a vector
w_opt <- as.vector(w_opt)



# 1. Simulación de rendimientos multivariados
R_sim <- mvrnorm(n = T * n_sim, mu = mu_vec, Sigma = Sigma)

# 2. Rendimiento del portafolio
R_port <- R_sim %*% w_opt

# 3. Reorganizar en matriz T x n_sim
R_port_mat <- matrix(R_port, nrow = T, ncol = n_sim)

# 4. Valor del portafolio (trayectorias)
V_port <- apply(R_port_mat, 2, function(r) {
  V0 * cumprod(1 + r)
})

# 5. Pasar a data frame para ggplot
V_df <- as.data.frame(V_port)
V_df$Time <- 1:T

V_long <- melt(V_df, id.vars = "Time")

# 6. Gráfico de trayectorias
ggplot(V_long, aes(x = Time, y = value, group = variable)) +
  geom_line(alpha = 0.15) +
  labs(title = "Trayectorias simuladas del valor del portafolio",
       x = "Tiempo",
       y = "Valor del portafolio") +
  theme_minimal()

#######################   continuando....

# Trayectoria media y percentiles
mean_path <- rowMeans(V_port)
p5 <- apply(V_port, 1, quantile, 0.05)
p95 <- apply(V_port, 1, quantile, 0.95)

traj_df <- data.frame(
  Time = 1:T,
  Mean = mean_path,
  P5 = p5,
  P95 = p95
)

ggplot(traj_df, aes(x = Time)) +
  geom_line(aes(y = Mean), linewidth = 1) +
  geom_ribbon(aes(ymin = P5, ymax = P95), alpha = 0.3) +
  labs(title = "Trayectoria media y bandas percentiles (5%–95%)",
       x = "Tiempo",
       y = "Valor del portafolio") +
  theme_minimal()


######################################### Representación gráfica de carteras
################################# Portafolios simulados + Frontera eficiente + Portafolio óptimo

library(ggplot2)

set.seed(123)

n_assets <- length(mu_vec)
n_port   <- 10000

W <- matrix(runif(n_assets * n_port),
            nrow = n_port,
            ncol = n_assets)

# Normalizar pesos
W <- W / rowSums(W)


# Riesgo y rendimiento del portafolio
risk_p  <- sqrt(diag(W %*% Sigma %*% t(W)))
ret_p   <- W %*% mu_vec

df_port <- data.frame(
  Riesgo = risk_p,
  Retorno = ret_p
)

df_opt <- data.frame(
  Riesgo = sqrt(t(w_opt) %*% Sigma %*% w_opt),
  Retorno = sum(w_opt * mu_vec)
)

ggplot(df_port, aes(Riesgo, Retorno)) +
  geom_point(alpha = 0.4, color = "gray50") +
  geom_point(data = df_opt, aes(Riesgo, Retorno),
             color = "red", size = 3) +
  labs(
    title = "Portafolios simulados y portafolio óptimo",
    x = "Riesgo (Desviación estándar)",
    y = "Rendimiento esperado"
  ) +
  theme_minimal()

########################################Gráfica de trayectorias del portafolio

# Serie de rendimientos del portafolio óptimo
returns_portfolio <- R %*% w_opt

capital <- cumprod(1 + returns_portfolio)

df_traj <- data.frame(
  Tiempo = 1:length(capital),
  Capital = capital
)

ggplot(df_traj, aes(Tiempo, Capital)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(
    title = "Trayectoria temporal del portafolio óptimo",
    x = "Tiempo",
    y = "Capital acumulado"
  ) +
  theme_minimal()

######################################### Gráfica de composición del portafolio

df_w <- data.frame(
  Activo = colnames(R),
  Peso   = as.numeric(w_opt)
)

df_w2 <- df_w
df_w2$Activo <- ifelse(df_w2$Activo == "BP BANCO  DESARROLLO DE LOS PUEBLOS  S.A., CODESARROLLO", "BP DESARROLLO", df_w2$Activo)

ggplot(df_w2, aes(x = reorder(Activo, Peso), y = Peso)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  labs(
    title = "Composición del portafolio óptimo",
    x = "Activo",
    y = "Peso"
  ) +
  theme_minimal()

ggplot(df_w2, aes(x = reorder(Activo, Peso), y = Peso, fill = Peso)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_viridis_c(option = "D") +
  labs(
    title = "Composición del portafolio óptimo",
    x = "Activo",
    y = "Peso",
    fill = "Peso"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
