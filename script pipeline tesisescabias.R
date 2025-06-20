library('fda')        # herramientas para tratar datos funcionales  
library('robustbase') # lmrob 
library(fda.usc)      # dataset + herramientas como fda (derivada, por ejemplo)
library(lattice)      # gráficos
library(gdata)        # uppertriangle
library(splines2)
library(fastmatrix)
library(spacetime)
library(sp)
library(dplyr)
library(glmnet)
library(FRK)
library(matlib)
library(FactoMineR)
library(factoextra)
library(funcharts)
library(pROC)

producto_interno_bases <- function(base1, base2 = NULL, a = 0, b = 1) {
  # Validaciones
  if (!is.matrix(base1)) stop("base1 debe ser una matriz")
  m <- nrow(base1)
  if (is.null(base2)) base2 <- base1
  if (nrow(base2) != m) stop("base1 y base2 deben tener el mismo número de filas (puntos de la grilla)")
  
  # Paso de cuadratura trapezoidal
  h <- (b - a) / (m - 1)
  w <- rep(h, m)
  w[1] <- h / 2
  w[m] <- h / 2
  W <- diag(w)
  
  # Producto interno entre bases
  t(base1) %*% W %*% base2
}

transformar_covariables <- function(X_true, basis = "bspline", K=10, t_grid) {
  limits <- c(min(t_grid),max(t_grid))
  # Lista de bases válidas (según fda.usc y fda)
  valid_basis <- c("fourier", "bspline", "fpca","exponential", "power", "constant")
  
  # Verificar si se pasa una base evaluada (matriz)
  if (is.matrix(basis)) {
    message("Se detectó una matriz evaluada como base. Se usará directamente como Phi_hat.")
    
    A_hat <- t(solve(t(basis) %*% basis) %*% t(basis) %*% X_true)
    Phi_hat <- basis
    Psi_hat <- producto_interno_bases(Phi_hat, a = limits[1], b = limits[2])
    
  } else if (is.character(basis) && length(basis) == 1) {
    # Verificar si la base es válida
    if (!(basis %in% valid_basis)) {
      stop(paste0("La base '", basis, "' no es válida. Opciones posibles: ",
                  paste(valid_basis, collapse = ", "), "."))
    }
    if (!is.matrix(X_true)) {
      message("X_true no  es una matriz esto va a dar error seguro")}
    
    # Crear base funcional
    if (basis == "bspline") {
      estimationBasis = create.bspline.basis(rangeval = limits, nbasis = K ,norder = 6)
    }
    if (basis == "fourier") {
      estimationBasis = create.fourier.basis(rangeval = limits, nbasis =  K)
    }
    if (basis == "fpca") {
      estimationBasis =  fpca.face(X_true, argvals=t_grid, knots=100,pve=0.99)$efunctions
    }
    
    
    if (is.basis(estimationBasis)) 
      {
      fdPar_obj <- fdPar(estimationBasis, 2, 0.01)
      
      # Suavizado de X_true
      smooth <- smooth.basis(argvals = t_grid, y = t(X_true), fdPar_obj)
      A_hat <- t(smooth$fd$coefs)[, 1:K]
      
      # Evaluar base
      Phi_hat <- eval.basis(t_grid, estimationBasis)[, 1:K]
    }
    else{
      
      # Paso 1: pesos de cuadratura (trapezoidal)
      m <- length(t_grid)
      h <- (max(t_grid) - min(t_grid)) / (m - 1)
      w <- rep(h, m)
      w[1] <- h / 2
      w[m] <- h / 2
      W <- diag(w)  # Matriz de pesos
      A_hat <- X_true %*% W %*% estimationBasis  # (n x m) %*% (m x m) %*% (m x K) → (n x K)
      A_hat <- A_hat
      
      Phi_hat <- estimationBasis
    }
    Psi_hat <- producto_interno_bases(Phi_hat, a = limits[1], b = limits[2])
    
  } else {
    stop("El argumento 'basis' debe ser un string válido o una matriz evaluada.")
  }
  
  K <- ncol(A_hat)
  # Crear B_hat a partir de productos cruzados simétricos de A_hat
  B_hat <- array(NA, dim = c(nrow(A_hat), K * (K + 1) / 2))
  for (i in 1:nrow(A_hat)) {
    contador <- 1
    for (j in 1:K) {
      for (k in 1:j) {
        B_hat[i, contador] <- A_hat[i, j] * A_hat[i, k] * 2
        contador <- contador + 1
      }
    }
  }
  
  # Índices para Omega_hat
  mat_idx <- matrix(1:(K * K), nrow = K, ncol = K)
  combos_p <- which(row(mat_idx) >= col(mat_idx), arr.ind = TRUE)
  combos_p <- combos_p[order(combos_p[, 1], combos_p[, 2]), ]
  
  Omega_hat <- matrix(0, nrow = nrow(combos_p), ncol = nrow(combos_p))
  for (m in 1:nrow(combos_p)) {
    for (n in 1:nrow(combos_p)) {
      i1 <- combos_p[m, 1]; i2 <- combos_p[m, 2]
      j1 <- combos_p[n, 1]; j2 <- combos_p[n, 2]
      Omega_hat[m, n] <- Psi_hat[i1, j1] * Psi_hat[i2, j2]
    }
  }
  
  return(list(
    A_hat     = A_hat,
    Psi_hat   = Psi_hat,
    B_hat     = B_hat,
    Omega_hat = Omega_hat,
    Phi_hat   = Phi_hat
  ))
}

realizar_pca_funcional <- function(A_hat, Psi_hat, B_hat, Omega_hat, var_threshold = 0.95) {
  # PCA lineal
  pca_APsi <- prcomp(A_hat %*% Psi_hat, center = TRUE, scale. = TRUE)
  varianza_lin <- pca_APsi$sdev^2
  var_acumulada_lin <- cumsum(varianza_lin) / sum(varianza_lin)
  p_lin <- which(var_acumulada_lin >= var_threshold)[1]
  
  # PCA cuadrático
  pca_BOmega <- prcomp(B_hat %*% Omega_hat, center = TRUE, scale. = TRUE)
  varianza_quad <- pca_BOmega$sdev^2
  var_acumulada_quad <- cumsum(varianza_quad) / sum(varianza_quad)
  p_quad <- which(var_acumulada_quad >= var_threshold)[1]

  
  if (var_threshold == 1 || p_lin > ncol(pca_APsi$rotation) || p_quad > ncol(pca_BOmega$rotation)) {
    p_lin = ncol(pca_APsi$rotation)
    p_quad = ncol(pca_BOmega$rotation)
  }
  
  
  sd_lin <- pca_APsi$scale[1:p_lin]
  centers_pca_APsi <-pca_APsi$center[1:p_lin]
  V_hat <- pca_APsi$rotation[, 1:p_lin]
  Z_lin <- pca_APsi$x[, 1:p_lin]
  
  sd_quad <- pca_BOmega$scale[1:p_quad]
  centers_pca_BOmega <- pca_BOmega$center[1:p_quad]
  W_hat <- pca_BOmega$rotation[, 1:p_quad]
  Z_quad <- pca_BOmega$x[, 1:p_quad]
  
  return(list(
    Z_lin = Z_lin,
    Z_quad = Z_quad,
    V_hat = V_hat,
    W_hat = W_hat,
    p_lin = p_lin,
    p_quad = p_quad,
    centers_lin = centers_pca_APsi,
    centers_quad = centers_pca_BOmega,
    sd_lin = sd_lin,
    sd_quad = sd_quad
  ))
}

fit_models <- function(A_hat, Psi_hat, B_hat, Omega_hat, Y_obs,
                       alpha = 0, crossValglmnet = FALSE,
                       variableSelectionMethod = c("none", "AIC", "LRT"),
                       var_threshold = 0.95) {
  
  variableSelectionMethod <- match.arg(variableSelectionMethod)
  
  # PCA de términos funcionales
  pca_result <- realizar_pca_funcional(A_hat, Psi_hat, B_hat, Omega_hat, var_threshold)
  
  Z_lin <- pca_result$Z_lin
  Z_quad <- pca_result$Z_quad
  Z_all <- cbind(Z_lin, Z_quad)
  
  # Modelos
  if (crossValglmnet) {
    fit_lin <- cv.glmnet(Z_lin, Y_obs, family = "binomial", alpha = alpha)
    fit_quad <- cv.glmnet(Z_all, Y_obs, family = "binomial", alpha = alpha)
  } else {
    Z_lin_df <- as.data.frame(Z_lin)
    names(Z_lin_df) <- paste0("Z_lin", seq_len(ncol(Z_lin_df)))
    
    Z_all_df <- as.data.frame(Z_all)
    names(Z_all_df) <- paste0("Z_all", seq_len(ncol(Z_all_df)))
    
    fit_lin <- glm(Y_obs ~ ., data = Z_lin_df, family = binomial)
    fit_quad <- glm(Y_obs ~ ., data = Z_all_df, family = binomial)
    
  }
  
  # Selección de variables
  if (!crossValglmnet && variableSelectionMethod != "none") {
    if (variableSelectionMethod == "AIC") {
      fit_lin <- step(fit_lin, direction = "forward", test = "AIC")
      fit_quad <- step(fit_quad, direction = "forward", test = "AIC")
    } else if (variableSelectionMethod == "LRT") {
      fit_lin <- step(fit_lin, direction = "forward", test = "LRT")
      fit_quad <- step(fit_quad, direction = "forward", test = "LRT")
    }
  }
  
  return(list(
    fit_lin      = fit_lin,
    fit_quad     = fit_quad,
    Z_lin_scaled = Z_lin,
    Z_all_scaled = Z_all,
    centers_lin  = pca_result$centers_lin,
    centers_all  = c(pca_result$centers_lin, pca_result$centers_quad),
    sd_lin       = pca_result$sd_lin,
    sds_all      = c(pca_result$sd_lin, pca_result$sd_quad),
    V_hat        = pca_result$V_hat,
    W_hat        = pca_result$W_hat,
    p_lin        = pca_result$p_lin,
    p_quad       = pca_result$p_quad
  ))
}



reconstruct_symmetric_matrix <- function(vec, K) {
  Gamma <- matrix(0, nrow = K, ncol = K)
  idx <- 1
  for (i in 1:K) {
    for (j in 1:i) {
      Gamma[i, j] <- vec[idx] / (2-(i==j))
      Gamma[j, i] <- vec[idx] / (2-(i==j)) # simetría
      idx <- idx + 1
    }
  }
  return(Gamma)
}


# 3. Función: reconstruir β_hat y γ_hat ----------------------------------------
reconstruct_functions <- function(Phi_hat, V_hat, W_hat, p_lin, p_quad, fit_lin, fit_quad, sd_lin, sds_all, centers_lin, centers_all) {
  if (p_lin < 1 || p_quad < 1) {
    print("LOS P NO EXISTEN PROBLEMAS DE DIMEN AYUDAAAA")
  }
  M <- nrow(Phi_hat)
  mat_idx  <- matrix(seq_len(ncol(Phi_hat)^2), nrow = ncol(Phi_hat), ncol = ncol(Phi_hat))
  combos_p <- which(row(mat_idx) >= col(mat_idx), arr.ind = TRUE)
  

  Xi_scaled_lin   <- coef(fit_lin)[-1]
  
  if (!is.matrix(Xi_scaled_lin)) {
    Xi_scaled_lin <- matrix(Xi_scaled_lin, ncol = 1)
  }
  
  
  Xi_unscaled_lin <- (Xi_scaled_lin / sd_lin) + centers_lin
  coef_beta_hat_lin_unscaled <- V_hat %*% Xi_unscaled_lin
  beta_hat_lin   <- Phi_hat %*% coef_beta_hat_lin_unscaled
  
  #V_hat_lin_unscaled <- sweep(V_hat /sd_lin, 2, centers_lin, "+")
  #coef_beta_hat_lin_unscaled <- V_hat_lin_unscaled%*%Xi_scaled_lin

  
  coef_scaled_quad    <- coef(fit_quad)[-1]
  Xi_scaled_quad      <- coef_scaled_quad[1:p_lin]
  if (!is.matrix(Xi_scaled_quad)) {
    Xi_scaled_quad <- matrix(Xi_scaled_quad, ncol = 1)
  }
  
  Xi_unscaled_quad <-( Xi_scaled_quad / sds_all[1:p_lin])+centers_all[1:p_lin]
  coef_beta_hat_quad_unscaled <- V_hat %*% Xi_unscaled_quad
  
  beta_hat_quad   <- Phi_hat %*% coef_beta_hat_quad_unscaled

  Kappa_scaled_quad   <- coef_scaled_quad[(p_lin+1):length(coef_scaled_quad)]
  Kappa_unscaled <- (Kappa_scaled_quad / sds_all[(p_lin + 1):length(sds_all)]) +
    centers_all[(p_lin + 1):length(centers_all)]
  
  coef_gamma_hat_quad_unscaled <- W_hat %*% Kappa_unscaled
  
  #W_hat_quad_unscaled <- sweep(W_hat /sds_all[(p_lin+1):length(sds_all)], 2, centers_all[(p_lin+1):length(centers_all)], "+")
  #coef_beta_hat_quad_unscaled <- W_hat_quad_unscaled%*%Kappa_scaled_quad
  gamma_coef <- reconstruct_symmetric_matrix(coef_gamma_hat_quad_unscaled,ncol(Phi_hat))
  
  gamma_hat_quad <- matrix(0, nrow = M, ncol = M)
  for (m in seq_len(nrow(combos_p))) {
    j      <- combos_p[m, 1]
    k      <- combos_p[m, 2]
    coef_t <- gamma_coef[j, k] 
    gamma_hat_quad <- gamma_hat_quad +
      coef_t * (Phi_hat[, j] %o% Phi_hat[, k])
    if (j != k) {
      gamma_hat_quad <- gamma_hat_quad +
        coef_t * (Phi_hat[, k] %o% Phi_hat[, j])
    }
  }
  
  gamma_hat_lin <- matrix(0, nrow = M, ncol = M)
  
  return(list(
    beta_hat_lin   = beta_hat_lin,
    beta_hat_quad  = beta_hat_quad,
    gamma_hat_lin  = gamma_hat_lin,
    gamma_hat_quad = gamma_hat_quad
  ))
}

# 4. Graficar (con poda al 90% interior) y evaluar métricas --------------------
# 4.1. Graficar β(t) en interior: podar 5% en cada extremo
plot_beta_interior <- function(t_grid, beta_true, beta_hat_lin, beta_hat_quad, title_suffix = "") {
  # Definir índices interiores
  interior_idx <- which(t_grid >= 0.05 & t_grid <= 0.95)
  t_int        <- t_grid[interior_idx]
  df_beta <- data.frame(
    t       = t_int,
    True    = beta_true[interior_idx],
    EstLin  = beta_hat_lin[interior_idx],
    EstQuad = beta_hat_quad[interior_idx]
  )
  df_beta_long <- reshape2::melt(df_beta, id.vars = "t",
                                 variable.name = "Modelo",
                                 value.name = "beta")
  
  p <- ggplot(df_beta_long, aes(x = t, y = beta, color = Modelo, linetype = Modelo)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("black", "blue", "red")) +
    scale_linetype_manual(values = c("solid", "dashed", "twodash")) +
    labs(x = "t", y = expression(beta(t)),
         title = paste0("β(t) – 90% interior – ", title_suffix)) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  print(p)
}

# 4.2. Graficar γ(s,t) con contour en interior (podar bordes)
plot_gamma_contour_interior <- function(t_grid, gamma_true, gamma_hat_lin, gamma_hat_quad, title_suffix = "") {
  # Identificar índices interiores
  interior_idx <- which(t_grid >= 0.05 & t_grid <= 0.95)
  
  # Recortar las matrices gamma a interior_idx × interior_idx
  gt_int_ll <- gamma_true[interior_idx, interior_idx, drop = FALSE]
  gl_int_ll <- gamma_hat_lin[interior_idx, interior_idx, drop = FALSE]
  gq_int_ll <- gamma_hat_quad[interior_idx, interior_idx, drop = FALSE]
  t_int     <- t_grid[interior_idx]
  
  prep_gamma <- function(mat, tipo) {
    df <- reshape2::melt(mat)
    colnames(df) <- c("s_idx", "t_idx", "gamma")
    df$Tipo <- tipo
    df$s    <- t_int[df$s_idx]
    df$t    <- t_int[df$t_idx]
    return(df)
  }
  
  df_gamma_true <- prep_gamma(gt_int_ll, paste0("True – ", title_suffix))
  df_gamma_lin  <- prep_gamma(gl_int_ll, paste0("Est. Lin. – ", title_suffix))
  df_gamma_quad <- prep_gamma(gq_int_ll, paste0("Est. Quad. – ", title_suffix))
  df_gamma_all  <- rbind(df_gamma_true, df_gamma_lin, df_gamma_quad)
  
  p <- ggplot(df_gamma_all, aes(x = t, y = s, z = gamma)) +
    geom_contour_filled(aes(fill = ..level..), bins = 20) +
    facet_wrap(~ Tipo, nrow = 1) +
    scale_fill_viridis_d(option = "plasma", name = expression(gamma(s, t))) +
    labs(x = "t", y = "s", title = paste0("Contour γ(s,t) – 90% interior – ", title_suffix)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          strip.text = element_text(size = 10),
          axis.text  = element_text(size = 6),
          axis.title = element_text(size = 8))
  
  print(p)
}

plot_gamma_3d_interior <- function(t_grid, gamma_true, gamma_hat_quad, title_suffix = "") {
  # Encontrar índices interiores
  interior_idx <- which(t_grid >= 0.05 & t_grid <= 0.95)
  t_int        <- t_grid[interior_idx]
  
  # Recortar matrices
  gt_int <- gamma_true[interior_idx, interior_idx, drop = FALSE]
  gq_int <- gamma_hat_quad[interior_idx, interior_idx, drop = FALSE]
  
  theta   <- 30
  phi_ang <- 30
  col_t   <- "lightblue"
  col_q   <- "lightcoral"
  
  # Calcular zmin y zmax combinados para las dos submatrices
  zmin <- min(c(gt_int, gq_int), na.rm = TRUE)
  zmax <- max(c(gt_int, gq_int), na.rm = TRUE)
  if (abs(zmax - zmin) < 1e-8) {
    zmin <- zmin - 1e-6
    zmax <- zmax + 1e-6
  }
  
  op <- par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
  
  persp(x      = t_int,
        y      = t_int,
        z      = gt_int,
        theta  = theta,
        phi    = phi_ang,
        expand = 1,
        col    = col_t,
        zlim   = c(zmin, zmax),
        ticktype = "detailed",
        xlab   = "t",
        ylab   = "s",
        zlab   = expression(gamma(s, t)),
        main   = paste0("γ_true – 90% interior – ", title_suffix),
        cex.main = 1.1, cex.lab = 0.9, cex.axis = 0.7)
  
  persp(x      = t_int,
        y      = t_int,
        z      = gq_int,
        theta  = theta,
        phi    = phi_ang,
        expand = 1,
        col    = col_q,
        zlim   = c(zmin, zmax),
        ticktype = "detailed",
        xlab   = "t",
        ylab   = "s",
        zlab   = expression(hat(gamma)(s, t)),
        main   = paste0("γ_hat_quad – 90% interior – ", title_suffix),
        cex.main = 1.1, cex.lab = 0.9, cex.axis = 0.7)
  
  par(op)
}

# 4.4. Matriz de confusión y CCR
calculate_confusion <- function(y_true, y_pred) {
  TP  <- sum(y_true == 1 & y_pred == 1)
  TN  <- sum(y_true == 0 & y_pred == 0)
  FP  <- sum(y_true == 0 & y_pred == 1)
  FN  <- sum(y_true == 1 & y_pred == 0)
  CCR <- (TP + TN) / length(y_true)
  return(list(TN = TN, FP = FP, FN = FN, TP = TP, CCR = CCR))
}

# 4.5. Métricas: matriz de confusión, CCR, ROC y AUC
compute_metrics <- function(fit_lin, fit_quad, Y_val, Z_lin_val, Z_all_val,
                            crossValglmnet, plot_rocs = TRUE) {
  if (crossValglmnet) {
    prob_lin   <- predict(fit_lin, newx = Z_lin_val, type = "response", s = "lambda.min")
    prob_quad  <- predict(fit_quad, newx = Z_all_val, type = "response", s = "lambda.min")
  } else {
    Z_lin_val_df <- as.data.frame(Z_lin_val)
    names(Z_lin_val_df) <- paste0("Z_lin", seq_len(ncol(Z_lin_val_df)))
    
    Z_all_val_df <- as.data.frame(Z_all_val)
    names(Z_all_val_df) <- paste0("Z_all", seq_len(ncol(Z_all_val_df)))
    
    prob_lin <- predict(fit_lin, newdata = Z_lin_val_df, type = "response")
    prob_quad <- predict(fit_quad, newdata = Z_all_val_df, type = "response")
  }
  
  pred_lin  <- ifelse(prob_lin > 0.5, 1, 0)
  pred_quad <- ifelse(prob_quad > 0.5, 1, 0)
  
  res_lin   <- calculate_confusion(Y_val, pred_lin)
  res_quad  <- calculate_confusion(Y_val, pred_quad)
  
  roc_lin   <- roc(response = Y_val, predictor = as.vector(prob_lin))
  roc_quad  <- roc(response = Y_val, predictor = as.vector(prob_quad))
  
  auc_lin   <- auc(roc_lin)
  auc_quad  <- auc(roc_quad)
  
  if (plot_rocs) {
    plot(roc_lin,  col = "blue", lwd = 2,
         main = "Curva ROC Validación",
         legacy.axes = TRUE)
    lines(roc_quad, col = "red", lwd = 2)
    abline(a = 0, b = 1, lty = 3, col = "gray40")
    legend("bottomright",
           legend = c(
             sprintf("Lineal  (AUC = %.3f)", auc_lin),
             sprintf("Cuadrático (AUC = %.3f)", auc_quad)
           ),
           col  = c("blue", "red"),
           lwd  = 2,
           bty  = "n")
  }
  
  return(list(
    confusion_lin   = res_lin,
    confusion_quad  = res_quad,
    auc_lin         = auc_lin,
    auc_quad        = auc_quad,
    roc_lin         = roc_lin,
    roc_quad        = roc_quad
  ))
}



simular_datos_funcionales_binarios <- function(n = 200, M = 101, K_true = 10,
                                               var_threshold = 0.95,
                                               lambda_fun = function(k) k^-2,
                                               beta_coef_custom = NULL,
                                               gamma_coef_custom = NULL) {
  
  # 1) Parámetros y base
  t_grid <- seq(0, 1, length.out = M)
  Delta_t <- t_grid[2] - t_grid[1]
  
  #Bs_basis <- create.bspline.basis(c(0, 1), K_true)
  #Phi_true <- eval.basis(t_grid, Bs_basis)
  
  Phi_true <- matrix(0, nrow = M, ncol = K_true)
  for (k in 1:K_true) {
    Phi_true[, k] <- sqrt(2) * sin((k - 0.5) * pi * t_grid)
  }
  
  
  
  lambda_k <- lambda_fun(1:K_true)
  means <- runif(K_true, min = -0.001, max = 0.001)
  stdvs <- lambda_k
  
  # 2) Simular coeficientes A
  A_true <- matrix(rnorm(n * K_true, mean = rep(means, each = n), sd = rep(stdvs, each = n)),
                   nrow = n, ncol = K_true)
  
  # 3) Producto interno Psi (matriz de Gram de la base)
  "  Psi_true <- inprod(Bs_basis,
                     Bs_basis)"
  Psi_true <- producto_interno_bases(Phi_true)
  # 4) Escalar A y reconstruir X
  A_means <- colMeans(A_true)
  A_sds <- apply(A_true, 2, sd)
  A_true_scaled <- scale(A_true, center = A_means, scale = A_sds)
  X_true <- A_true_scaled %*% t(Phi_true)
  
  # 5) Construir B_true_scaled para términos cuadráticos
  B_true_scaled <- matrix(NA, nrow = n, ncol = K_true * (K_true + 1) / 2)
  contador <- 1
  for (j in 1:K_true) {
    for (k in 1:j) {
      B_true_scaled[, contador] <- (2 - (j == k)) * A_true_scaled[, j] * A_true_scaled[, k]
      contador <- contador + 1
    }
  }
  
  # 6) Matriz Omega_true (producto tensorial cuadrático)
  mat_idx <- matrix(1:(K_true * K_true), nrow = K_true, ncol = K_true)
  combos <- which(row(mat_idx) >= col(mat_idx), arr.ind = TRUE)
  combos <- combos[order(combos[, 1], combos[, 2]), ]
  
  Omega_true <- matrix(0, nrow = ncol(B_true_scaled), ncol = ncol(B_true_scaled))
  for (m in 1:nrow(Omega_true)) {
    for (ñ in 1:ncol(Omega_true)) {
      i1 <- combos[m, 1]
      i2 <- combos[m, 2]
      j1 <- combos[ñ, 1]
      j2 <- combos[ñ, 2]
      Omega_true[m, ñ] <- Psi_true[i1, j1] * Psi_true[i2, j2]
    }
  }
  
  # 7) Predictor lineal + cuadrático
  gamma_lower_by_row <- gamma_coef_A[lower.tri(gamma_coef_A, diag = TRUE)]
  Lin_part <- A_true_scaled %*% Psi_true %*% beta_coef_A
  Quad_part <- B_true_scaled %*% Omega_true %*% gamma_lower_by_row
  eta_true <- Lin_part + Quad_part
  
  # 8) Simulación de respuestas binarias
  prob_true <- exp(eta_true) / (1 + exp(eta_true))
  Y_obs <- rbinom(n, size = 1, prob = prob_true)
  
  # 9) beta(t)
  beta_true <- Phi_true %*% beta_coef_custom
  
  # 10) gamma(s,t)
  gamma_true <- matrix(0, nrow = M, ncol = M)
  for (m in seq_len(nrow(combos))) {
    j <- combos[m, 1]
    k <- combos[m, 2]
    coef_t <- gamma_coef_custom[j, k] * (2 - ifelse(j == k, 1, 0))
    gamma_true <- gamma_true +
      coef_t * (Phi_true[, j] %o% Phi_true[, k])
    if (j != k) {
      gamma_true <- gamma_true +
        coef_t * (Phi_true[, k] %o% Phi_true[, j])
    }
  }
  
  return(list(
    t_grid = t_grid,
    Delta_t = Delta_t,
    K_true = K_true,
    Phi_true = Phi_true,
    X_true = X_true,
    Y_obs = Y_obs,
    beta_true = beta_true,
    gamma_true = gamma_true
  ))
}

# =====================================
# Script: Grid Search con Validación Cruzada para Modelos Funcionales
# =====================================

# Asume que las funciones necesarias están definidas en el entorno:
# - transformar_covariables
# - fit_models
# - reconstruct_functions
# - compute_metrics
# - plot_beta_interior
# - plot_gamma_contour_interior
# - plot_gamma_3d_interior

ajustar_modelo_funcional_escabias <- function(t_grid, X_true, Y_obs, beta_true, gamma_true, Phi_true,
                                     var_threshold = 0.99,
                                     K,
                                     basis4smoothing = "bspline",
                                     alpha = 0,
                                     limits = c(0,1),
                                     crossValglmnet = FALSE,
                                     variableSelection = "LRT",
                                     escenario = "modelo funcional cuadrático",
                                     plot_results = TRUE) {
  covariables_estimadas <- transformar_covariables(X_true = X_true, basis = basis4smoothing, K = K, t_grid = t_grid)
  
  fits <- fit_models(
    covariables_estimadas$A_hat,
    covariables_estimadas$Psi_hat,
    covariables_estimadas$B_hat,
    covariables_estimadas$Omega_hat,
    Y_obs,
    alpha = alpha,
    crossValglmnet = crossValglmnet,
    variableSelection = variableSelection,
    var_threshold = var_threshold
  )
  
  if ((fits$p_lin < 1)|| (fits$p_quad < 1)) {
    print((fits$p_lin < 1)|| (fits$p_quad < 1))
  }
  
  rec <- reconstruct_functions(
    Phi_hat       = covariables_estimadas$Phi_hat,
    V_hat         = fits$V_hat,
    W_hat         = fits$W_hat,
    p_lin         = fits$p_lin,
    p_quad        = fits$p_quad,
    fit_lin       = fits$fit_lin,
    fit_quad      = fits$fit_quad,
    sd_lin        = fits$sd_lin,
    sds_all       = fits$sds_all,
    centers_lin   = fits$centers_lin,
    centers_all   = fits$centers_all
  )
  
  if (plot_results) {
    plot_beta_interior(t_grid, beta_true, rec$beta_hat_lin, rec$beta_hat_quad,
                       title_suffix = paste0(" (p=", fits$p_lin, ",", fits$p_quad, ")"))
    plot_gamma_contour_interior(t_grid, gamma_true, rec$gamma_hat_lin, rec$gamma_hat_quad,
                                title_suffix = paste0(" (p=", fits$p_lin, ",", fits$p_quad, ")"))
    plot_gamma_3d_interior(t_grid, gamma_true, rec$gamma_hat_quad,
                           title_suffix = paste0(" (p=", fits$p_lin, ",", fits$p_quad, ")"))
    par(mfrow = c(1, 1))
  }
  
  return(list(
    Z_lin_scaled = fits$Z_lin_scaled,
    Z_all_scaled = fits$Z_all_scaled,
    fit_lin = fits$fit_lin,
    fit_quad = fits$fit_quad,
    crossValglmnet = crossValglmnet,
    p_lin         = fits$p_lin,
    p_quad        = fits$p_quad
  ))
}



evaluar_modelo_funcional <- function(modelo, X, Y, basis, K, t_grid, p_lin, p_quad, plot = FALSE) {
  covariables_estimadas <- transformar_covariables(X_true =  X, basis = basis, K = K, t_grid = t_grid)
  
  pca_result <- realizar_pca_funcional(A_hat = covariables_estimadas$A_hat,
                                       Psi_hat = covariables_estimadas$Psi_hat,
                                       B_hat = covariables_estimadas$B_hat,
                                       Omega_hat = covariables_estimadas$Omega_hat,
                                       var_threshold = 1)
  
  p_lin_valid <- min(p_lin,ncol(pca_result$Z_lin))
  p_quad_valid <- min(p_quad,ncol(pca_result$Z_quad))
  
  if (p_lin > ncol(pca_result$Z_lin) || p_quad > ncol(pca_result$Z_quad)) {
    print("Esto da error por que la cantidad de parametros en validación es menor a la cantidad de PC's necesarias")
  }
  Z_lin_val <- pca_result$Z_lin[, 1:p_lin_valid, drop = FALSE]
  Z_quad_val <- pca_result$Z_quad[, 1:p_quad_valid, drop = FALSE]
  Z_all_val <- cbind(Z_lin_val, Z_quad_val)
  
  
  res <- compute_metrics(
    fit_lin        = modelo$fit_lin,
    fit_quad       = modelo$fit_quad,
    Y_val          = Y,
    Z_lin_val   = Z_lin_val,
    Z_all_val   = Z_all_val,
    crossValglmnet = modelo$crossValglmnet,
    plot_rocs      = plot
  )

  return(res)
}

#En este grid_search_Cv hay que tener en cuenta que la cantidad de datos (filas) en 
#los conjuntos de validación tienen que ser mayores que la cantidad de p_lin y p_quad
#parametros explicativos (o, en el mejor caso, más datos que columnas en general).
#Sino, al hacer PCA en los conjuntos de validación vamos a tener problemas para predecir
# ya que el modelo armado va a esperar parametros que no están
grid_search_cv_Escabias <- function(X, Y, t_grid, beta_true, gamma_true, Phi_true,
                           var_threshold_vals, K_vals, basis_vals,
                           variableSelection_vals, crossVal_vals, alpha_vals,
                           nfolds = 3, seed = 123, limits=c(0,1)) {
  
  
  set.seed(seed)
  n <- nrow(X)
  
  # Dividir en conjunto de entrenamiento y conjunto de test final (20%)
  idx_total <- 1:n
  idx_test_final <- sample(idx_total, size = round(0.2 * n))
  idx_train_total <- setdiff(idx_total, idx_test_final)
  
  X_train_total <- X[idx_train_total, ]
  Y_train_total <- Y[idx_train_total]
  X_test_final <- X[idx_test_final, ]
  Y_test_final <- Y[idx_test_final]
  
  # Crear folds de cross-validation sobre el conjunto de entrenamiento
  folds <- sample(rep(1:nfolds, length.out = length(idx_train_total)))
  
  resultados_grid <- list()
  contador <- 1
  
  for (var_threshold in var_threshold_vals) {
    for (K in K_vals) {
      for (basis in basis_vals) {
        for (cv in crossVal_vals) {
          for (vs in variableSelection_vals) {
            for (a in alpha_vals) {
              if ((cv && vs != "none") || (!cv && a != 0)) next
              
              cat(sprintf("Parámetros #%d: vt=%.2f, K=%d, basis=%s, CV=%s, VS=%s, alpha=%.2f\n",
                          contador, var_threshold, K, basis, cv, vs, a))
              
              metricas_fold <- list()
              
              for (k in 1:nfolds) {
                idx_valid <- which(folds == k)
                idx_train <- setdiff(1:length(idx_train_total), idx_valid)
                
                X_train <- X_train_total[idx_train, ]
                Y_train <- Y_train_total[idx_train]
                X_valid <- X_train_total[idx_valid, ]
                Y_valid <- Y_train_total[idx_valid]
                
                
                res <- tryCatch({
                  ajustar_modelo_funcional_escabias(
                    t_grid = t_grid,
                    X_true = X_train,
                    Y_obs = Y_train,
                    beta_true = beta_true,
                    gamma_true = gamma_true,
                    Phi_true = Phi_true,
                    var_threshold = var_threshold,
                    K = K,
                    basis4smoothing = basis,
                    crossValglmnet = cv,
                    variableSelection = vs,
                    alpha = a,
                    plot_results = FALSE
                  )
                }, error = function(e) {
                  message(sprintf("Error en fold %d: %s", k, e$message))
                  return(NULL)
                })
                
                # ⛔ Si hay error en el modelo, pasar al siguiente fold
                if (is.null(res)) {
                  next
                }
                
                metrica <- evaluar_modelo_funcional(modelo = res,
                                                    X = X_valid,
                                                    Y = Y_valid,
                                                    basis = basis, 
                                                    K = K, 
                                                    t_grid = t_grid,
                                                    p_lin = res$p_lin,
                                                    p_quad = res$p_quad,
                                                    plot = FALSE)
                
                metricas_fold[[k]] <- metrica  
                
              }
              
              
              resultados_grid[[contador]] <- list(
                var_threshold = var_threshold,
                K = K,
                basis = basis,
                crossValglmnet = cv,
                variableSelection = vs,
                alpha = a,
                metricas = metricas_fold
              )
              contador <- contador + 1
            }
          }
        }
      }
    }
  }
  
  
  # Crear un vector para guardar el promedio de AUROC de cada combinación
  auroc_promedios <- numeric(length(resultados_grid))
  
  for (i in seq_along(resultados_grid)) {
    metricas_folds <- resultados_grid[[i]]$metricas
    
    # Filtrar folds que no sean NULL (en caso de errores)
    metricas_folds <- Filter(Negate(is.null), metricas_folds)
    
    # Extraer los AUROC (podés cambiar auc_quad por auc_lin si querés comparar el modelo lineal)
    aurocs <- sapply(metricas_folds, function(m) m$auc_quad)
    
    # Calcular el promedio
    auroc_promedios[i] <- mean(aurocs)
  }
  
  # Encontrar el índice del mejor modelo
  mejor_indice <- which.max(auroc_promedios)
  
  # Extraer la mejor combinación de parámetros
  mejor_modelo <- resultados_grid[[mejor_indice]]
  
  # Mostrar el mejor AUROC promedio
  cat("Mejor AUROC promedio:", auroc_promedios[mejor_indice], "\n")
  cat("Parámetros del mejor modelo:\n")
  print(mejor_modelo[c("var_threshold", "K", "basis", "crossValglmnet", "variableSelection", "alpha")])
  
  cat("\nGraficando mejor modelo...\n")

  res <- ajustar_modelo_funcional_escabias(
    t_grid = t_grid,
    X_true = X_train_total,
    Y_obs = Y_train_total,
    beta_true = beta_true,
    gamma_true = gamma_true,
    Phi_true = Phi_true,
    var_threshold = mejor_modelo$var_threshold,
    K = mejor_modelo$K,
    basis4smoothing = mejor_modelo$basis,
    crossValglmnet = mejor_modelo$crossValglmnet,
    variableSelection = mejor_modelo$variableSelection,
    alpha = mejor_modelo$alpha,
    plot_results = TRUE
  )
  
  evaluar_modelo_funcional(modelo = res,
                           X = X_test_final,
                           Y = Y_test_final,
                           basis = mejor_modelo$basis, 
                           K = mejor_modelo$K, 
                           t_grid = t_grid,
                           p_lin = res$p_lin,
                           p_quad = res$p_quad,
                           plot = TRUE)
  
  
  
  return(resultados_grid)
}


# ===================
# simulo los datos
#====================
sim <- simular_datos_funcionales_binarios(gamma_coef_custom = gamma_coef_A,
                                          beta_coef_custom = beta_coef_A,
                                          K_true = K_true,
                                          n = 3000, M=201)

var_threshold_vals = c(0.75,0.8, 0.85, 0.9, 0.95, 0.99)
variableSelection_vals = c("none", "LRT", "AIC")
crossVal_vals = c(FALSE, TRUE)
alpha_vals = c(0,0.5,1)


# ===================
# Definición de parámetros del grid
# ===================
var_threshold_vals     <- c(0.85,0.9,0.95,0.99)
K_vals                 <-  seq(25,37,3) #25:35#c(1:10) #4:30
basis_vals             <- "bspline" #c("fourier", "bspline")
variableSelection_vals <- c("none", "LRT", "AIC")
crossVal_vals          <- FALSE#c(FALSE, TRUE)
alpha_vals             <- 0#c(0,0.25,0.5,0.75, 1)
#debug(reconstruct_functions)


grid_caso_1 <- grid_search_cv_Escabias(X = sim$X_true,
                                         Y = sim$Y_obs,
                                         beta_true =  sim$beta_true,
                                         gamma_true = sim$gamma_true,
                                         Phi_true = sim$Phi_true,
                                         var_threshold_vals =  var_threshold_vals,
                                         K_vals =K_vals,
                                         basis_vals =  basis_vals,
                                         crossVal_vals= crossVal_vals,
                                         alpha_vals =  alpha_vals,
                                         variableSelection_vals = variableSelection_vals,
                                         t_grid = sim$t_grid)



K_true <- 10
beta_coef_A <- numeric(10)
beta_coef_A[1] <- 2
beta_coef_A[2] <- -2.5
gamma_coef_A <- matrix(0, nrow = 10, ncol = 10)
gamma_coef_A[1, 1] <- 6.5
gamma_coef_A[1, 3] <- -5.0
gamma_coef_A[2, 3] <- -4.3
gamma_coef_A[2, 2] <- -7.0
gamma_coef_A[3, 5] <- 10.0

sim <- simular_datos_funcionales_binarios(gamma_coef_custom = gamma_coef_A,
                                          beta_coef_custom = beta_coef_A,
                                          K_true = K_true,
                                          n = 3000, M=201)


grid_caso_2 <- grid_search_cv_Escabias(X = sim$X_true,
                                       Y = sim$Y_obs,
                                       beta_true =  sim$beta_true,
                                       gamma_true = sim$gamma_true,
                                       Phi_true = sim$sPhi_true,
                                       var_threshold_vals =  var_threshold_vals,
                                       K_vals =K_vals,
                                       basis_vals =  basis_vals,
                                       crossVal_vals= crossVal_vals,
                                       alpha_vals =  alpha_vals,
                                       variableSelection_vals = variableSelection_vals,
                                       t_grid = sim$t_grid)

