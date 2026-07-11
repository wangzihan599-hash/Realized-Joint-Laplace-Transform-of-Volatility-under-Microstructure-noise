get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

parse_uv_grid <- function(env_name, default_grid) {
  raw_value <- Sys.getenv(env_name, unset = "")
  if (raw_value == "") {
    return(default_grid)
  }

  pairs <- strsplit(raw_value, ",")[[1]]
  uv_mat <- do.call(
    rbind,
    lapply(pairs, function(pair) {
      parts <- strsplit(trimws(pair), ":")[[1]]
      if (length(parts) != 2L) {
        stop("Each UV pair must have format u:v")
      }
      c(u = as.numeric(parts[1]), v = as.numeric(parts[2]))
    })
  )

  data.frame(
    u = as.numeric(uv_mat[, "u"]),
    v = as.numeric(uv_mat[, "v"])
  )
}

discover_energy_onesecond_files <- function(root_dir) {
  files <- list.files(
    root_dir,
    pattern = "_onesecond\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  rel_paths <- substring(files, nchar(normalizePath(root_dir, winslash = "/")) + 2L)
  parts <- strsplit(rel_paths, .Platform$file.sep, fixed = TRUE)

  out <- lapply(seq_along(files), function(i) {
    folder <- parts[[i]][1L]
    file_name <- basename(files[i])
    ticker <- sub("data.*$", "", folder)
    date_value <- strsplit(file_name, "_", fixed = TRUE)[[1]][2L]
    period <- if (date_value < "2010-01-01") "2007_2008" else "2016_2017"

    data.frame(
      path = normalizePath(files[i], winslash = "/"),
      folder = folder,
      ticker = ticker,
      date = date_value,
      period = period,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

load_onesecond_logprice <- function(path) {
  dat <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("date", "time", "price") %in% names(dat))) {
    stop("Unexpected onesecond file format: ", path)
  }
  log(as.numeric(dat$price))
}

prepare_single_preaverage <- function(y,
                                      theta,
                                      times = NULL,
                                      noise_sd = NA_real_,
                                      clip_multiplier = Inf) {
  g_fun <- function(x) pmin(x, 1 - x)

  n_y <- length(y)
  if (is.null(times)) {
    times <- seq_len(n_y) / n_y
  }
  k_n <- max(3L, as.integer(theta * sqrt(n_y)))
  l_n <- max(3L, as.integer(theta * (n_y^0.6)))
  weight <- g_fun((1:(k_n - 1L)) / k_n)
  y_diff <- diff(y)

  if (is.finite(clip_multiplier) && is.finite(noise_sd)) {
    y_diff[abs(y_diff) > clip_multiplier * noise_sd] <- 0
  }

  n_blocks <- n_y - l_n
  y_pre <- numeric(n_blocks)
  omegahat <- numeric(n_blocks)

  for (i in seq_len(n_blocks)) {
    y_pre[i] <- sum(weight * y_diff[i:(i + k_n - 2L)])
    omegahat[i] <- sum(y_diff[i:(i + l_n - 1L)]^2) / (2 * l_n)
  }

  delta_n <- max(times) / n_y
  phi2 <- sum(weight^2) / k_n
  phi1 <- 1
  upper_i <- as.integer((n_y - l_n) / k_n)
  idx_vec <- seq.int(1L, by = k_n, length.out = upper_i + 1L)

  list(
    y_pre = y_pre,
    omegahat = omegahat,
    k_n = k_n,
    l_n = l_n,
    delta_n = delta_n,
    phi2 = phi2,
    phi1 = phi1,
    idx_vec = idx_vec,
    theta = theta
  )
}

evaluate_single_preaverage <- function(u, prep) {
  vals <-
    cos(sqrt(2 * u) * prep$y_pre[prep$idx_vec] /
      sqrt(prep$phi2 * prep$delta_n * prep$k_n)) *
    exp(u * prep$phi1 * prep$omegahat[prep$idx_vec] / (prep$phi2 * prep$theta^2))

  prep$k_n * prep$delta_n * sum(vals)
}

prepare_joint_preaverage <- function(y1,
                                     y2,
                                     theta,
                                     times = NULL,
                                     noise_sd1 = NA_real_,
                                     noise_sd2 = NA_real_,
                                     clip_multiplier = Inf) {
  g_fun <- function(x) pmin(x, 1 - x)

  n_y <- length(y1)
  if (is.null(times)) {
    times <- seq_len(n_y) / n_y
  }
  k_n <- max(3L, as.integer(theta * sqrt(n_y)))
  l_n <- max(3L, as.integer(theta * (n_y^0.6)))
  weight <- g_fun((1:(k_n - 1L)) / k_n)

  y1_diff <- diff(y1)
  y2_diff <- diff(y2)

  if (is.finite(clip_multiplier) && is.finite(noise_sd1)) {
    y1_diff[abs(y1_diff) > clip_multiplier * noise_sd1] <- 0
  }
  if (is.finite(clip_multiplier) && is.finite(noise_sd2)) {
    y2_diff[abs(y2_diff) > clip_multiplier * noise_sd2] <- 0
  }

  n_blocks <- n_y - l_n
  y1_preave <- numeric(n_blocks)
  y2_preave <- numeric(n_blocks)
  omegahat1_sq <- numeric(n_blocks)
  omegahat2_sq <- numeric(n_blocks)

  for (i in seq_len(n_blocks)) {
    y1_preave[i] <- sum(weight * y1_diff[i:(i + k_n - 2L)])
    y2_preave[i] <- sum(weight * y2_diff[i:(i + k_n - 2L)])
    omegahat1_sq[i] <- sum(y1_diff[i:(i + l_n - 1L)]^2) / (2 * l_n)
    omegahat2_sq[i] <- sum(y2_diff[i:(i + l_n - 1L)]^2) / (2 * l_n)
  }

  delta_n <- max(times) / n_y
  phi2 <- sum(weight^2) / k_n
  phi1 <- 1
  upper_i <- as.integer((n_y - k_n - l_n) / (2 * k_n))
  idx1 <- 2L * seq_len(upper_i) * k_n
  idx2 <- idx1 + k_n

  list(
    y1preave = y1_preave,
    y2preave = y2_preave,
    omegahat1 = omegahat1_sq,
    omegahat2 = omegahat2_sq,
    k_n = k_n,
    l_n = l_n,
    delta_n = delta_n,
    phi2 = phi2,
    phi1 = phi1,
    idx1 = idx1,
    idx2 = idx2,
    theta = theta
  )
}

evaluate_joint_preaverage <- function(u, v, prep) {
  cosine_terms <-
    cos(sqrt(2 * u) * prep$y1preave[prep$idx1] /
      sqrt(prep$phi2 * prep$delta_n * prep$k_n)) *
    cos(sqrt(2 * v) * prep$y2preave[prep$idx2] /
      sqrt(prep$phi2 * prep$delta_n * prep$k_n))

  bias_terms <-
    exp(u * prep$phi1 * prep$omegahat1[prep$idx1] / (prep$phi2 * prep$theta^2)) *
    exp(v * prep$phi1 * prep$omegahat2[prep$idx2] / (prep$phi2 * prep$theta^2))

  2 * prep$k_n * prep$delta_n * sum(cosine_terms * bias_terms)
}

hac_variance <- function(x,
                         bandwidth = max(1L, as.integer(floor(length(x)^(1 / 4))))) {
  x_centered <- x - mean(x)
  n_x <- length(x_centered)
  gamma0 <- sum(x_centered * x_centered) / n_x
  out <- gamma0

  if (bandwidth >= 1L) {
    for (lag in seq_len(bandwidth)) {
      gamma_lag <- sum(x_centered[(lag + 1L):n_x] * x_centered[1L:(n_x - lag)]) / n_x
      weight <- 1 - lag / (bandwidth + 1)
      out <- out + 2 * weight * gamma_lag
    }
  }

  max(out, .Machine$double.eps)
}

compute_dependence_test <- function(L12, L1, L2, bandwidth = NULL) {
  stopifnot(length(L12) == length(L1), length(L1) == length(L2))
  T_days <- length(L12)
  if (is.null(bandwidth)) {
    bandwidth <- max(1L, as.integer(floor(T_days^(1 / 4))))
  }

  mean_L12 <- mean(L12)
  mean_L1 <- mean(L1)
  mean_L2 <- mean(L2)
  diff_mean <- mean_L12 - mean_L1 * mean_L2
  psi <- (L12 - mean_L12) -
    mean_L2 * (L1 - mean_L1) -
    mean_L1 * (L2 - mean_L2)
  eta2_hat <- hac_variance(psi, bandwidth = bandwidth)
  test_stat <- sqrt(T_days) * diff_mean / sqrt(eta2_hat)
  p_value <- 2 * (1 - pnorm(abs(test_stat)))

  data.frame(
    T_days = T_days,
    bandwidth = bandwidth,
    mean_L12 = mean_L12,
    mean_L1 = mean_L1,
    mean_L2 = mean_L2,
    diff_mean = diff_mean,
    eta2_hat = eta2_hat,
    test_stat = test_stat,
    p_value = p_value,
    reject05 = as.integer(abs(test_stat) > qnorm(0.975)),
    reject10 = as.integer(abs(test_stat) > qnorm(0.95))
  )
}

render_pairwise_matrix_tex <- function(matrix_df, value_col, file_path, digits = 4) {
  tickers <- colnames(matrix_df)[-1L]
  col_spec <- paste0("l", paste(rep("c", length(tickers)), collapse = ""))
  fmt <- function(x) {
    ifelse(
      is.na(x),
      "--",
      ifelse(abs(x) < 10^(-digits), paste0("<", sprintf(paste0("%.", digits, "f"), 10^(-digits))), sprintf(paste0("%.", digits, "f"), x))
    )
  }

  lines <- c(
    paste0("\\begin{tabular}{", col_spec, "}"),
    "\\toprule",
    paste0(paste(c("", tickers), collapse = " & "), " \\\\"),
    "\\midrule"
  )

  for (i in seq_len(nrow(matrix_df))) {
    row_vals <- fmt(as.numeric(matrix_df[i, tickers, drop = TRUE]))
    lines <- c(lines, paste0(paste(c(matrix_df$ticker[i], row_vals), collapse = " & "), " \\\\"))
  }

  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, con = file_path)
}
