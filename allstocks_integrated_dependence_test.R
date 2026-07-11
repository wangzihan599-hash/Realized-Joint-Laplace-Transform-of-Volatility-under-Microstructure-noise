# Integrated full-sample all-stock volatility-dependence test.
# The test aggregates G_{n,T}(u,v)^2 over K=[0.8,1.3]^2 using a 10x10 grid.

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
source(file.path(script_dir, "energy_dependence_utils.R"))

parse_char_grid <- function(env_name, default_values) {
  raw_value <- Sys.getenv(env_name, unset = "")
  if (raw_value == "") {
    return(default_values)
  }
  trimws(strsplit(raw_value, ",")[[1]])
}

parse_sector_roots <- function() {
  roots <- list(
    energy = file.path(script_dir, "energy-onesecond"),
    health = file.path(script_dir, "health-onesecond")
  )
  raw_value <- Sys.getenv("ALLSTOCKS_INT_ROOTS", unset = "")
  if (raw_value == "") {
    return(roots)
  }

  items <- strsplit(raw_value, ",")[[1]]
  out <- list()
  for (item in items) {
    parts <- strsplit(trimws(item), ":", fixed = TRUE)[[1]]
    if (length(parts) != 2L) {
      stop("ALLSTOCKS_INT_ROOTS entries must have format sector:path")
    }
    out[[parts[1]]] <- parts[2]
  }
  out
}

parse_folder_ticker <- function(folder) {
  if (grepl("data", folder, fixed = TRUE)) {
    return(sub("data.*$", "", folder))
  }
  sub("-.*$", "", folder)
}

discover_onesecond_files_multi <- function(root_map) {
  rows <- list()
  row_id <- 1L
  for (sector in names(root_map)) {
    root_dir <- normalizePath(root_map[[sector]], winslash = "/", mustWork = TRUE)
    files <- list.files(
      root_dir,
      pattern = "_onesecond\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )

    for (path in files) {
      folder <- basename(dirname(path))
      ticker <- parse_folder_ticker(folder)
      date_value <- strsplit(basename(path), "_", fixed = TRUE)[[1]][2L]
      period <- if (date_value < "2010-01-01") "2007_2008" else "2016_2017"
      rows[[row_id]] <- data.frame(
        path = normalizePath(path, winslash = "/"),
        folder = folder,
        sector = sector,
        ticker = ticker,
        date = date_value,
        period = period,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }
  do.call(rbind, rows)
}

sanitize_value <- function(x) {
  sub("\\.$", "", sub("0+$", "", sprintf("%.4f", x)))
}

period_label <- function(period) {
  ifelse(period == "2007_2008", "2007--2008", "2016--2017")
}

pair_type <- function(ticker1, ticker2, sector_lookup) {
  s1 <- sector_lookup[[ticker1]]
  s2 <- sector_lookup[[ticker2]]
  if (identical(s1, s2)) {
    paste0("within_", s1)
  } else {
    "cross_sector"
  }
}

format_pvalue <- function(x) {
  ifelse(
    is.na(x),
    "",
    ifelse(x < 0.1, "$<0.1$", sprintf("%.2f", x))
  )
}

matrix_from_long <- function(long_df, tickers, value_col) {
  out <- data.frame(ticker = tickers, stringsAsFactors = FALSE)
  for (ticker in tickers) {
    out[[ticker]] <- NA_real_
  }
  for (i in seq_len(nrow(long_df))) {
    a <- long_df$ticker1[i]
    b <- long_df$ticker2[i]
    out[out$ticker == a, b] <- long_df[[value_col]][i]
    out[out$ticker == b, a] <- long_df[[value_col]][i]
  }
  out
}

render_sector_pvalue_tex <- function(matrix_df,
                                     tickers,
                                     sector_lookup,
                                     caption,
                                     label,
                                     file_path) {
  sectors <- vapply(tickers, function(ticker) sector_lookup[[ticker]], character(1))
  sector_runs <- rle(sectors)
  sector_labels <- ifelse(sector_runs$values == "energy", "Energy", "Health")
  col_blocks <- vapply(sector_runs$lengths, function(n) paste(rep("c", n), collapse = ""), character(1))
  col_spec <- paste0("l|", paste(col_blocks, collapse = "|"))
  sector_header <- vapply(seq_along(sector_runs$lengths), function(i) {
    align <- if (i < length(sector_runs$lengths)) "c|" else "c"
    sprintf("\\multicolumn{%d}{%s}{%s}", sector_runs$lengths[i], align, sector_labels[i])
  }, character(1))

  lines <- c(
    "\\begin{table}[htbp!]",
    "\\centering",
    "\\scriptsize",
    paste0("\\caption{\\textcolor{red}{", caption, "}}"),
    paste0("\\label{", label, "}"),
    "{\\color{red}",
    paste0("\\begin{tabular}{", col_spec, "}"),
    "\\toprule",
    paste0("Sector & ", paste(sector_header, collapse = " & "), " \\\\"),
    paste0("Symbol & ", paste(tickers, collapse = " & "), " \\\\"),
    "\\midrule"
  )

  for (i in seq_along(tickers)) {
    row_ticker <- tickers[i]
    values <- character(length(tickers))
    for (j in seq_along(tickers)) {
      if (j <= i) {
        values[j] <- ""
      } else {
        values[j] <- format_pvalue(matrix_df[matrix_df$ticker == row_ticker, tickers[j]])
      }
    }
    lines <- c(lines, paste0(row_ticker, " & ", paste(values, collapse = " & "), " \\\\"))
    if (i %in% cumsum(sector_runs$lengths)[-length(sector_runs$lengths)]) {
      lines <- c(lines, "\\midrule")
    }
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    "\\end{table}",
    ""
  )
  writeLines(lines, file_path)
}

evaluate_joint_from_single_preps <- function(u, v, prep1, prep2) {
  k_n <- prep1$k_n
  n_y <- length(prep1$y_pre) + prep1$l_n
  upper_i <- as.integer((n_y - k_n - prep1$l_n) / (2 * k_n))
  idx1 <- 2L * seq_len(upper_i) * k_n
  idx2 <- idx1 + k_n

  cosine_terms <-
    cos(sqrt(2 * u) * prep1$y_pre[idx1] /
      sqrt(prep1$phi2 * prep1$delta_n * prep1$k_n)) *
    cos(sqrt(2 * v) * prep2$y_pre[idx2] /
      sqrt(prep2$phi2 * prep2$delta_n * prep2$k_n))

  bias_terms <-
    exp(u * prep1$phi1 * prep1$omegahat[idx1] / (prep1$phi2 * prep1$theta^2)) *
    exp(v * prep2$phi1 * prep2$omegahat[idx2] / (prep2$phi2 * prep2$theta^2))

  2 * k_n * prep1$delta_n * sum(cosine_terms * bias_terms)
}

hac_covariance_matrix <- function(psi_matrix, bandwidth) {
  x <- scale(psi_matrix, center = TRUE, scale = FALSE)
  n_t <- nrow(x)
  out <- crossprod(x) / n_t
  if (bandwidth >= 1L) {
    for (lag in seq_len(bandwidth)) {
      gamma_lag <- crossprod(x[(lag + 1L):n_t, , drop = FALSE], x[1L:(n_t - lag), , drop = FALSE]) / n_t
      weight <- 1 - lag / (bandwidth + 1)
      out <- out + weight * (gamma_lag + t(gamma_lag))
    }
  }
  (out + t(out)) / 2
}

simulate_mixture_pvalue <- function(observed_stat,
                                    eigenvalues,
                                    scale_factor,
                                    n_sim,
                                    chunk_size = 5000L) {
  eigenvalues <- pmax(Re(eigenvalues), 0)
  eigenvalues <- eigenvalues[eigenvalues > 1e-12]
  if (length(eigenvalues) == 0L) {
    return(NA_real_)
  }

  n_done <- 0L
  n_ge <- 0L
  while (n_done < n_sim) {
    b <- min(chunk_size, n_sim - n_done)
    draws <- matrix(rchisq(b * length(eigenvalues), df = 1), nrow = b)
    sim_values <- scale_factor * as.vector(draws %*% eigenvalues)
    n_ge <- n_ge + sum(sim_values >= observed_stat)
    n_done <- n_done + b
  }
  (n_ge + 1) / (n_sim + 1)
}

compute_integrated_test <- function(L12_mat,
                                    L1_mat,
                                    L2_mat,
                                    grid_df,
                                    area,
                                    bandwidth = NULL,
                                    n_sim = 50000L) {
  T_days <- nrow(L12_mat)
  M <- nrow(grid_df)
  if (is.null(bandwidth)) {
    bandwidth <- max(1L, as.integer(floor(T_days^(1 / 4))))
  }

  mean_L12 <- colMeans(L12_mat)
  mean_L1 <- colMeans(L1_mat)
  mean_L2 <- colMeans(L2_mat)
  diff_vec <- mean_L12 - mean_L1 * mean_L2
  G_vec <- sqrt(T_days) * diff_vec
  observed <- area / M * sum(G_vec^2)

  psi_mat <- L12_mat
  for (r in seq_len(M)) {
    psi_mat[, r] <- (L12_mat[, r] - mean_L12[r]) -
      mean_L2[r] * (L1_mat[, r] - mean_L1[r]) -
      mean_L1[r] * (L2_mat[, r] - mean_L2[r])
  }

  cov_hat <- hac_covariance_matrix(psi_mat, bandwidth = bandwidth)
  eig_vals <- eigen(cov_hat, symmetric = TRUE, only.values = TRUE)$values
  p_value <- simulate_mixture_pvalue(
    observed_stat = observed,
    eigenvalues = eig_vals,
    scale_factor = area / M,
    n_sim = n_sim
  )

  data.frame(
    T_days = T_days,
    bandwidth = bandwidth,
    grid_points = M,
    area = area,
    statistic = observed,
    p_value = p_value,
    reject05 = as.integer(!is.na(p_value) && p_value < 0.05),
    reject10 = as.integer(!is.na(p_value) && p_value < 0.10),
    min_eigen = min(eig_vals),
    positive_eigen = sum(eig_vals > 1e-12),
    stringsAsFactors = FALSE
  )
}

out_dir <- file.path(
  script_dir,
  Sys.getenv("ALLSTOCKS_INT_OUTDIR", "allstocks_integrated_outputs")
)
matrix_dir <- file.path(out_dir, "pairwise_matrices")
latex_dir <- file.path(out_dir, "latex_tables")
cache_dir <- file.path(out_dir, "cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(latex_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

tickers <- parse_char_grid(
  "ALLSTOCKS_INT_TICKERS",
  c("MRO", "NOV", "OXY", "RIG", "VLO", "ESRX", "GILD", "HOLX", "MCK", "PFE")
)
periods_to_run <- parse_char_grid("ALLSTOCKS_INT_PERIODS", c("2007_2008", "2016_2017"))
theta_value <- as.numeric(Sys.getenv("ALLSTOCKS_INT_THETA", as.character(1 / 3)))
clip_multiplier <- as.numeric(Sys.getenv("ALLSTOCKS_INT_CLIP", "Inf"))
mc_cores <- as.integer(Sys.getenv("ALLSTOCKS_INT_CORES", "4"))
n_sim <- as.integer(Sys.getenv("ALLSTOCKS_INT_NSIM", "50000"))
set.seed(as.integer(Sys.getenv("ALLSTOCKS_INT_SEED", "20260710")))

u_min <- as.numeric(Sys.getenv("ALLSTOCKS_INT_UMIN", "0.8"))
u_max <- as.numeric(Sys.getenv("ALLSTOCKS_INT_UMAX", "1.3"))
v_min <- as.numeric(Sys.getenv("ALLSTOCKS_INT_VMIN", "0.8"))
v_max <- as.numeric(Sys.getenv("ALLSTOCKS_INT_VMAX", "1.3"))
n_grid_side <- as.integer(Sys.getenv("ALLSTOCKS_INT_GRID_SIDE", "10"))
u_values <- seq(u_min, u_max, length.out = n_grid_side)
v_values <- seq(v_min, v_max, length.out = n_grid_side)
uv_grid <- expand.grid(u = u_values, v = v_values)
uv_grid <- uv_grid[order(uv_grid$u, uv_grid$v), ]
area_K <- (u_max - u_min) * (v_max - v_min)

root_map <- parse_sector_roots()
metadata <- discover_onesecond_files_multi(root_map)
metadata <- metadata[metadata$ticker %in% tickers & metadata$period %in% periods_to_run, ]
metadata <- metadata[order(metadata$period, metadata$ticker, metadata$date), ]
if (nrow(metadata) == 0L) {
  stop("No matching one-second files found.")
}

sector_lookup <- as.list(tapply(metadata$sector, metadata$ticker, function(x) unique(x)[1]))
missing_tickers <- setdiff(tickers, names(sector_lookup))
if (length(missing_tickers) > 0L) {
  stop("Missing tickers in one-second data: ", paste(missing_tickers, collapse = ", "))
}

common_dates_list <- list()
common_dates_rows <- list()
for (period in periods_to_run) {
  period_meta <- metadata[metadata$period == period, ]
  ticker_sets <- lapply(tickers, function(ticker) {
    sort(unique(period_meta$date[period_meta$ticker == ticker]))
  })
  names(ticker_sets) <- tickers
  common_dates <- sort(Reduce(intersect, ticker_sets))
  if (length(common_dates) == 0L) {
    stop("No common dates found for period ", period)
  }
  common_dates_list[[period]] <- common_dates
  common_dates_rows[[period]] <- data.frame(period = period, date = common_dates)
}
write.csv(do.call(rbind, common_dates_rows), file.path(out_dir, "allstocks_integrated_common_dates.csv"), row.names = FALSE)

path_lookup <- setNames(
  metadata$path,
  paste(metadata$period, metadata$ticker, metadata$date, sep = "|")
)
pairs <- combn(tickers, 2L, simplify = FALSE)
names(pairs) <- vapply(pairs, paste, collapse = "_", character(1))

compute_period_day <- function(task) {
  period <- task$period
  date <- task$date
  cache_path <- file.path(cache_dir, paste0("daily_", period, "_", date, ".rds"))
  if (file.exists(cache_path)) {
    return(readRDS(cache_path))
  }

  date_paths <- setNames(path_lookup[paste(period, tickers, date, sep = "|")], tickers)
  day_prices <- lapply(date_paths, load_onesecond_logprice)
  n_day <- length(day_prices[[1]])
  local_times <- seq_len(n_day) / n_day

  preps <- lapply(tickers, function(ticker) {
    prepare_single_preaverage(
      y = day_prices[[ticker]],
      theta = theta_value,
      times = local_times,
      clip_multiplier = clip_multiplier
    )
  })
  names(preps) <- tickers

  marginal_values <- lapply(tickers, function(ticker) {
    vapply(u_values, evaluate_single_preaverage, prep = preps[[ticker]], numeric(1))
  })
  names(marginal_values) <- tickers

  joint_values <- lapply(pairs, function(pair) {
    vapply(seq_len(nrow(uv_grid)), function(k) {
      evaluate_joint_from_single_preps(
        u = uv_grid$u[k],
        v = uv_grid$v[k],
        prep1 = preps[[pair[1]]],
        prep2 = preps[[pair[2]]]
      )
    }, numeric(1))
  })

  out <- list(
    period = period,
    date = date,
    marginal = marginal_values,
    joint = joint_values
  )
  saveRDS(out, cache_path)
  out
}

all_day_tasks <- do.call(
  rbind,
  lapply(periods_to_run, function(period) {
    data.frame(period = period, date = common_dates_list[[period]], stringsAsFactors = FALSE)
  })
)

message("Computing daily transforms for integrated test...")
day_results <- if (mc_cores > 1L) {
  parallel::mclapply(
    split(all_day_tasks, seq_len(nrow(all_day_tasks))),
    compute_period_day,
    mc.cores = mc_cores
  )
} else {
  lapply(split(all_day_tasks, seq_len(nrow(all_day_tasks))), compute_period_day)
}

names(day_results) <- paste(
  vapply(day_results, `[[`, character(1), "period"),
  vapply(day_results, `[[`, character(1), "date"),
  sep = "|"
)

build_integrated_panel <- function(period, pair) {
  dates <- common_dates_list[[period]]
  M <- nrow(uv_grid)
  L12_mat <- matrix(NA_real_, nrow = length(dates), ncol = M)
  L1_mat <- matrix(NA_real_, nrow = length(dates), ncol = M)
  L2_mat <- matrix(NA_real_, nrow = length(dates), ncol = M)
  pair_name <- paste(pair, collapse = "_")

  for (i in seq_along(dates)) {
    day_obj <- day_results[[paste(period, dates[i], sep = "|")]]
    L12_mat[i, ] <- day_obj$joint[[pair_name]]
    L1_u <- day_obj$marginal[[pair[1]]]
    L2_v <- day_obj$marginal[[pair[2]]]
    L1_mat[i, ] <- L1_u[match(uv_grid$u, u_values)]
    L2_mat[i, ] <- L2_v[match(uv_grid$v, v_values)]
  }

  list(L12 = L12_mat, L1 = L1_mat, L2 = L2_mat)
}

message("Computing integrated full-sample tests...")
result_rows <- list()
row_id <- 1L
for (period in periods_to_run) {
  for (pair in pairs) {
    panel <- build_integrated_panel(period, pair)
    stats <- compute_integrated_test(
      L12_mat = panel$L12,
      L1_mat = panel$L1,
      L2_mat = panel$L2,
      grid_df = uv_grid,
      area = area_K,
      n_sim = n_sim
    )
    result_rows[[row_id]] <- cbind(
      data.frame(
        period = period,
        ticker1 = pair[1],
        ticker2 = pair[2],
        pair_type = pair_type(pair[1], pair[2], sector_lookup),
        K = paste0("[", u_min, ",", u_max, "]x[", v_min, ",", v_max, "]"),
        start_date = min(common_dates_list[[period]]),
        end_date = max(common_dates_list[[period]]),
        stringsAsFactors = FALSE
      ),
      stats
    )
    row_id <- row_id + 1L
  }
}

integrated_results <- do.call(rbind, result_rows)
integrated_results <- integrated_results[order(integrated_results$period, integrated_results$ticker1, integrated_results$ticker2), ]
write.csv(integrated_results, file.path(out_dir, "allstocks_integrated_full_sample_results.csv"), row.names = FALSE)

message("Rendering integrated p-value matrices...")
for (period in periods_to_run) {
  sub <- integrated_results[integrated_results$period == period, ]
  p_matrix <- matrix_from_long(sub, tickers, "p_value")
  stat_matrix <- matrix_from_long(sub, tickers, "statistic")
  base_name <- paste0("allstocks_integrated_pairwise_", period)
  write.csv(p_matrix, file.path(matrix_dir, paste0(base_name, "_pvalue.csv")), row.names = FALSE)
  write.csv(stat_matrix, file.path(matrix_dir, paste0(base_name, "_statistic.csv")), row.names = FALSE)
  render_sector_pvalue_tex(
    p_matrix,
    tickers,
    sector_lookup,
    paste0(
      "$p$-value of the integrated independence test for volatilities, ",
      period_label(period),
      ", $(u,v)\\in[",
      u_min,
      ",",
      u_max,
      "]\\times[",
      v_min,
      ",",
      v_max,
      "]$"
    ),
    paste0("tab:integrated_pvalue_", period),
    file.path(latex_dir, paste0(base_name, "_pvalue.tex"))
  )
}

summary_by_period <- aggregate(
  cbind(reject05, reject10) ~ period + pair_type,
  data = integrated_results,
  FUN = sum
)
summary_counts <- aggregate(
  reject05 ~ period + pair_type,
  data = integrated_results,
  FUN = length
)
names(summary_counts)[names(summary_counts) == "reject05"] <- "n_pairs"
summary_by_period <- merge(summary_by_period, summary_counts, by = c("period", "pair_type"))
summary_by_period$reject05_share <- summary_by_period$reject05 / summary_by_period$n_pairs
summary_by_period$reject10_share <- summary_by_period$reject10 / summary_by_period$n_pairs
write.csv(summary_by_period, file.path(out_dir, "allstocks_integrated_group_summary.csv"), row.names = FALSE)

writeLines(
  c(
    "Integrated test design summary",
    "------------------------------",
    paste0("Tickers: ", paste(tickers, collapse = ", ")),
    paste0("Periods: ", paste(periods_to_run, collapse = ", ")),
    paste0("K: [", u_min, ",", u_max, "] x [", v_min, ",", v_max, "]"),
    paste0("Grid side: ", n_grid_side),
    paste0("Grid points: ", nrow(uv_grid)),
    paste0("Theta: ", theta_value),
    paste0("H_T: floor(T^(1/4))"),
    paste0("Mixture simulations: ", n_sim),
    paste0("Clip multiplier: ", ifelse(is.finite(clip_multiplier), clip_multiplier, "disabled"))
  ),
  con = file.path(out_dir, "allstocks_integrated_analysis_notes.txt")
)

message("Completed integrated full-sample all-stock dependence analysis.")
message("Outputs saved under: ", out_dir)
