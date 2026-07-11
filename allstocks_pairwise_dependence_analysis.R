get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
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

sanitize_value <- function(x) {
  sub("\\.$", "", sub("0+$", "", sprintf("%.4f", x)))
}

parse_sector_roots <- function() {
  roots <- list(
    energy = file.path(script_dir, "energy-onesecond"),
    health = file.path(script_dir, "health-onesecond")
  )
  raw_value <- Sys.getenv("ALLSTOCKS_DEP_ROOTS", unset = "")
  if (raw_value == "") {
    return(roots)
  }

  items <- strsplit(raw_value, ",")[[1]]
  out <- list()
  for (item in items) {
    parts <- strsplit(trimws(item), ":", fixed = TRUE)[[1]]
    if (length(parts) != 2L) {
      stop("ALLSTOCKS_DEP_ROOTS entries must have format sector:path")
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

pair_type <- function(ticker1, ticker2, sector_lookup) {
  s1 <- sector_lookup[[ticker1]]
  s2 <- sector_lookup[[ticker2]]
  if (identical(s1, s2)) {
    paste0("within_", s1)
  } else {
    "cross_sector"
  }
}

out_dir <- file.path(
  script_dir,
  Sys.getenv("ALLSTOCKS_DEP_OUTDIR", "allstocks_dependence_outputs")
)
matrix_dir <- file.path(out_dir, "pairwise_matrices")
latex_dir <- file.path(out_dir, "latex_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(latex_dir, recursive = TRUE, showWarnings = FALSE)

default_uv_grid <- data.frame(
  u = c(1, 1),
  v = c(1, 0.5)
)
uv_grid <- parse_uv_grid("ALLSTOCKS_DEP_UV_GRID", default_uv_grid)
tickers <- parse_char_grid(
  "ALLSTOCKS_DEP_TICKERS",
  c("MRO", "NOV", "OXY", "RIG", "VLO", "ESRX", "GILD", "HOLX", "MCK", "PFE")
)
periods_to_run <- parse_char_grid("ALLSTOCKS_DEP_PERIODS", c("2007_2008", "2016_2017"))
window_length <- as.integer(Sys.getenv("ALLSTOCKS_DEP_WINDOW", "44"))
theta_value <- as.numeric(Sys.getenv("ALLSTOCKS_DEP_THETA", as.character(1 / 3)))
clip_multiplier <- as.numeric(Sys.getenv("ALLSTOCKS_DEP_CLIP", "Inf"))
mc_cores <- as.integer(Sys.getenv("ALLSTOCKS_DEP_CORES", "4"))

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
  common_dates_rows[[period]] <- data.frame(
    period = period,
    date = common_dates,
    stringsAsFactors = FALSE
  )
}

common_dates_df <- do.call(rbind, common_dates_rows)
write.csv(common_dates_df, file.path(out_dir, "allstocks_common_dates.csv"), row.names = FALSE)

path_lookup <- setNames(
  metadata$path,
  paste(metadata$period, metadata$ticker, metadata$date, sep = "|")
)

param_values <- sort(unique(c(uv_grid$u, uv_grid$v)))
marginal_tasks <- do.call(
  rbind,
  lapply(periods_to_run, function(period) {
    expand.grid(
      period = period,
      ticker = tickers,
      date = common_dates_list[[period]],
      stringsAsFactors = FALSE
    )
  })
)
marginal_tasks$path <- path_lookup[paste(marginal_tasks$period, marginal_tasks$ticker, marginal_tasks$date, sep = "|")]

compute_marginal_day <- function(i) {
  row <- marginal_tasks[i, ]
  y <- load_onesecond_logprice(row$path)
  prep <- prepare_single_preaverage(
    y = y,
    theta = theta_value,
    times = seq_len(length(y)) / length(y),
    clip_multiplier = clip_multiplier
  )
  values <- vapply(param_values, function(param) evaluate_single_preaverage(param, prep), numeric(1))
  data.frame(
    period = row$period,
    ticker = row$ticker,
    date = row$date,
    param = param_values,
    estimate = values,
    stringsAsFactors = FALSE
  )
}

message("Computing daily marginal transforms for all stocks...")
marginal_list <- if (mc_cores > 1L) {
  parallel::mclapply(seq_len(nrow(marginal_tasks)), compute_marginal_day, mc.cores = mc_cores)
} else {
  lapply(seq_len(nrow(marginal_tasks)), compute_marginal_day)
}
marginal_df <- do.call(rbind, marginal_list)
marginal_df <- marginal_df[order(marginal_df$period, marginal_df$ticker, marginal_df$date, marginal_df$param), ]
write.csv(marginal_df, file.path(out_dir, "allstocks_daily_marginals.csv"), row.names = FALSE)

pairs <- combn(tickers, 2L, simplify = FALSE)
joint_tasks <- do.call(
  rbind,
  lapply(periods_to_run, function(period) {
    data.frame(
      period = period,
      date = common_dates_list[[period]],
      stringsAsFactors = FALSE
    )
  })
)

compute_joint_day <- function(i) {
  row <- joint_tasks[i, ]
  date_paths <- setNames(
    path_lookup[paste(row$period, tickers, row$date, sep = "|")],
    tickers
  )
  day_prices <- lapply(date_paths, load_onesecond_logprice)
  n_day <- length(day_prices[[1]])
  local_times <- seq_len(n_day) / n_day

  out_rows <- vector("list", length(pairs))
  for (j in seq_along(pairs)) {
    pair <- pairs[[j]]
    prep <- prepare_joint_preaverage(
      y1 = day_prices[[pair[1]]],
      y2 = day_prices[[pair[2]]],
      theta = theta_value,
      times = local_times,
      clip_multiplier = clip_multiplier
    )
    values <- vapply(seq_len(nrow(uv_grid)), function(k) {
      evaluate_joint_preaverage(uv_grid$u[k], uv_grid$v[k], prep)
    }, numeric(1))
    out_rows[[j]] <- data.frame(
      period = row$period,
      date = row$date,
      ticker1 = pair[1],
      ticker2 = pair[2],
      pair_type = pair_type(pair[1], pair[2], sector_lookup),
      u = uv_grid$u,
      v = uv_grid$v,
      estimate = values,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out_rows)
}

message("Computing daily joint transforms for all pairs...")
joint_list <- if (mc_cores > 1L) {
  parallel::mclapply(seq_len(nrow(joint_tasks)), compute_joint_day, mc.cores = mc_cores)
} else {
  lapply(seq_len(nrow(joint_tasks)), compute_joint_day)
}
joint_df <- do.call(rbind, joint_list)
joint_df <- joint_df[order(joint_df$period, joint_df$ticker1, joint_df$ticker2, joint_df$date, joint_df$u, joint_df$v), ]
write.csv(joint_df, file.path(out_dir, "allstocks_daily_joint_transforms.csv"), row.names = FALSE)

build_pair_panel <- function(period, ticker1, ticker2, u, v) {
  joint_sub <- joint_df[
    joint_df$period == period &
      joint_df$ticker1 == ticker1 &
      joint_df$ticker2 == ticker2 &
      joint_df$u == u &
      joint_df$v == v,
    c("date", "estimate")
  ]
  names(joint_sub)[2] <- "L12"

  l1_sub <- marginal_df[
    marginal_df$period == period &
      marginal_df$ticker == ticker1 &
      marginal_df$param == u,
    c("date", "estimate")
  ]
  names(l1_sub)[2] <- "L1"

  l2_sub <- marginal_df[
    marginal_df$period == period &
      marginal_df$ticker == ticker2 &
      marginal_df$param == v,
    c("date", "estimate")
  ]
  names(l2_sub)[2] <- "L2"

  out <- merge(joint_sub, l1_sub, by = "date", all = FALSE)
  out <- merge(out, l2_sub, by = "date", all = FALSE)
  out[order(out$date), ]
}

message("Computing full-sample statistics...")
full_rows <- list()
row_id <- 1L
for (period in periods_to_run) {
  for (pair in pairs) {
    ptype <- pair_type(pair[1], pair[2], sector_lookup)
    for (k in seq_len(nrow(uv_grid))) {
      panel <- build_pair_panel(period, pair[1], pair[2], uv_grid$u[k], uv_grid$v[k])
      stats <- compute_dependence_test(panel$L12, panel$L1, panel$L2)
      full_rows[[row_id]] <- cbind(
        data.frame(
          period = period,
          ticker1 = pair[1],
          ticker2 = pair[2],
          pair_type = ptype,
          u = uv_grid$u[k],
          v = uv_grid$v[k],
          start_date = min(panel$date),
          end_date = max(panel$date),
          stringsAsFactors = FALSE
        ),
        stats
      )
      row_id <- row_id + 1L
    }
  }
}
full_results <- do.call(rbind, full_rows)
full_results <- full_results[order(full_results$period, full_results$u, full_results$v, full_results$ticker1, full_results$ticker2), ]
write.csv(full_results, file.path(out_dir, "allstocks_full_sample_results.csv"), row.names = FALSE)

message("Rendering full-sample pairwise matrices...")
for (period in periods_to_run) {
  for (k in seq_len(nrow(uv_grid))) {
    u_value <- uv_grid$u[k]
    v_value <- uv_grid$v[k]
    subset_rows <- full_results[
      full_results$period == period &
        full_results$u == u_value &
        full_results$v == v_value,
    ]

    p_matrix <- matrix_from_long(subset_rows, tickers, "p_value")
    stat_matrix <- matrix_from_long(subset_rows, tickers, "test_stat")

    base_name <- paste0(
      "allstocks_pairwise_",
      period,
      "_u", sanitize_value(u_value),
      "_v", sanitize_value(v_value)
    )

    write.csv(p_matrix, file.path(matrix_dir, paste0(base_name, "_pvalue.csv")), row.names = FALSE)
    write.csv(stat_matrix, file.path(matrix_dir, paste0(base_name, "_teststat.csv")), row.names = FALSE)
    render_pairwise_matrix_tex(p_matrix, "p_value", file.path(latex_dir, paste0(base_name, "_pvalue.tex")))
    render_pairwise_matrix_tex(stat_matrix, "test_stat", file.path(latex_dir, paste0(base_name, "_teststat.tex")))
  }
}

message("Computing rolling-window statistics...")
rolling_rows <- list()
row_id <- 1L
for (period in periods_to_run) {
  for (pair in pairs) {
    ptype <- pair_type(pair[1], pair[2], sector_lookup)
    for (k in seq_len(nrow(uv_grid))) {
      panel <- build_pair_panel(period, pair[1], pair[2], uv_grid$u[k], uv_grid$v[k])
      n_obs <- nrow(panel)
      if (n_obs < window_length) {
        next
      }
      for (start_idx in seq_len(n_obs - window_length + 1L)) {
        end_idx <- start_idx + window_length - 1L
        window_panel <- panel[start_idx:end_idx, ]
        stats <- compute_dependence_test(window_panel$L12, window_panel$L1, window_panel$L2)
        rolling_rows[[row_id]] <- cbind(
          data.frame(
            period = period,
            ticker1 = pair[1],
            ticker2 = pair[2],
            pair_type = ptype,
            u = uv_grid$u[k],
            v = uv_grid$v[k],
            window_start = window_panel$date[1L],
            window_end = window_panel$date[nrow(window_panel)],
            window_index = start_idx,
            stringsAsFactors = FALSE
          ),
          stats
        )
        row_id <- row_id + 1L
      }
    }
  }
}
rolling_results <- do.call(rbind, rolling_rows)
rolling_results <- rolling_results[order(
  rolling_results$period,
  rolling_results$u,
  rolling_results$v,
  rolling_results$ticker1,
  rolling_results$ticker2,
  rolling_results$window_start
), ]
write.csv(rolling_results, file.path(out_dir, "allstocks_rolling_window_results.csv"), row.names = FALSE)

rolling_summary <- aggregate(
  cbind(reject05, reject10, p_value) ~ period + ticker1 + ticker2 + pair_type + u + v,
  data = rolling_results,
  FUN = function(x) c(mean = mean(x), min = min(x), median = median(x))
)

expand_summary <- function(df) {
  numeric_cols <- c("reject05", "reject10", "p_value")
  out <- df[, c("period", "ticker1", "ticker2", "pair_type", "u", "v")]
  for (col in numeric_cols) {
    out[[paste0(col, "_mean")]] <- df[[col]][, "mean"]
    out[[paste0(col, "_min")]] <- df[[col]][, "min"]
    out[[paste0(col, "_median")]] <- df[[col]][, "median"]
  }
  out
}

rolling_summary_df <- expand_summary(rolling_summary)
rolling_summary_df <- rolling_summary_df[order(
  rolling_summary_df$period,
  rolling_summary_df$u,
  rolling_summary_df$v,
  rolling_summary_df$ticker1,
  rolling_summary_df$ticker2
), ]
write.csv(rolling_summary_df, file.path(out_dir, "allstocks_rolling_window_summary.csv"), row.names = FALSE)

message("Rendering rolling rejection-share matrices...")
for (period in periods_to_run) {
  for (k in seq_len(nrow(uv_grid))) {
    u_value <- uv_grid$u[k]
    v_value <- uv_grid$v[k]
    subset_rows <- rolling_summary_df[
      rolling_summary_df$period == period &
        rolling_summary_df$u == u_value &
        rolling_summary_df$v == v_value,
    ]

    reject05_matrix <- matrix_from_long(subset_rows, tickers, "reject05_mean")
    reject10_matrix <- matrix_from_long(subset_rows, tickers, "reject10_mean")
    base_name <- paste0(
      "allstocks_rolling_rejectshare_",
      period,
      "_u", sanitize_value(u_value),
      "_v", sanitize_value(v_value)
    )

    write.csv(reject05_matrix, file.path(matrix_dir, paste0(base_name, "_05.csv")), row.names = FALSE)
    write.csv(reject10_matrix, file.path(matrix_dir, paste0(base_name, "_10.csv")), row.names = FALSE)
    render_pairwise_matrix_tex(reject05_matrix, "reject05_mean", file.path(latex_dir, paste0(base_name, "_05.tex")))
    render_pairwise_matrix_tex(reject10_matrix, "reject10_mean", file.path(latex_dir, paste0(base_name, "_10.tex")))
  }
}

full_group_summary <- aggregate(
  cbind(reject05, reject10) ~ period + pair_type + u + v,
  data = full_results,
  FUN = sum
)
full_group_counts <- aggregate(
  reject05 ~ period + pair_type + u + v,
  data = full_results,
  FUN = length
)
names(full_group_counts)[names(full_group_counts) == "reject05"] <- "n_pairs"
full_group_summary <- merge(full_group_summary, full_group_counts, by = c("period", "pair_type", "u", "v"))
full_group_summary$reject05_share <- full_group_summary$reject05 / full_group_summary$n_pairs
full_group_summary$reject10_share <- full_group_summary$reject10 / full_group_summary$n_pairs
write.csv(full_group_summary, file.path(out_dir, "allstocks_full_sample_group_summary.csv"), row.names = FALSE)

rolling_group_summary <- aggregate(
  cbind(reject05_mean, reject10_mean, p_value_median) ~ period + pair_type + u + v,
  data = rolling_summary_df,
  FUN = mean
)
write.csv(rolling_group_summary, file.path(out_dir, "allstocks_rolling_group_summary.csv"), row.names = FALSE)

notes_path <- file.path(out_dir, "allstocks_dependence_analysis_notes.txt")
writeLines(
  c(
    "Design summary",
    "--------------",
    paste0("Tickers: ", paste(tickers, collapse = ", ")),
    paste0("Periods: ", paste(periods_to_run, collapse = ", ")),
    paste0("UV grid: ", paste(apply(uv_grid, 1, function(x) paste0("(", x[1], ",", x[2], ")")), collapse = ", ")),
    paste0("Rolling window length: ", window_length, " trading days"),
    paste0("Theta: ", theta_value),
    paste0("Clip multiplier: ", ifelse(is.finite(clip_multiplier), clip_multiplier, "disabled")),
    "",
    "Implementation choices",
    "----------------------",
    "1. Combine the five energy stocks and the five health stocks into one 10-stock pairwise analysis.",
    "2. Analyze 2007--2008 and 2016--2017 separately.",
    "3. Use the common-date intersection across all 10 stocks within each period.",
    "4. Compute daily joint transforms with V_n(Y1,Y2,1,u,v).",
    "5. Compute daily marginal transforms with the univariate estimator from equation (24).",
    "6. Full-sample inference uses all common trading days in each period.",
    "7. Rolling-window inference uses 44-day windows.",
    "8. Two-sided p-values use the asymptotic N(0,1) approximation of the studentized statistic."
  ),
  con = notes_path
)

message("Completed all-stock pairwise dependence analysis.")
message("Outputs saved under: ", out_dir)
