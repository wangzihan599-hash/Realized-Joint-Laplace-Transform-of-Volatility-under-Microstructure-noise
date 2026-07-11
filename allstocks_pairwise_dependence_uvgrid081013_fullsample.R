# Full-sample all-stock pairwise dependence tests for the UV grid {0.8, 1, 1.3}^2.
# This script intentionally skips rolling-window calculations.

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
  raw_value <- Sys.getenv("ALLSTOCKS_UV081013_ROOTS", unset = "")
  if (raw_value == "") {
    return(roots)
  }

  items <- strsplit(raw_value, ",")[[1]]
  out <- list()
  for (item in items) {
    parts <- strsplit(trimws(item), ":", fixed = TRUE)[[1]]
    if (length(parts) != 2L) {
      stop("ALLSTOCKS_UV081013_ROOTS entries must have format sector:path")
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

period_label <- function(period) {
  ifelse(period == "2007_2008", "2007--2008", "2016--2017")
}

pair_type_label <- function(pair_type) {
  switch(
    pair_type,
    cross_sector = "Cross-sector",
    within_energy = "Within energy",
    within_health = "Within health",
    pair_type
  )
}

format_pvalue <- function(x) {
  ifelse(
    is.na(x),
    "",
    ifelse(x < 0.01, "$<0.01$", sprintf("%.2f", x))
  )
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

  col_spec <- paste0("l|", paste(rep("c", length(tickers)), collapse = ""))
  lines <- c(
    "\\begin{table}[htbp!]",
    "\\centering",
    "\\scriptsize",
    paste0("\\caption{\\textcolor{blue}{", caption, "}}"),
    paste0("\\label{", label, "}"),
    "{\\color{blue}",
    paste0("\\begin{tabular}{", col_spec, "}"),
    "\\toprule",
    paste0(
      "Sector & ",
      paste(
        sprintf("\\multicolumn{%d}{c}{%s}", sector_runs$lengths, sector_labels),
        collapse = " & "
      ),
      " \\\\"
    ),
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

render_group_table_tex <- function(full_summary,
                                   uv_grid,
                                   out_dir) {
  uv_names <- paste0("uv_", sanitize_value(uv_grid$u), "_", sanitize_value(uv_grid$v))
  uv_labels <- paste0("$(u,v)=(", uv_grid$u, ",", uv_grid$v, ")$")
  row_order <- do.call(
    rbind,
    lapply(c("2007_2008", "2016_2017"), function(period) {
      data.frame(
        period = period,
        pair_type = c("cross_sector", "within_energy", "within_health"),
        stringsAsFactors = FALSE
      )
    })
  )

  rows <- lapply(seq_len(nrow(row_order)), function(i) {
    period <- row_order$period[i]
    ptype <- row_order$pair_type[i]
    row <- data.frame(
      Period = period_label(period),
      `Pair group` = pair_type_label(ptype),
      `No. pairs` = full_summary$n_pairs[
        full_summary$period == period &
          full_summary$pair_type == ptype
      ][1],
      check.names = FALSE
    )
    for (k in seq_len(nrow(uv_grid))) {
      sub <- full_summary[
        full_summary$period == period &
          full_summary$pair_type == ptype &
          abs(full_summary$u - uv_grid$u[k]) < 1e-12 &
          abs(full_summary$v - uv_grid$v[k]) < 1e-12,
      ]
      row[[uv_names[k]]] <- paste0(sub$reject10, "/", sub$n_pairs)
    }
    row
  })
  full_table <- do.call(rbind, rows)
  write.csv(full_table, file.path(out_dir, "table7_full_sample_uv081013_10pct.csv"), row.names = FALSE)

  align <- paste0("lll", paste(rep("c", length(uv_names)), collapse = ""))
  lines <- c(
    "\\begin{table}[htbp!]",
    "\\centering",
    "\\caption{\\textcolor{blue}{Full-sample pairwise volatility-dependence tests for ten stocks}}",
    "\\label{tab:emp_allstocks_fullsample_uv081013}",
    "{\\color{blue}",
    "\\scriptsize",
    "\\resizebox{0.98\\textwidth}{!}{%",
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste(c("Period", "Pair group", "No. pairs", uv_labels), collapse = " & "),
    "\\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(full_table))) {
    values <- c(
      full_table$Period[i],
      full_table[["Pair group"]][i],
      full_table[["No. pairs"]][i],
      as.character(full_table[i, uv_names])
    )
    lines <- c(lines, paste(values, collapse = " & "), "\\\\")
  }
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}}",
    "}",
    "\\vspace{0.1cm}",
    "\\parbox{0.96\\textwidth}{\\color{blue}\\scriptsize Note: Each entry reports the number of pairwise tests rejecting the null of volatility independence at the 10\\% significance level. The tests use the original, uncorrected realized joint and marginal Laplace transforms.}",
    "\\end{table}",
    ""
  )
  writeLines(lines, file.path(out_dir, "table7_full_sample_uv081013_10pct.tex"))
}

out_dir <- file.path(
  script_dir,
  Sys.getenv("ALLSTOCKS_UV081013_OUTDIR", "allstocks_dependence_outputs_uv081013_raw")
)
matrix_dir <- file.path(out_dir, "pairwise_matrices")
latex_dir <- file.path(out_dir, "latex_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(latex_dir, recursive = TRUE, showWarnings = FALSE)

uv_values <- c(0.8, 1, 1.3)
uv_grid <- expand.grid(u = uv_values, v = uv_values)
uv_grid <- uv_grid[order(uv_grid$u, uv_grid$v), ]
tickers <- parse_char_grid(
  "ALLSTOCKS_UV081013_TICKERS",
  c("MRO", "NOV", "OXY", "RIG", "VLO", "ESRX", "GILD", "HOLX", "MCK", "PFE")
)
periods_to_run <- parse_char_grid("ALLSTOCKS_UV081013_PERIODS", c("2007_2008", "2016_2017"))
theta_value <- as.numeric(Sys.getenv("ALLSTOCKS_UV081013_THETA", as.character(1 / 3)))
clip_multiplier <- as.numeric(Sys.getenv("ALLSTOCKS_UV081013_CLIP", "Inf"))
mc_cores <- as.integer(Sys.getenv("ALLSTOCKS_UV081013_CORES", "4"))

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
    estimate = values
  )
}

message("Computing daily marginal transforms...")
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
    data.frame(period = period, date = common_dates_list[[period]], stringsAsFactors = FALSE)
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
      estimate = values
    )
  }
  do.call(rbind, out_rows)
}

message("Computing daily joint transforms...")
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
          end_date = max(panel$date)
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
    render_pairwise_matrix_tex(p_matrix, "p_value", file.path(latex_dir, paste0(base_name, "_pvalue_simple.tex")))
    render_sector_pvalue_tex(
      p_matrix,
      tickers,
      sector_lookup,
      paste0(
        "p-value of the independence test for volatilities, ",
        period_label(period),
        ", $(u,v)=(",
        u_value,
        ",",
        v_value,
        ")$"
      ),
      paste0("tab:pvalue_", period, "_u", sanitize_value(u_value), "_v", sanitize_value(v_value)),
      file.path(latex_dir, paste0(base_name, "_pvalue_sector.tex"))
    )
    render_pairwise_matrix_tex(stat_matrix, "test_stat", file.path(latex_dir, paste0(base_name, "_teststat.tex")))
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
render_group_table_tex(full_group_summary, uv_grid, out_dir)

writeLines(
  c(
    "Design summary",
    "--------------",
    paste0("Tickers: ", paste(tickers, collapse = ", ")),
    paste0("Periods: ", paste(periods_to_run, collapse = ", ")),
    paste0("UV grid: ", paste(apply(uv_grid, 1, function(x) paste0("(", x[1], ",", x[2], ")")), collapse = ", ")),
    paste0("Theta: ", theta_value),
    paste0("Clip multiplier: ", ifelse(is.finite(clip_multiplier), clip_multiplier, "disabled")),
    "Full-sample analysis only; rolling-window tests are skipped in this script."
  ),
  con = file.path(out_dir, "allstocks_dependence_analysis_notes.txt")
)

message("Completed full-sample all-stock pairwise dependence analysis.")
message("Outputs saved under: ", out_dir)
