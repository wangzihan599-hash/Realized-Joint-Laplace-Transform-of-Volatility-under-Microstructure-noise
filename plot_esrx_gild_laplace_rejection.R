get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
input_dir <- file.path(
  script_dir,
  Sys.getenv("ALLSTOCKS_DEP_OUTDIR", "allstocks_dependence_outputs")
)
figure_dir <- file.path(input_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

u_value <- as.numeric(Sys.getenv("PLOT_U", "1"))
v_value <- as.numeric(Sys.getenv("PLOT_V", "1"))
alpha_value <- as.numeric(Sys.getenv("PLOT_ALPHA", "0.10"))
asset1 <- Sys.getenv("PLOT_ASSET1", "ESRX")
asset2 <- Sys.getenv("PLOT_ASSET2", "GILD")

joint_df <- read.csv(file.path(input_dir, "allstocks_daily_joint_transforms.csv"),
                     stringsAsFactors = FALSE)
marginal_df <- read.csv(file.path(input_dir, "allstocks_daily_marginals.csv"),
                        stringsAsFactors = FALSE)
rolling_df <- read.csv(file.path(input_dir, "allstocks_rolling_window_results.csv"),
                       stringsAsFactors = FALSE)
full_df <- read.csv(file.path(input_dir, "allstocks_full_sample_results.csv"),
                    stringsAsFactors = FALSE)

fmt_p <- function(p) {
  ifelse(p < 0.0001, "<0.0001", sprintf("%.4f", p))
}

sanitize_value <- function(x) {
  sub("\\.$", "", sub("0+$", "", sprintf("%.4f", x)))
}

period_label <- function(period) {
  if (period == "2007_2008") {
    return("2007--2008")
  }
  if (period == "2016_2017") {
    return("2016--2017")
  }
  period
}

get_pair_rows <- function(df, period, keep_cols = names(df)) {
  out <- df[
    df$period == period &
      df$ticker1 == asset1 &
      df$ticker2 == asset2 &
      abs(df$u - u_value) < 1e-12 &
      abs(df$v - v_value) < 1e-12,
    keep_cols,
    drop = FALSE
  ]
  if (nrow(out) == 0L) {
    out <- df[
      df$period == period &
        df$ticker1 == asset2 &
        df$ticker2 == asset1 &
        abs(df$u - v_value) < 1e-12 &
        abs(df$v - u_value) < 1e-12,
      keep_cols,
      drop = FALSE
    ]
  }
  out
}

build_daily_panel <- function(period) {
  joint_sub <- get_pair_rows(
    joint_df,
    period,
    c("date", "estimate")
  )
  if (nrow(joint_sub) == 0L) {
    stop("No joint transform rows found for ", period, ", ",
         asset1, "-", asset2, ", u=", u_value, ", v=", v_value)
  }
  names(joint_sub)[2] <- "L12"

  l1_sub <- marginal_df[
    marginal_df$period == period &
      marginal_df$ticker == asset1 &
      abs(marginal_df$param - u_value) < 1e-12,
    c("date", "estimate")
  ]
  names(l1_sub)[2] <- "L1"

  l2_sub <- marginal_df[
    marginal_df$period == period &
      marginal_df$ticker == asset2 &
      abs(marginal_df$param - v_value) < 1e-12,
    c("date", "estimate")
  ]
  names(l2_sub)[2] <- "L2"

  panel <- merge(joint_sub, l1_sub, by = "date", all = FALSE)
  panel <- merge(panel, l2_sub, by = "date", all = FALSE)
  panel$date <- as.Date(panel$date)
  panel[order(panel$date), ]
}

draw_one_period <- function(period, output_prefix) {
  daily_panel <- build_daily_panel(period)
  rolling_sub <- get_pair_rows(rolling_df, period)
  full_sub <- get_pair_rows(full_df, period)

  if (nrow(rolling_sub) == 0L) {
    stop("No rolling-window rows found for ", period, ", ",
         asset1, "-", asset2, ", u=", u_value, ", v=", v_value)
  }
  if (nrow(full_sub) == 0L) {
    stop("No full-sample rows found for ", period, ", ",
         asset1, "-", asset2, ", u=", u_value, ", v=", v_value)
  }

  rolling_sub$window_start <- as.Date(rolling_sub$window_start)
  rolling_sub <- rolling_sub[order(rolling_sub$window_start), ]

  transform_values <- c(daily_panel$L12, daily_panel$L1, daily_panel$L2)
  y_pad <- diff(range(transform_values)) * 0.08
  if (!is.finite(y_pad) || y_pad == 0) {
    y_pad <- 0.001
  }
  y_lim <- range(transform_values) + c(-y_pad, y_pad)
  z90 <- qnorm(0.95)
  rolling_sub$diff_se <- sqrt(pmax(rolling_sub$eta2_hat, 0) / rolling_sub$T_days)
  rolling_sub$diff_low90 <- rolling_sub$diff_mean - z90 * rolling_sub$diff_se
  rolling_sub$diff_high90 <- rolling_sub$diff_mean + z90 * rolling_sub$diff_se

  diff_values <- c(
    rolling_sub$diff_mean,
    rolling_sub$diff_low90,
    rolling_sub$diff_high90,
    0
  )
  diff_pad <- diff(range(diff_values, na.rm = TRUE)) * 0.12
  if (!is.finite(diff_pad) || diff_pad == 0) {
    diff_pad <- max(abs(diff_values), na.rm = TRUE) * 0.12
  }
  if (!is.finite(diff_pad) || diff_pad == 0) {
    diff_pad <- 1e-6
  }
  diff_lim <- range(diff_values, na.rm = TRUE) + c(-diff_pad, diff_pad)

  rolling_sub$window_end <- as.Date(rolling_sub$window_end)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  layout(matrix(c(1, 2), ncol = 1), heights = c(2.15, 1.25))
  par(oma = c(0.6, 0, 0, 0))

  par(mar = c(1.4, 5.1, 4.3, 2.0))
  plot(
    daily_panel$date, daily_panel$L12,
    type = "n",
    ylim = y_lim,
    xlab = "",
    ylab = "Laplace transform",
    xaxt = "n",
    main = paste0(asset1, "-", asset2, " realized Laplace transforms, ",
                  period_label(period))
  )
  grid(col = "gray88", lty = "dotted")
  lines(daily_panel$date, daily_panel$L12, col = "#1B9E77", lwd = 2.2)
  lines(daily_panel$date, daily_panel$L1, col = "#386CB0", lwd = 2.0)
  lines(daily_panel$date, daily_panel$L2, col = "#F0027F", lwd = 2.0)
  legend_labels <- as.expression(list(
    substitute(paste("Joint: ", hat(U)["12,t"], "(u,v)")),
    substitute(paste(ASSET, " marginal: ", hat(U)["1,t"], "(u,v)"),
               list(ASSET = asset1)),
    substitute(paste(ASSET, " marginal: ", hat(U)["2,t"], "(u,v)"),
               list(ASSET = asset2))
  ))
  legend(
    "bottomleft",
    legend = legend_labels,
    col = c("#1B9E77", "#386CB0", "#F0027F"),
    lwd = c(2.2, 2.0, 2.0),
    bty = "n",
    cex = 0.86
  )

  par(mar = c(4.2, 5.1, 0.7, 2.0))
  reject_col <- ifelse(rolling_sub$reject10 == 1L, "#D73027", "gray55")
  plot(
    rolling_sub$window_end, rolling_sub$diff_mean,
    type = "n",
    ylim = diff_lim,
    xlim = range(daily_panel$date, na.rm = TRUE),
    xlab = "Rolling-window end date",
    ylab = expression(hat(U)["12,t"] - hat(U)["1,t"] %*% hat(U)["2,t"])
  )
  grid(col = "gray88", lty = "dotted")
  polygon(
    c(rolling_sub$window_end, rev(rolling_sub$window_end)),
    c(rolling_sub$diff_low90, rev(rolling_sub$diff_high90)),
    border = NA,
    col = adjustcolor("#80B1D3", alpha.f = 0.32)
  )
  abline(h = 0, col = "gray30", lwd = 1.2, lty = 2)
  lines(rolling_sub$window_end, rolling_sub$diff_mean, col = "#1B9E77", lwd = 1.4)
  points(
    rolling_sub$window_end, rolling_sub$diff_mean,
    pch = 16,
    cex = 0.72,
    col = reject_col
  )
  legend(
    "topleft",
    legend = c("Estimate", "90% CI", "Zero", "Reject H0 at 10%", "Do not reject H0"),
    col = c("#1B9E77", adjustcolor("#80B1D3", alpha.f = 0.65), "gray30", "#D73027", "gray55"),
    pch = c(NA, 15, NA, 16, 16),
    lty = c(1, NA, 2, NA, NA),
    lwd = c(1.4, NA, 1.2, NA, NA),
    bty = "n",
    cex = 0.76
  )
}

save_period_plot <- function(period) {
  pair_name <- paste0(tolower(asset1), "_", tolower(asset2))
  base_name <- paste0(
    pair_name,
    "_laplace_rejection_",
    period,
    "_u", sanitize_value(u_value),
    "_v", sanitize_value(v_value)
  )
  png_path <- file.path(figure_dir, paste0(base_name, ".png"))
  pdf_path <- file.path(figure_dir, paste0(base_name, ".pdf"))

  png(png_path, width = 1800, height = 1350, res = 200)
  draw_one_period(period, base_name)
  dev.off()

  pdf(pdf_path, width = 9, height = 6.8)
  draw_one_period(period, base_name)
  dev.off()

  data.frame(period = period, png = png_path, pdf = pdf_path, stringsAsFactors = FALSE)
}

outputs <- do.call(rbind, lapply(c("2007_2008", "2016_2017"), save_period_plot))
print(outputs)
