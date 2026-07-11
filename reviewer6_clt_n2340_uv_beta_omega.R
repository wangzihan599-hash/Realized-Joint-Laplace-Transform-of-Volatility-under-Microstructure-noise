# CLT diagnostics for n = 2340 under several values of v, beta, and omega.

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

parse_numeric_grid <- function(env_name, default_values) {
  raw_value <- Sys.getenv(env_name, unset = "")
  if (raw_value == "") {
    return(default_values)
  }
  as.numeric(strsplit(raw_value, ",")[[1]])
}

script_dir <- get_script_dir()
source(file.path(script_dir, "reviewer_sim_utils.R"))

set.seed(as.integer(Sys.getenv("CLT_N2340_SEED", "20260625")))

out_dir <- file.path(script_dir, "reviewer_outputs")
tables_dir <- file.path(out_dir, "latex_tables")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
if (!dir.exists(tables_dir)) {
  dir.create(tables_dir, recursive = TRUE)
}

n_rep <- as.integer(Sys.getenv("CLT_N2340_REPS", "300"))
n_value <- as.integer(Sys.getenv("CLT_N2340_N", "2340"))
u_value <- as.numeric(Sys.getenv("CLT_N2340_U", "1"))
v_grid <- parse_numeric_grid("CLT_N2340_V_GRID", c(0.5, 1, 1.5))
beta_grid <- parse_numeric_grid("CLT_N2340_BETA_GRID", c(0.5, 0.9))
omega_grid <- parse_numeric_grid("CLT_N2340_OMEGA_GRID", c(0.01, 0.03, 0.05))
theta_value <- as.numeric(Sys.getenv("CLT_N2340_THETA", as.character(1 / 3)))
gam_value <- as.numeric(Sys.getenv("CLT_N2340_GAM", "0"))

grid <- expand.grid(
  beta = beta_grid,
  omega = omega_grid,
  v = v_grid
)
grid <- grid[order(grid$beta, grid$omega, grid$v), ]

results <- vector("list", nrow(grid))
start_time <- Sys.time()

for (row_id in seq_len(nrow(grid))) {
  beta_value <- grid$beta[row_id]
  omega_value <- grid$omega[row_id]
  v_value <- grid$v[row_id]

  message(
    "Running CLT diagnostics ",
    row_id,
    "/",
    nrow(grid),
    ": n=",
    n_value,
    ", beta=",
    beta_value,
    ", u=",
    u_value,
    ", v=",
    v_value,
    ", omega=",
    omega_value,
    ", reps=",
    n_rep
  )

  clt_summary <- run_clt_replications(
    n_rep = n_rep,
    n = n_value,
    u = u_value,
    v = v_value,
    alpha1 = beta_value,
    alpha2 = beta_value,
    omega1 = omega_value,
    omega2 = omega_value,
    theta = theta_value,
    gam = gam_value,
    noise_model = "iid",
    progress = FALSE
  )

  results[[row_id]] <- cbind(
    data.frame(
      n = n_value,
      beta1 = beta_value,
      beta2 = beta_value,
      u = u_value,
      v = v_value,
      omega = omega_value,
      reps = n_rep
    ),
    clt_summary$summary
  )
}

result_table <- do.call(rbind, results)

output_suffix <- Sys.getenv("CLT_N2340_OUTPUT_SUFFIX", "")

csv_file <- file.path(out_dir, paste0("reviewer6_clt_n2340_uv_beta_omega", output_suffix, ".csv"))
write.csv(result_table, file = csv_file, row.names = FALSE)

latex_file <- file.path(tables_dir, paste0("reviewer6_clt_n2340_uv_beta_omega", output_suffix, "_table.tex"))
lines <- c(
  "\\begin{table}[htbp!]",
  "\\centering",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{3.5pt}",
  "\\caption{Feasible CLT diagnostics for $n=2340$ under different jump activities, tuning values, and noise levels.}",
  "\\label{tab:reviewer6_clt_n2340_uv_beta_omega}",
  "\\begin{tabular}{cccccccccc}",
  "\\toprule",
  "$\\beta$ & $\\omega$ & $v$ & relbias.mean & rel.sd & rel.mse & stud.mean & stud.sd & coverage95 & reject95 \\\\",
  "\\midrule"
)

current_beta <- NA
for (i in seq_len(nrow(result_table))) {
  if (is.na(current_beta) || result_table$beta1[i] != current_beta) {
    if (!is.na(current_beta)) {
      lines <- c(lines, "\\addlinespace")
    }
    lines <- c(
      lines,
      sprintf(
        "\\multicolumn{10}{l}{Panel: $\\beta_1=\\beta_2=%.1f$} \\\\",
        result_table$beta1[i]
      )
    )
    current_beta <- result_table$beta1[i]
  }

  lines <- c(
    lines,
    paste(
      c(
        sprintf("%.1f", result_table$beta1[i]),
        sprintf("%.2f", result_table$omega[i]),
        sprintf("%.1f", result_table$v[i]),
        sprintf("%.4f", result_table$relbias_mean[i]),
        sprintf("%.4f", result_table$relbias_sd[i]),
        sprintf("%.4f", result_table$relbias_mse[i]),
        sprintf("%.4f", result_table$studentized_mean[i]),
        sprintf("%.4f", result_table$studentized_sd[i]),
        sprintf("%.4f", result_table$coverage95[i]),
        sprintf("%.4f", result_table$reject95[i])
      ),
      collapse = " & "
    )
  )
  lines[length(lines)] <- paste0(lines[length(lines)], " \\\\")
}

lines <- c(
  lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\vspace{0.1cm}",
  "\\begin{minipage}{0.98\\textwidth}",
  "\\scriptsize",
  "Note: The sample size and first tuning parameter are fixed at $n=2340$ and $u=1$, respectively. Each design uses 300 Monte Carlo replications under iid microstructure noise with standard deviation $\\omega$.",
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(lines, latex_file)

message("CLT diagnostics completed.")
message("Elapsed time: ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2), " minutes")
message("Saved CSV: ", csv_file)
message("Saved LaTeX table: ", latex_file)
