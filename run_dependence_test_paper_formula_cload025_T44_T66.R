get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
wrapper_path <- file.path(script_dir, "dependence_test_sim_uvn_grid_paper_formula_cload025.R")
out_dir <- file.path(script_dir, "dependence_test_outputs_uvn_grid")
latex_dir <- file.path(out_dir, "latex_tables")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
if (!dir.exists(latex_dir)) {
  dir.create(latex_dir, recursive = TRUE)
}

format_omega_tag <- function(omega_value) {
  gsub("\\.", "", as.character(omega_value))
}

make_tag <- function(T_days,
                     reps = "150",
                     omega = "0.01") {
  paste0(
    "paperformula_cload025_T", T_days,
    "_uv1_v081013_n2340_4680_23400_rhoall_omega",
    format_omega_tag(omega),
    "_reps", reps
  )
}

run_one_T <- function(T_days) {
  reps_value <- Sys.getenv("DEPENDENCE_TEST_UVN_REPS", "150")
  cores_value <- Sys.getenv("DEPENDENCE_TEST_UVN_CORES", "8")
  omega_value <- Sys.getenv("DEPENDENCE_TEST_OMEGA", "0.01")
  output_tag <- make_tag(T_days, reps = reps_value, omega = omega_value)

  env_values <- c(
    paste0("DEPENDENCE_TEST_UVN_REPS=", reps_value),
    paste0("DEPENDENCE_TEST_UVN_CORES=", cores_value),
    paste0("DEPENDENCE_TEST_OMEGA=", omega_value),
    paste0("DEPENDENCE_TEST_UVN_T=", T_days),
    "DEPENDENCE_TEST_UVN_N_GRID=2340,4680,23400",
    "DEPENDENCE_TEST_UVN_RHO_GRID=0,0.2,0.5,-0.5,0.8",
    "DEPENDENCE_TEST_UVN_UV_GRID=1:0.8,1:1,1:1.3",
    paste0("DEPENDENCE_TEST_UVN_TAG=", output_tag)
  )

  env_names <- sub("=.*$", "", env_values)
  env_payload <- sub("^[^=]*=", "", env_values)
  old_env <- Sys.getenv(env_names, unset = NA_character_)
  names(old_env) <- env_names
  on.exit({
    for (env_name in names(old_env)) {
      if (is.na(old_env[[env_name]])) {
        Sys.unsetenv(env_name)
      } else {
        do.call(Sys.setenv, as.list(setNames(old_env[[env_name]], env_name)))
      }
    }
  }, add = TRUE)

  do.call(Sys.setenv, as.list(setNames(env_payload, env_names)))

  message("Running c_load=0.25 dependence-test simulation for T=", T_days)
  source(wrapper_path, local = new.env(parent = globalenv()))

  file.path(out_dir, paste0("dependence_test_uvn_grid_", output_tag, ".csv"))
}

render_combined_reject_table <- function(results, file_path) {
  fmt <- function(x) sprintf("%.2f", x)

  T_values <- c(44, 66)
  n_values <- c(2340, 4680, 23400)
  uv_grid <- data.frame(
    u = c(1, 1, 1),
    v = c(0.8, 1.0, 1.3)
  )
  rho_values <- c(0, 0.2, 0.5, -0.5, 0.8)

  lines <- c(
    "\\begin{table}[htbp!]",
    "\\centering",
    "{\\color{blue}",
    "\\caption{Finite-sample rejection frequencies at the nominal $5\\%$ level for the pointwise dependence test.}",
    "\\label{tab:dep-test-cload025}",
    "\\setlength{\\tabcolsep}{1.2mm}",
    "\\renewcommand{\\arraystretch}{0.95}",
    "\\resizebox{0.92\\textwidth}{!}{",
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}llcccccc}",
    "\\toprule",
    "&&& $T=44$ && &$T=66$& \\\\",
    "\\cmidrule{3-5}\\cmidrule{6-8}",
    "$u$ & $v$ & $n=2340$ & $n=4680$ & $n=23400$  &$n=2340$ & $n=4680$ & $n=23400$ \\\\",
    "\\midrule"
  )

  for (rho_value in rho_values) {
    lines <- c(lines, paste0("\\multicolumn{8}{l}{$\\rho_\\tau=", rho_value, "$} \\\\"))
    for (i in seq_len(nrow(uv_grid))) {
      row_values <- character(length(T_values) * length(n_values))
      value_id <- 1L

      for (T_value in T_values) {
        for (n_value in n_values) {
          match_row <- results[
            results$T_days == T_value &
              results$n_day == n_value &
              results$rho_dep == rho_value &
              results$u == uv_grid$u[i] &
              results$v == uv_grid$v[i],
          ]

          if (nrow(match_row) != 1L) {
            stop(
              "Cannot find exactly one row for T=", T_value,
              ", n=", n_value,
              ", rho=", rho_value,
              ", u=", uv_grid$u[i],
              ", v=", uv_grid$v[i],
              call. = FALSE
            )
          }

          row_values[value_id] <- fmt(match_row$reject05)
          value_id <- value_id + 1L
        }
      }

      lines <- c(
        lines,
        paste0(
          uv_grid$u[i], " & ", uv_grid$v[i], " & ",
          paste(row_values, collapse = " & "),
          " \\\\"
        )
      )
    }

    if (!identical(rho_value, tail(rho_values, 1))) {
      lines <- c(lines, "\\midrule")
    }
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular*}}",
    "}",
    "\\end{table}"
  )

  writeLines(lines, con = file_path)
}

T_grid <- c(44, 66)
csv_paths <- vapply(T_grid, run_one_T, character(1))
result_list <- lapply(csv_paths, read.csv)
combined_results <- do.call(rbind, result_list)
combined_results <- combined_results[
  order(
    combined_results$rho_dep,
    combined_results$u,
    combined_results$v,
    combined_results$T_days,
    combined_results$n_day
  ),
]

reps_value <- Sys.getenv("DEPENDENCE_TEST_UVN_REPS", "150")
omega_value <- Sys.getenv("DEPENDENCE_TEST_OMEGA", "0.01")
combined_tag <- paste0(
  "paperformula_cload025_T44_T66_uv1_v081013_n2340_4680_23400_rhoall_omega",
  format_omega_tag(omega_value),
  "_reps", reps_value
)

combined_csv_path <- file.path(
  out_dir,
  paste0("dependence_test_uvn_grid_", combined_tag, "_combined.csv")
)
combined_latex_path <- file.path(
  latex_dir,
  paste0("dependence_test_uvn_grid_", combined_tag, "_reject05_table.tex")
)

write.csv(combined_results, file = combined_csv_path, row.names = FALSE)
render_combined_reject_table(combined_results, file_path = combined_latex_path)

message("Combined results saved: ", combined_csv_path)
message("Combined LaTeX table saved: ", combined_latex_path)
