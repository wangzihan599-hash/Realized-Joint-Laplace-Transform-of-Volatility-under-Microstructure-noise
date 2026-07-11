# Reviewer 2, Comment 8:
# Report finite-sample coverage probabilities and rejection frequencies for
# nominal 95% feasible confidence intervals.

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

set.seed(as.integer(Sys.getenv("REVIEWER_SEED", "20260426")))

out_dir <- file.path(script_dir, "reviewer_outputs")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

n_rep <- as.integer(Sys.getenv("REVIEWER8_REPS", "300"))

omega_grid <- parse_numeric_grid("REVIEWER8_OMEGA_GRID", c(0.01, 0.03, 0.05))
v_grid <- parse_numeric_grid("REVIEWER8_V_GRID", c(0.5, 1.0, 1.5))

n_value <- as.integer(Sys.getenv("REVIEWER8_N", "23400"))
beta_value <- as.numeric(Sys.getenv("REVIEWER8_BETA", "0.5"))
u_value <- as.numeric(Sys.getenv("REVIEWER8_U", "1"))
theta_value <- as.numeric(Sys.getenv("REVIEWER8_THETA", as.character(1 / 3)))
gam_value <- as.numeric(Sys.getenv("REVIEWER8_GAM", "0"))

summary_rows <- list()
draw_rows <- list()
row_id <- 1L

for (omega_value in omega_grid) {
  for (v_value in v_grid) {
    clt_run <- run_clt_replications(
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

    summary_rows[[row_id]] <- cbind(
      data.frame(
        omega = omega_value,
        u = u_value,
        v = v_value,
        n = n_value,
        beta = beta_value,
        reps = n_rep
      ),
      clt_run$summary
    )

    draw_rows[[row_id]] <- cbind(
      data.frame(
        omega = omega_value,
        u = u_value,
        v = v_value
      ),
      clt_run$draws
    )

    row_id <- row_id + 1L
  }
}

summary_table <- do.call(rbind, summary_rows)
draws_table <- do.call(rbind, draw_rows)

write.csv(
  summary_table,
  file = file.path(out_dir, "reviewer8_coverage_summary.csv"),
  row.names = FALSE
)
write.csv(
  draws_table,
  file = file.path(out_dir, "reviewer8_studentized_draws.csv"),
  row.names = FALSE
)

message("Reviewer 8 simulations completed.")
message("Saved: ", file.path(out_dir, "reviewer8_coverage_summary.csv"))
message("Saved: ", file.path(out_dir, "reviewer8_studentized_draws.csv"))
