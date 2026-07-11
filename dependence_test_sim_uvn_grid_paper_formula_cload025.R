get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()

# Paper-formula DGP, but with the stronger volatility loading
# sigma_{l,t} = exp(-0.3125 + 0.25 tau_{l,t}), l = 1,2.
#
# The remaining DGP is the same as the manuscript simulation:
# tau is generated day by day from the AR(1) specification, and the two
# price processes share the common Brownian shock W_t.
Sys.setenv(DEPENDENCE_TEST_TAU_MODE = "daily")
Sys.setenv(DEPENDENCE_TEST_PRICE_SHOCK_MODE = "common")
Sys.setenv(DEPENDENCE_TEST_CLOAD = "0.25")

set_default_env <- function(name, value) {
  if (Sys.getenv(name, unset = "") == "") {
    Sys.setenv(structure(value, names = name))
  }
}

set_default_env("DEPENDENCE_TEST_UVN_REPS", "150")
set_default_env("DEPENDENCE_TEST_UVN_CORES", "8")
set_default_env("DEPENDENCE_TEST_UVN_T", "44")
set_default_env("DEPENDENCE_TEST_UVN_N_GRID", "2340,4680,23400")
set_default_env("DEPENDENCE_TEST_UVN_RHO_GRID", "0,0.2,0.5,-0.5,0.8")
set_default_env("DEPENDENCE_TEST_UVN_UV_GRID", "1:0.8,1:1,1:1.3")
set_default_env("DEPENDENCE_TEST_OMEGA", "0.01")

if (Sys.getenv("DEPENDENCE_TEST_UVN_TAG", unset = "") == "") {
  omega_tag <- gsub("\\.", "", Sys.getenv("DEPENDENCE_TEST_OMEGA", "0.01"))
  default_tag <- paste0(
    "paperformula_cload025_T", Sys.getenv("DEPENDENCE_TEST_UVN_T", "44"),
    "_uv1_v081013_n2340_4680_23400_rhoall_omega", omega_tag,
    "_reps", Sys.getenv("DEPENDENCE_TEST_UVN_REPS", "150")
  )
  Sys.setenv(DEPENDENCE_TEST_UVN_TAG = default_tag)
}

source(file.path(script_dir, "dependence_test_sim_uvn_grid.R"))
