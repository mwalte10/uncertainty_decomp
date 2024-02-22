## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: Maggie Walters
##
## Email: mwalters@ic.ac.uk
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

rnorm_alt <- function(mean_x, sd.mult = 0.2, n = 1) {
  rnorm(n, mean = mean_x, sd = sd.mult * mean_x)
}


vary_param <- function(n, param, sd.mult = 0.2, draw = 1) {
  if (length(dim(param)) == 0) {
    if(length(param) > 1){
      lapply(param, rnorm_alt, sd.mult = 0.2, n = 1) %>% unlist()
    }else{
      rnorm_alt(mean_x = param[[1]], sd.mult)
    }
  } else {
    apply(param, MARGIN = c(1:length(dim(param))), FUN = rnorm_alt, sd.mult, n)
  }
}

vary_param_normal_sd <- function(index, param_mean, param_sd, quantile) {

  x <- quantile(rnorm(n = 1, mean = param_mean[index], sd = param_sd[index]), quant)
  return(unname(x))
}


get_draw_values <- function(param, draw, sd.mult, run) {
  var <- attr(param, "var")
  print(var)
  attr(param, "var") <- NULL

  if (var %in% c("mtct", "pmtct_mtct", "paed_art_num")) {
    param <- vary_param(param, n = 1, sd.mult)
  } else if (var == "ctx_effect") {
    param <- vary_param(param[[1]], n = 1, sd.mult)
  }

  if (run == "/01_recreate_baseline/" & var == "mat_hiv_births") {
    param <- as.vector(as.numeric(mat_prev[Iteration == draw, -(1:2)]))
  } else if (run == "/02_expand_baseline/" & var == "mat_hiv_births") {
    param <- vary_param(param, n = 1, sd.mult)
  }

  return(param)
}

prepare_hivp_sample <- function(projp, parms, parms_list, sim_number = 300, prop_gte350, prop_lt200, sd.mult = 0.05, run) {
  uncertain_params <- lapply(parms_list, get_draw_values, draw = sim_number, sd.mult, run)

  for (i in 1:length(parms)) {
    projp[[parms[i]]] <- uncertain_params[[i]]
  }

  projp$prop_gte350 <- prop_gte350[iteration == sim_number, prop_gte350]
  projp$prop_lt200 <- prop_lt200[iteration == sim_number, prop_lt200]

  qs::qsave(projp, paste0(output_dir, "/hivp/", sim_number, ".qs"))
}


run_leapfrog <- function(sim_number, outdir, dem) {
  hivp <- qs::qread(paste0(outdir, "/hivp/", sim_number, ".qs"))
  out <- leapfrog::leapfrogR(dem, hivp)
  qs::qsave(out, paste0(outdir, "/results/", sim_number, ".qs"))
}
