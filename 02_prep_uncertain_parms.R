#!/usr/bin/env Rscript
# args <- commandArgs(trailingOnly = TRUE)
# run_name <- args[1]
# uncertainty_params <- args[2]
# uncertainty_params <- unlist(strsplit(uncertainty_params, split = " "))

rnorm_alt <- function(mean_x, sd.mult = 0.2, n = 1) {
  rnorm(n, mean = mean_x, sd = sd.mult * mean_x)
}

vary_param <- function(n, param, sd.mult = 0.2, draw = 1) {
  if (length(dim(param)) == 0) {
    if(length(param) > 1){
      unlist(lapply(param, rnorm_alt, sd.mult = 0.2, n = 1))
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

get_quantile <- function(mean_x, sd_scalar = 0.2, quant){
  out <- quantile(rnorm(1000, mean = mean_x, sd = mean_x * sd_scalar), quant)
  out <- unname(out)
  return(out)
}

get_quantile_noscalar <- function(mean_x, sd_x = 0.2, quant){
  out <- quantile(rnorm(1000, mean = mean_x, sd = sd_x), quant)
  out <- unname(out)
  return(out)
}

prep_mtct_draws <- function(){
  root <- getwd()

  new_mtct <- fread(paste0(root, '/new_trans_rates_posterior.csv'))
  mtct <- array(data = NA, dim = c(7,2), dimnames = list(cd4 = unique(row.names(hivp$mtct)), bf = unique(colnames(hivp$mtct))))
  mtct_sd <- array(data = NA, dim = c(7,2), dimnames = list(cd4 = unique(row.names(hivp$mtct)), bf = unique(colnames(hivp$mtct))))
  mtct[1:3,1] <- new_mtct[cd4_range == '>350' & trans_timing == 'Perinatal', mean]
  mtct_sd[1:3,1] <- new_mtct[cd4_range == '>350' & trans_timing == 'Perinatal', sd]
  mtct[1:3,2] <- new_mtct[cd4_range == '>350' & trans_timing == 'Breastfeeding', mean]
  mtct_sd[1:3,2] <- new_mtct[cd4_range == '>350' & trans_timing == 'Breastfeeding', sd]

  mtct[4:5,1] <- new_mtct[cd4_range == '[200-350)' & trans_timing == 'Perinatal', mean]
  mtct_sd[4:5,1] <- new_mtct[cd4_range == '[200-350)' & trans_timing == 'Perinatal', sd]
  mtct[4:5,2] <- new_mtct[cd4_range == '[200-350)' & trans_timing == 'Breastfeeding', mean]
  mtct_sd[4:5,2] <- new_mtct[cd4_range == '[200-350)' & trans_timing == 'Breastfeeding', sd]

  mtct[6:7,1] <- new_mtct[cd4_range == '[0-200)' & trans_timing == 'Perinatal', mean]
  mtct_sd[6:7,1] <- new_mtct[cd4_range == '[0-200)' & trans_timing == 'Perinatal', sd]
  mtct[6:7,2] <- new_mtct[cd4_range == '[0-200)' & trans_timing == 'Breastfeeding', mean]
  mtct_sd[6:7,2] <- new_mtct[cd4_range == '[0-200)' & trans_timing == 'Breastfeeding', sd]

  pmtct_mtct <- array(data = NA, dim = c(7,7,2), dimnames = list(cd4 = unique(row.names(hivp$mtct)), pmtct_reg = colnames(hivp$pmtct_mtct[,,1]), bf = unique(colnames(hivp$mtct))))
  pmtct_mtct_sd <- array(data = NA, dim = c(7,7,2), dimnames = list(cd4 = unique(row.names(hivp$mtct)), pmtct_reg = colnames(hivp$pmtct_mtct[,,1]),  bf = unique(colnames(hivp$mtct))))
  pmtct_mtct[5:7,1,1] <- new_mtct[cd4_range == 'Option A' & trans_timing == 'Perinatal', mean]
  pmtct_mtct_sd[5:7,1,1] <- new_mtct[cd4_range == 'Option A' & trans_timing == 'Perinatal', sd]
  pmtct_mtct[5:7,1,2] <- new_mtct[cd4_range == 'Option A' & trans_timing == 'Breastfeeding', mean]
  pmtct_mtct_sd[5:7,1,2] <- new_mtct[cd4_range == 'Option A' & trans_timing == 'Breastfeeding', sd]

  pmtct_mtct[5:7,2,1] <- new_mtct[cd4_range == 'Option B' & trans_timing == 'Perinatal', mean]
  pmtct_mtct_sd[5:7,2,1] <- new_mtct[cd4_range == 'Option B' & trans_timing == 'Perinatal', sd]
  pmtct_mtct[5:7,2,2] <- new_mtct[cd4_range == 'Option B' & trans_timing == 'Breastfeeding', mean]
  pmtct_mtct_sd[5:7,2,2] <- new_mtct[cd4_range == 'Option B' & trans_timing == 'Breastfeeding', sd]

  pmtct_mtct[1:7,3,1] <- new_mtct[cd4_range == 'SDNVP' & trans_timing == 'Perinatal', mean]
  pmtct_mtct_sd[1:7,3,1] <- new_mtct[cd4_range == 'SDNVP' & trans_timing == 'Perinatal', sd]
  pmtct_mtct[1:7,3,2] <- new_mtct[cd4_range == 'SDNVP' & trans_timing == 'Breastfeeding', mean]
  pmtct_mtct_sd[1:7,3,2] <- new_mtct[cd4_range == 'SDNVP' & trans_timing == 'Breastfeeding', sd]

  pmtct_mtct[,4,1] <- new_mtct[cd4_range == 'Dual ARV' & trans_timing == 'Perinatal', mean]
  pmtct_mtct_sd[,4,1] <- new_mtct[cd4_range == 'Dual ARV' & trans_timing == 'Perinatal', sd]
  pmtct_mtct[,4,2] <- new_mtct[cd4_range == 'Dual ARV' & trans_timing == 'Breastfeeding', mean]
  pmtct_mtct_sd[,4,2] <- new_mtct[cd4_range == 'Dual ARV' & trans_timing == 'Breastfeeding', sd]

  pmtct_mtct[1:7,5,1] <- new_mtct[cd4_range == 'Option B+, on ART preconception' & trans_timing == 'Perinatal', mean]
  pmtct_mtct_sd[1:7,5,1] <- new_mtct[cd4_range == 'Option B+, on ART preconception' & trans_timing == 'Perinatal', sd]
  pmtct_mtct[1:7,5,2] <- new_mtct[cd4_range == 'Option B+, on ART preconception' & trans_timing == 'Breastfeeding', mean]
  pmtct_mtct_sd[1:7,5,2] <- new_mtct[cd4_range == 'Option B+, on ART preconception' & trans_timing == 'Breastfeeding', sd]

  pmtct_mtct[1:7,6,1] <- new_mtct[cd4_range == 'Option B+, on ART 5-39 weeks' & trans_timing == 'Perinatal', mean]
  pmtct_mtct_sd[1:7,6,1] <- new_mtct[cd4_range == 'Option B+, on ART 5-39 weeks' & trans_timing == 'Perinatal', sd]
  pmtct_mtct[1:7,6,2] <- new_mtct[cd4_range == 'Option B+, on ART 5-39 weeks' & trans_timing == 'Breastfeeding', mean]
  pmtct_mtct_sd[1:7,6,2] <- new_mtct[cd4_range == 'Option B+, on ART 5-39 weeks' & trans_timing == 'Breastfeeding', sd]

  pmtct_mtct[1:7,7,1] <- new_mtct[cd4_range == 'Option B+, on ART <4 weeks' & trans_timing == 'Perinatal', mean]
  pmtct_mtct_sd[1:7,7,1] <- new_mtct[cd4_range == 'Option B+, on ART <4 weeks' & trans_timing == 'Perinatal', sd]
  pmtct_mtct[1:7,7,2] <- new_mtct[cd4_range == 'Option B+, on ART <4 weeks' & trans_timing == 'Breastfeeding', mean]
  pmtct_mtct_sd[1:7,7,2] <- new_mtct[cd4_range == 'Option B+, on ART <4 weeks' & trans_timing == 'Breastfeeding', sd]
  pmtct_mtct[is.na(pmtct_mtct)] <- 0
  pmtct_mtct_sd[is.na(pmtct_mtct_sd)] <- 0
}

prep_input_draws <- function(run_name, uncertainty_params){
  root <- getwd()

  demp <- readRDS(paste0(root, "/demp.RDS"))
  hivp <- readRDS(paste0(root, "/hivp.RDS"))
  inc_rate <- readRDS(paste0(root, "/inc_rate_draws.RDS"))
  mat_prev <- readRDS(paste0(root, "/hiv_births_draws.RDS"))
  mat_prev_sd <- apply(mat_prev[, -(1:2)], MARGIN = 2, FUN = sd)
  out_prop_gte350 <- readRDS(paste0(root, "/prop_gte350_draws.RDS"))
  out_prop_lt200 <- readRDS(paste0(root, "/prop_lt200_draws.RDS"))

  output_dir <- paste0(root, "/runs/", run_name, "/")
  lapply(c(output_dir, paste0(output_dir, "/hivp/"), paste0(output_dir, "/results/")), dir.create, recursive = T)

  input_dir <- paste0(root, "/runs/", run_name, '/uncertainty_files/')
  lapply(paste0(input_dir, '/', uncertainty_params, '/'),
         dir.create, recursive = T)

  stored_hivp_mtct <- hivp$mtct
  stored_hivp_pmtct_mtct <- hivp$pmtct_mtct
  stored_hivp_ctx_effect <- hivp$ctx_effect
  stored_hivp_paed_art_val <- hivp$paed_art_val
  saved_hivp <- hivp

 lapply(1:300, function(i) {
  if(run_name == '/01_recreate_baseline/'){
    mtct <- as.array(vary_param(param = stored_hivp_mtct, sd.mult = 0.05, n = 1))
    pmtct_mtct <- as.array(vary_param(param = stored_hivp_pmtct_mtct, sd.mult = 0.05, n = 1))
    mtct[mtct < 0] <- 0
    pmtct_mtct[pmtct_mtct < 0] <- 0
    qs::qsave(mtct, file = paste0(input_dir, '/mtct/', i, '.qs'))
    qs::qsave(pmtct_mtct, file = paste0(input_dir, '/pmtct_mtct/', i, '.qs'))

    ctx_effect <- array(rnorm(n = 1, mean = stored_hivp_ctx_effect, sd = 0.1))
    qs::qsave(ctx_effect, file = paste0(input_dir, '/ctx_effect/', i, '.qs'))

    scalar = rnorm(1, mean = 1, sd = 0.05)
    if(scalar > 2){
      scalar = 2
    }
    paed_art_val <- array(stored_hivp_paed_art_val * scalar)
    qs::qsave(paed_art_val, file = paste0(input_dir, '/paed_art_val/', i, '.qs'))

  }

  if(run_name == '/02_expand_baseline/'){

    quant <- runif(n = 1, min = 0 , max = 1)
    mat_hiv_births <- unlist(lapply(hivp$mat_hiv_births, get_quantile, sd_scalar = 0.15, quant))
    qs::qsave(mat_hiv_births, file = paste0(input_dir, '/mat_hiv_births/', i, '.qs'))

    mtct <- as.array(vary_param(param = stored_hivp_mtct, sd.mult = 0.15, n = 1))
    pmtct_mtct <- as.array(vary_param(param = stored_hivp_pmtct_mtct, sd.mult = 0.15, n = 1))
    mtct[mtct < 0] <- 0
    pmtct_mtct[pmtct_mtct < 0] <- 0
    qs::qsave(mtct, file = paste0(input_dir, '/mtct/', i, '.qs'))
    qs::qsave(pmtct_mtct, file = paste0(input_dir, '/pmtct_mtct/', i, '.qs'))

    scalar = rnorm(1, mean = 1, sd = 0.15)
    ctx_effect <- array(stored_hivp_ctx_effect * scalar)
    qs::qsave(ctx_effect, file = paste0(input_dir, '/ctx_effect/', i, '.qs'))


    scalar = rnorm(1, mean = 1, sd = 0.15)
    if(scalar > 2){
      scalar = 2
    }
    paed_art_val <- array(stored_hivp_paed_art_val * scalar)
    qs::qsave(paed_art_val, file = paste0(input_dir, '/paed_art_val/', i, '.qs'))


  }

  if(run_name == '/03_new_sd_mtct/'){
    quant <- runif(n = 1, min = 0 , max = 1)
    mtct <- array(unlist(lapply(1:length(mtct), vary_param_normal_sd, param_mean = mtct, param_sd = mtct_sd, quant)),
                  dim = c(7,2),
                  dimnames = list(cd4 = unique(row.names(hivp$mtct)), bf = unique(colnames(hivp$mtct))))
    pmtct_mtct <- array(unlist(lapply(1:length(pmtct_mtct), vary_param_normal_sd, param_mean = pmtct_mtct, param_sd = pmtct_mtct_sd, quant)),
                        dim = c(7,7,2),
                        dimnames = list(cd4 = unique(row.names(hivp$pmtct_mtct)), pmtct_reg = colnames(hivp$pmtct_mtct[,,1]), bf = unique(colnames(hivp$mtct))))

    mtct[mtct < 0] <- 0
    pmtct_mtct[pmtct_mtct < 0] <- 0
    qs::qsave(mtct, file = paste0(input_dir, '/mtct/', i, '.qs'))
    qs::qsave(pmtct_mtct, file = paste0(input_dir, '/pmtct_mtct/', i, '.qs'))


    ctx_effect <- array(rnorm(n = 1, mean = stored_hivp_ctx_effect, sd = 0.1))
    qs::qsave(ctx_effect, file = paste0(input_dir, '/ctx_effect/', i, '.qs'))


    scalar = rnorm(1, mean = 1, sd = 0.05)
    if(scalar > 2){
      scalar = 2
    }
    paed_art_val <- array(stored_hivp_paed_art_val * scalar)
    qs::qsave(paed_art_val, file = paste0(input_dir, '/paed_art_val/', i, '.qs'))
  }

  if(run_name == '/04_all_hivp_sources/'){

    #########################################################################
    ###Incidence params
    #########################################################################
    #mtct_rates
    mtct <- as.array(vary_param(param = stored_hivp_mtct, sd.mult = 0.15, n = 1))
    pmtct_mtct <- as.array(vary_param(param = stored_hivp_pmtct_mtct, sd.mult = 0.15, n = 1))
    mtct[mtct < 0] <- 0
    pmtct_mtct[pmtct_mtct < 0] <- 0
    qs::qsave(list(mtct = mtct, pmtct_mtct = pmtct_mtct), file = paste0(input_dir, '/mtct_rates/', i, '.qs'))

    #mat_hiv_births
    quant <- runif(n = 1, min = 0 , max = 1)
    mat_hiv_births <- unlist(lapply(hivp$mat_hiv_births, get_quantile, sd_scalar = 0.15, quant))
    qs::qsave(mat_hiv_births, file = paste0(input_dir, '/mat_hiv_births/', i, '.qs'))

    #pmtct
    scalar = rnorm(1, mean = 1, sd = 0.15)
    if(scalar > 1){
      scalar = 1
    }
    pmtct <- array(hivp$pmtct * scalar)
    qs::qsave(pmtct, file = paste0(input_dir, '/pmtct/', i, '.qs'))


    #bf_duration
    bf_duration <- as.array(vary_param(param = hivp$bf_duration, sd.mult = 0.15, n = 1))
    bf_duration[bf_duration > 1] <- 1
    qs::qsave(bf_duration, file = paste0(input_dir, '/bf_duration/', i, '.qs'))


    #########################################################################
    ###Treatment params
    #########################################################################
    #paed_art_val
    paed_art_val<- as.array(vary_param(param = stored_hivp_paed_art_val, sd.mult = 0.15, n = 1))
    qs::qsave(paed_art_val, file = paste0(input_dir, '/paed_art_val/', i, '.qs'))

    #ctx_effect
    scalar = rnorm(1, mean = 1, sd = 0.15)
    ctx_effect <- array(stored_hivp_ctx_effect * scalar)
    qs::qsave(ctx_effect, file = paste0(input_dir, '/ctx_effect/', i, '.qs'))

    #ctx_val
    ctx_val<- as.array(vary_param(param = hivp$ctx_val, sd.mult = 0.15, n = 1))
    qs::qsave(ctx_val, file = paste0(input_dir, '/ctx_val/', i, '.qs'))


    #########################################################################
    ###Natural history params
    ########################################################################
    #paed_cd4_prog
    paed_cd4_prog <- as.array(vary_param(param = hivp$paed_cd4_prog, sd.mult = 0.15, n = 1))
    paed_cd4_prog[paed_cd4_prog < 0] <- 0

    #adol_cd4_prog
    adol_cd4_prog <- as.array(vary_param(param = hivp$adol_cd4_prog, sd.mult = 0.15, n = 1))
    adol_cd4_prog[adol_cd4_prog < 0] <- 0

    #paed_art_mort
    paed_art_mort <- as.array(vary_param(param = hivp$paed_art_mort, sd.mult = 0.15, n = 1))
    paed_art_mort[paed_art_mort < 0] <- 0

    #adol_art_mort
    adol_art_mort <- as.array(vary_param(param = hivp$adol_art_mort, sd.mult = 0.15, n = 1))
    adol_art_mort[adol_art_mort < 0] <- 0

    #paed_cd4_mort
    paed_cd4_mort <- as.array(vary_param(param = hivp$paed_cd4_mort, sd.mult = 0.15, n = 1))
    paed_cd4_mort[paed_cd4_mort < 0] <- 0

    #adol_cd4_mort
    adol_cd4_mort <- as.array(vary_param(param = hivp$adol_cd4_mort, sd.mult = 0.15, n = 1))
    adol_cd4_mort[adol_cd4_mort < 0] <- 0

    qs::qsave(list(paed_cd4_prog = paed_cd4_prog, adol_cd4_prog = adol_cd4_prog,
                   paed_cd4_mort = paed_cd4_mort, adol_cd4_mort = adol_cd4_mort), paste0(input_dir, '/natural_history_noart_params/', i, '.qs'))

    qs::qsave(list(paed_art_mort = paed_art_mort, adol_art_mort = adol_art_mort), paste0(input_dir, '/natural_history_art_params/', i, '.qs'))


  }}
)
}



