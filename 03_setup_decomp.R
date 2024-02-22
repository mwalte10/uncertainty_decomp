#!/usr/bin/env Rscript
rm(list = ls())
library(eppasm)
devtools::load_all("C:/Users/mwalters/frogger/")
devtools::load_all("C:/Users/mwalters/leapfrog/")
library(data.table)
library(dplyr)
library(tidyr)
library(parallel)
library(snow)
source("C:/Users/mwalters/Documents/Projects/paediatric_uncertainty/workflow/central_functions/recreate_uncertainty_functions.R")
source("C:/Users/mwalters/ihme-imperial/gets_weighted_pairs.R")

args <- commandArgs(trailingOnly = TRUE)
run_name <- args[1]
print(run_name)

output_dir <- paste0("C:/Users/mwalters/Documents/Projects/paediatric_uncertainty/uncertainty_calc/runs", run_name)
base_dir <- "C:/Users/mwalters/Documents/Projects/paediatric_uncertainty/workflow/"
input_dir <- paste0(base_dir, "/runs/", run_name)
lapply(c(output_dir), dir.create, recursive = T)

pjnz <- paste0(base_dir, "/from_spectrum/Botswana2023v4 WPP 02_03_2023 KOS.PJNZ")
demp <- readRDS(paste0(base_dir, "demp.RDS"))
hivp <- readRDS(paste0(base_dir, "hivp.RDS"))
inc_rate <- readRDS(paste0(base_dir, "inc_rate_draws.RDS"))
mat_prev <- readRDS(paste0(base_dir, "hiv_births_draws.RDS"))
mat_prev_sd <- apply(mat_prev[, -(1:2)], MARGIN = 2, FUN = sd)
out_prop_gte350 <- readRDS(paste0(base_dir, "prop_gte350_draws.RDS"))
out_prop_gte350[is.nan(prop_gte350), prop_gte350 := 0]
out_prop_lt200 <- readRDS(paste0(base_dir, "prop_lt200_draws.RDS"))
out_prop_lt200[is.nan(prop_lt200), prop_lt200 := 0]

stored_hivp_mtct <- hivp$mtct
stored_hivp_pmtct_mtct <- hivp$pmtct_mtct
stored_hivp_ctx_effect <- hivp$ctx_effect
stored_hivp_paed_art_val <- hivp$paed_art_val
saved_hivp <- hivp

uncertainty_params <- readRDS(paste0(output_dir, "uncertainty_params.RDS"))
uncertainty_params <- unique(c(setdiff(uncertainty_params, 'paed_art_num'), 'paed_art_val'))
params_list <- readRDS(paste0(output_dir, "/params_list.RDS"))

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

##Note: will need to make this more generalizable when I expand uncertainty 
  run_lmod_new_input <- function(run.name, iter, hivp, demp_x,
                                 input.table, scenario){
    # library(data.table)
    # hivp$tfr <- demp_x$tfr
  
  if(run.name == '/01_recreate_baseline/'){
    mat_hiv_births_flag <- input.table[scenario, ],
    mtct_flag,
    ctx_effect_flag,
    paed_art_val_flag,
    if (mat_hiv_births_flag) {
      hivp$mat_hiv_births <- as.numeric(mat_prev[Iteration ==iter, -(1:2)])
      hivp$prop_gte350 <- out_prop_gte350[iteration ==iter, prop_gte350]
      hivp$prop_lt200 <- out_prop_lt200[iteration ==iter, prop_lt200]
    }

    if (mtct_flag) {
      hivp$mtct <- as.array(mtct_list[[iter]])
      hivp$pmtct_mtct <- as.array(pmtct_mtct_list[[i]])

    }

    if (ctx_effect_flag) {
      hivp$ctx_effect <- as.array(rep(ctx_effect_list[[iter]], 61))
    }

    if (paed_art_val_flag) {
      hivp$paed_art_val <- as.array(paed_art_val_list[[iter]])
    }
  }

  if(run.name == '/02_expand_baseline/'){
    if (mat_hiv_births_flag) {
      hivp$mat_hiv_births <- as.array(mat_hiv_births_list[[iter]])
      prop_gte350_mean <- unique(out_prop_gte350[,prop_gte350 := mean(prop_gte350), by = 'year']$prop_gte350)
      prop_lt200_mean <- unique(out_prop_lt200[,prop_lt200 := mean(prop_lt200), by = 'year']$prop_lt200)
      hivp$prop_gte350 <-  prop_gte350_mean
      hivp$prop_lt200 <-  prop_lt200_mean
    }

    if (mtct_flag) {
      hivp$mtct <- as.array(mtct_list[[iter]])
      hivp$pmtct_mtct <- as.array(pmtct_mtct_list[[iter]])

    }

    if (ctx_effect_flag) {
      hivp$ctx_effect <- as.array(rep(ctx_effect_list[[iter]], 61))
    }

    if (paed_art_val_flag) {
      hivp$paed_art_val <- as.array(paed_art_val_list[[iter]])
    }

  }

  if(run.name == '/03_new_sd_mtct/'){
    if (mat_hiv_births_flag) {
      hivp$mat_hiv_births <- as.numeric(mat_prev[Iteration == iter, -(1:2)])
      hivp$prop_gte350 <- out_prop_gte350[iteration == iter, prop_gte350]
      hivp$prop_lt200 <- out_prop_lt200[iteration == iter, prop_lt200]
    }

    if (mtct_flag) {
      hivp$mtct <- as.array(mtct_list[[iter]])
      hivp$pmtct_mtct <- as.array(pmtct_mtct_list[[iter]])

    }

    if (ctx_effect_flag) {
      hivp$ctx_effect <- as.array(rep(ctx_effect_list[[iter]], 61))
    }

    if (paed_art_val_flag) {
      hivp$paed_art_val <- as.array(paed_art_val_list[[iter]])
    }
  }
  
  lmod <- leapfrog::leapfrogR(demp_x, hivp)
  

  return(lmod)

}

## decomp pipeline is located here: C:\Users\mwalters\ihme-imperial\pipeline
input <- lapply(1, get_weighted_pairs, var_vec = setdiff(uncertainty_params, "pmtct_mtct")) ## removing pmtct_mtct as that will be varied alongside mtct
input <- input[[1]]$total_combos

unparms <- paste0(input_dir, '/uncertainty_files/', uncertainty_params, '/') 
if(run_name == '/03_new_sd_mtct/'){
  unparms <- unparms[1:4]
}

for(dir in  unparms){
  assign(paste0(uncertainty_params[which(dir == unparms)], '_list'), lapply(paste0(dir,1:300, '.qs'), qs::qread))
  print(dir)
}


start <- Sys.time()
for (scenario in 1:16) {
  if (!dir.exists(paste0(output_dir, "/scenario_", scenario))) {
    dir.create(paste0(output_dir, "/scenario_", scenario), recursive = T)
  }

  under5incidence <- array(dim = c(5, 2, 61, 300))
  tot_pop <- array(dim = c(15, 2, 61, 300))
  artpopu15 <- array(dim = c(3, 7, 15, 2, 61, 300))
  hivpopu15 <- array(dim = c(7, 4, 15, 2, 61, 300))
  deaths_artpop <- array(dim = c(15, 2, 61, 300))
  deaths_noartpop <- array(dim = c(15, 2, 61, 300))
  
  for(i in 1:300){
    lmod <- run_lmod_new_input(iter = i, run.name = run_name, 
                       mat_hiv_births_flag = input[scenario,mat_hiv_births],
                       mtct_flag = input[scenario,mtct],
                       ctx_effect_flag = input[scenario,ctx_effect],
                       paed_art_val_flag = input[scenario,paed_art_val],
                       hivp = saved_hivp, 
                       demp_x = demp)
    
    deaths_artpop[, , , i] <- apply(lmod$aidsdeaths_art_paed, MARGIN = c(3,4,5), FUN = sum) 
    deaths_noartpop[, , , i] <- apply(lmod$aidsdeaths_noart_paed, MARGIN = c(3,4,5), FUN = sum)
    under5incidence[, , , i] <- lmod$infections[1:5, , ]
    artpopu15[, , , , , i] <- lmod$artstrat_paeds
    hivpopu15[, , , , , i] <- lmod$hivstrat_paeds
    tot_pop[, , , i] <- lmod$totpop1[1:15, , ]
    print((i + ((scenario-1) * 300)) / (300 * 16))
  }
  
  deaths <- deaths_artpop + deaths_noartpop
  deaths_x <- data.table(melt(deaths))
  deaths_x <- deaths_x[, .(age = Var1 - 1, sex = ifelse(Var2 == 1, "male", "female"), year = Var3 + 1969, draw = Var4, value)]
  deaths_x <- deaths_x[, .(value = sum(value)), by = c("age", "year", "draw")]
  # ggplot(deaths_x[age == 0], aes(year, value, group = as.factor(draw))) + geom_line() + facet_wrap(~age, scales = 'free')

  hivpop <- data.table(melt(hivpopu15))
  hivpop <- hivpop[, .(age = Var3 - 1, sex = ifelse(Var4 == 1, "male", "female"), year = Var5 + 1969, draw = Var6, value)]
  hivpop <- hivpop[, .(value = sum(value)), by = c("age", "year", "draw")]
  # ggplot(hivpop[age == 0], aes(year, value, group = as.factor(draw))) + geom_line() + facet_wrap(~age, scales = 'free')
  
  totpop <- data.table(melt(tot_pop))
  totpop <- totpop[, .(age = Var1 - 1, sex = ifelse(Var4 == 2, "male", "female"), year = Var3 + 1969, draw = Var4, value)]
  totpop <- totpop[, .(value = sum(value)), by = c("age", "year", "draw")]

  hivnpop <- merge(hivpop, totpop, by = c("age", "year", "draw"))
  hivnpop <- hivnpop[, .(age, year, draw, value = value.y - value.x)]

  ## maybe need to collapse to under 3?
  inf <- data.table(melt(under5incidence))
  inf <- inf[, .(age = Var1 - 1, sex = ifelse(Var2 == 1, "male", "female"), year = Var3 + 1969, draw = Var4, value)]
  inf <- inf[, .(value = sum(value)), by = c("age", "year", "draw")]
  inf_rate <- merge(inf, hivnpop, by = c("year", "draw", "age"))
  inf_rate <- inf_rate[, .(value = value.x / value.y, year, draw, age)]
  
  ##Switching to PAF of deaths due to HIV among the whole population to remove sensitivity to small numbers of HIV positive pop
  death_rate <- merge(deaths_x, totpop, by = c("age", "year", 'draw'), allow.cartesian = T)
  death_rate <- death_rate[, .(value = value.x / value.y, age, year, draw)]

  out_sd <- list(
    deaths = death_rate,
    deaths_art = deaths_artpop,
    deaths_noart = deaths_noartpop,
    inf = inf_rate,
    hivpop = hivpop,
    hivnpop = hivnpop
  )
  saveRDS(out_sd, file = paste0(output_dir, "/scenario_", scenario, "/results.RDS"))
  # print(scenario)
}
end <- Sys.time()
diff <- start - end
diff
