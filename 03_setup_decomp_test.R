root <- getwd()
source(paste0(root, '/central_functions/recreate_uncertainty_functions.R'))
source(paste0(root, '/central_functions/gets_weighted_pairs.R'))

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

run_lmod_new_input <- function(iter, scenario, scenario_row, hivp_x, demp_x,
                               input_list,
                                mat.prev, gte350, lt200, output.dir){
  if(length(names(scenario_row)[which(scenario_row == 1)]) > 0){
    special_vec = c("mtct_rates", "natural_history_noart_params", "natural_history_art_params" )
    if(any(special_vec %in% names(scenario_row)[which(scenario_row == 1)])){
      x = lapply(input_list[intersect(special_vec, names(scenario_row)[which(scenario_row == 1)])], '[[', iter)
      x <- unlist(x, recursive=  F)
      names(x) <- unlist(lapply(strsplit(names(x),split = '\\.'), '[[', 2))
    }

    y = lapply(input_list[setdiff(names(scenario_row)[which(scenario_row == 1)], special_vec)], '[[', iter)
    # names(x) <- unlist(lapply(strsplit(names(x),split = '\\.'), '[[', 2))
    x <-  unlist(list(x,y), recursive = F)

    hivp_x[names(x)] <- x
    if('mat_hiv_births' %in% names(x)){
      hivp_x$mat_hiv_births <- as.numeric(mat.prev[iter, -(1:2)])
      hivp_x$prop_gte350 <- as.array(unlist(gte350[iter, -(1)]))
      hivp_x$prop_lt200 <- as.array(unlist(lt200[iter,  -(1)]))
    }
  }

  hivp_x$tfr = demp_x$tfr
  demp_x$projection_period <- 'calendar'
  hivp_x$incrate <- rep(0,61)
  ##save out the smaller thing here
  lmod <-  leapfrog::leapfrogR(demp_x, hivp_x, hiv_strat = "full", hiv_steps_per_year = 10L)
  return(lmod)
 #qs::qsave(lmod, file = paste0(output.dir, "/lmod/scenario_", scenario, "/",iter,".qs"))


}


get_index <- function(x,index){
  y <- x[index]
  names(y) <- names(x)
  return(y)
}


sub_hivp <- function(input.list, hivp){
  hivp[names(input.list)] <- input.list
  return(hivp)
}


run_decomp <- function(run_name, scenario, uncertainty_params, input_table){

    root <- getwd()

    output_dir <- paste0(root, "/runs", run_name)
    lapply(c(output_dir, paste0(output_dir, '/lmod/')), dir.create, recursive = T)

    demp <- readRDS(paste0(root, "/demp.RDS"))
    hivp <- readRDS(paste0(root, "/hivp.RDS"))
    inc_rate <- readRDS(paste0(root, "/inc_rate_draws.RDS"))
    mat_prev <- readRDS(paste0(root, "/hiv_births_draws.RDS"))
    mat_prev_sd <- apply(mat_prev[, -(1:2)], MARGIN = 2, FUN = sd)
    out_prop_gte350 <- readRDS(paste0(root, "/prop_gte350_draws.RDS"))
    out_prop_lt200 <- readRDS(paste0(root, "/prop_lt200_draws.RDS"))
    out_prop_gte350 <- data.table::dcast(out_prop_gte350, iteration ~ year, value.var = 'prop_gte350')
    out_prop_lt200 <- data.table::dcast(out_prop_lt200, iteration ~ year, value.var = 'prop_lt200')
    stored_hivp_mtct <- hivp$mtct
    stored_hivp_pmtct_mtct <- hivp$pmtct_mtct
    stored_hivp_ctx_effect <- hivp$ctx_effect
    stored_hivp_paed_art_val <- hivp$paed_art_val
    saved_hivp <- hivp

    input <- input_table

    unparms <- paste0(output_dir, '/uncertainty_files/', uncertainty_params, '/')
    if(run_name %in% c('/03_new_sd_mtct/')){
    #if(run_name %in% c('/03_new_sd_mtct/', '/04_all_hivp_sources/')){
      unparms <- unparms[1:4]
    }

    input.x <- list()
    for(dir in  unparms){
      # assign(paste0(uncertainty_params[which(dir == unparms)], '_list'), lapply(paste0(dir,1:300, '.qs'), qs::qread))
      input.x[[which(dir == unparms)]] <- lapply(paste0(dir,1:300, '.qs'), qs::qread)
      print(dir)
    }
    names(input.x) <- uncertainty_params

    if(!dir.exists(paste0(output_dir, "/lmod/scenario_", scenario))){
      dir.create(paste0(output_dir, "/lmod/scenario_", scenario))
    }


    lmod <- lapply(1:300, run_lmod_new_input,  hivp_x = saved_hivp, demp_x = demp, scenario_row = input[scenario,],
           scenario = scenario, input_list = input.x,
           mat.prev = mat_prev, gte350 = out_prop_gte350, lt200 = out_prop_lt200,
           output.dir = output_dir)


    aidsdeaths_art_paed <- lapply(lmod, function(x) as.array(apply(x$aidsdeaths_art_paed, MARGIN = c(3,4,5), FUN = sum)))
    aidsdeaths_noart_paed <- lapply(lmod, function(x) as.array(apply(x$aidsdeaths_noart_paed, MARGIN = c(3,4,5), FUN = sum)))
    incidence <- lapply(lmod, function(x) x$infections[1:5,,])
    artpop <- lapply(lmod, function(x) x$artstrat_paeds)
    hivpop <- lapply(lmod, function(x) x$hivstrat_paeds)
    totpop <- lapply(lmod, function(x) x$totpop1[1:15,,])

    deaths_art <- data.table(melt(aidsdeaths_art_paed))
    deaths_art <- deaths_art[,.(age = Var1 - 1, sex = ifelse(Var2 == 1, 'male', 'female'), year = Var3 + 1969, draw = L1, value)]
    deaths_noart <- data.table(melt(aidsdeaths_noart_paed))
    deaths_noart <- deaths_noart[,.(age = Var1 - 1, sex = ifelse(Var2 == 1, 'male', 'female'), year = Var3 + 1969, draw = L1, value)]
    deaths_x <- merge(deaths_noart, deaths_art, by = c('age', 'sex', 'year', 'draw'))
    deaths_x <- deaths_x[,.(value = value.x + value.y), by = c('age', 'sex',  'year', 'draw')]

    hivpop <- data.table(melt(hivpop))
    hivpop <- hivpop[, .(age = Var3 - 1, sex = ifelse(Var4 == 1, "male", "female"), year = Var5 + 1969, draw = L1, value)]
    hivpop <- hivpop[, .(value = sum(value)), by = c("age", "year", "draw")]

    totpop <- data.table(melt(totpop))
    totpop <- totpop[, .(age = Var1 - 1, sex = ifelse(Var2 == 1, "male", "female"), year = Var3 + 1969, draw = L1, value)]
    totpop <- totpop[, .(value = sum(value)), by = c("age", "year", "draw")]

    hivnpop <- merge(hivpop, totpop, by = c("age", "year", "draw"))
    hivnpop <- hivnpop[, .(age, year, draw, value = value.y - value.x)]

    inf <- data.table(melt(incidence))
    inf <- inf[, .(age = Var1 - 1, sex = ifelse(Var2 == 1, "male", "female"), year = Var3 + 1969, draw = L1, value)]
    inf <- inf[, .(value = sum(value)), by = c("age", "year", "draw")]
    inf_rate <- merge(inf, hivnpop, by = c("year", "draw", "age"))
    inf_rate <- inf_rate[, .(value = value.x / value.y, year, draw, age)]

    death_rate <- merge(deaths_x, totpop, by = c("age", "year", 'draw'), allow.cartesian = T)
    death_rate <- death_rate[, .(value = value.x / value.y, age, year, draw)]

    out_sd <- list(
      deaths = death_rate,
      deaths_art = deaths_art,
      deaths_noart = deaths_noart,
      inf = inf_rate,
      hivpop = hivpop,
      hivnpop = hivnpop
    )
    qs::qsave(out_sd, file = paste0(root, '/runs/', run_name, "/scenario_", scenario, ".qs"))



}

##Not using any of these
extract_lmod <- function(run_name.x, scenario.x, uncertainty_params.x, input_table.x){
  root <- getwd()
  #lmod <- lapply(paste0(root, '/runs/', run_name, '/lmod/scenario_', scenario, '/', 1:300, '.qs'), qs::qread)
  lmod <-  run_decomp(run_name = run_name.x,
                       scenario = scenario.x,
                       uncertainty_params = uncertainty_params.x,
                       input_table = input_table.x)

  aidsdeaths_art_paed <- lapply(lmod, function(x) as.array(apply(x$aidsdeaths_art_paed, MARGIN = c(3,4,5), FUN = sum)))
  aidsdeaths_noart_paed <- lapply(lmod, function(x) as.array(apply(x$aidsdeaths_noart_paed, MARGIN = c(3,4,5), FUN = sum)))
  incidence <- lapply(lmod, function(x) x$infections[1:5,,])
  artpop <- lapply(lmod, function(x) x$artstrat_paeds)
  hivpop <- lapply(lmod, function(x) x$hivstrat_paeds)
  totpop <- lapply(lmod, function(x) x$totpop1[1:15,,])

  deaths_art <- data.table(melt(aidsdeaths_art_paed))
  deaths_art <- deaths_art[,.(age = Var1 - 1, sex = ifelse(Var2 == 1, 'male', 'female'), year = Var3 + 1969, draw = L1, value)]
  deaths_noart <- data.table(melt(aidsdeaths_noart_paed))
  deaths_noart <- deaths_noart[,.(age = Var1 - 1, sex = ifelse(Var2 == 1, 'male', 'female'), year = Var3 + 1969, draw = L1, value)]
  deaths_x <- merge(deaths_noart, deaths_art, by = c('age', 'sex', 'year', 'draw'))
  deaths_x <- deaths_x[,.(value = value.x + value.y), by = c('age', 'sex',  'year', 'draw')]

  hivpop <- data.table(melt(hivpop))
  hivpop <- hivpop[, .(age = Var3 - 1, sex = ifelse(Var4 == 1, "male", "female"), year = Var5 + 1969, draw = L1, value)]
  hivpop <- hivpop[, .(value = sum(value)), by = c("age", "year", "draw")]

  totpop <- data.table(melt(totpop))
  totpop <- totpop[, .(age = Var1 - 1, sex = ifelse(Var2 == 1, "male", "female"), year = Var3 + 1969, draw = L1, value)]
  totpop <- totpop[, .(value = sum(value)), by = c("age", "year", "draw")]

  hivnpop <- merge(hivpop, totpop, by = c("age", "year", "draw"))
  hivnpop <- hivnpop[, .(age, year, draw, value = value.y - value.x)]

  inf <- data.table(melt(incidence))
  inf <- inf[, .(age = Var1 - 1, sex = ifelse(Var2 == 1, "male", "female"), year = Var3 + 1969, draw = L1, value)]
  inf <- inf[, .(value = sum(value)), by = c("age", "year", "draw")]
  inf_rate <- merge(inf, hivnpop, by = c("year", "draw", "age"))
  inf_rate <- inf_rate[, .(value = value.x / value.y, year, draw, age)]

  death_rate <- merge(deaths_x, totpop, by = c("age", "year", 'draw'), allow.cartesian = T)
  death_rate <- death_rate[, .(value = value.x / value.y, age, year, draw)]

  out_sd <- list(
    deaths = death_rate,
    deaths_art = deaths_art,
    deaths_noart = deaths_noart,
    inf = inf_rate,
    hivpop = hivpop,
    hivnpop = hivnpop
  )
  qs::qsave(out_sd, file = paste0(root, '/runs/', run_name, "/scenario_", scenario, ".qs"))

}

unlink_scenario_files <- function(run_name, scenario){
  root <- getwd()

  dir_delete <-  paste0(root, "/runs/", run_name, "/lmod/scenario_", scenario)
  do.call(unlink,list(list.files(dir_link,full.names=TRUE)))

 #unlink(dir_delete, recursive = T)
}

run_total <- function(run_name, scenario, uncertainty_params, input_table) {

  run_decomp(run_name = run_name,
             scenario = scenario,
             uncertainty_params = uncertainty_params,
             input_table = input_table)
  extract_lmod(run_name = run_name,
               scenario = scenario)
  unlink_scenario_files(run_name = run_name,
                       scenario = scenario)

}

