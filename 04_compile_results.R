library(data.table)
source(paste0(root, '/central_functions/recreate_uncertainty_functions.R'))
source(paste0(root, '/central_functions/gets_weighted_pairs.R'))

## decomp pipeline is located here: C:\Users\mwalters\ihme-imperial\pipeline
extract_lmod <- function(run_name, scenario){
  root <- getwd()
  lmod <- lapply(paste0(root, '/runs/', run_name, '/lmod/scenario_', scenario, '/', 1:300, '.qs'), qs::qread)
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

get_pair <- function(var, w_pairs, iter, run_name) {
  w_var_pairs <- w_pairs[variable == var, ]
  a <- qs::qread(paste0(root, '/runs/', run_name, "/scenario_", w_var_pairs[pair == iter, `standard var index`], ".qs"))
  b <- qs::qread(paste0(root, '/runs/', run_name, "/scenario_", w_var_pairs[pair == iter, `no var index`], ".qs"))

  return(list(a, b))
}

extract_deaths <- function(a_deaths, b_deaths, weight) {
  a_deaths <- dcast(a_deaths, age + year ~ draw, value.var = "value")
  keep <- a_deaths[, 1:2]
  range <- apply(a_deaths[, 3:102], MARGIN = 1, FUN = function(x) quantile(x, 0.975, na.rm = T) - quantile(x, 0.025, na.rm = T))
  a_deaths <- cbind(keep, range)

  b_deaths <- dcast(b_deaths, age + year ~ draw, value.var = "value")
  keep <- b_deaths[, 1:2]
  range <- apply(b_deaths[, 3:102], MARGIN = 1, FUN = function(x) quantile(x, 0.975, na.rm = T) - quantile(x, 0.025, na.rm = T))
  b_deaths <- cbind(keep, range)

  deaths_pair <- merge(a_deaths, b_deaths, by = c("age", "year"))
  deaths_pair[, range_diff := (range.x - range.y) * weight]

  return(deaths_pair)
}

extract_deaths_art_spec <- function(a_deaths, b_deaths, weight) {
  a_deaths <- a_deaths[,.(value = sum(value)), by = c('age', 'year', 'draw')]
  a_deaths <- dcast(a_deaths, age + year ~ draw, value.var = "value")
  keep <- a_deaths[, 1:2]
  range <- apply(a_deaths[, 3:102], MARGIN = 1, FUN = function(x) quantile(x, 0.975, na.rm = T) - quantile(x, 0.025, na.rm = T))
  a_deaths <- cbind(keep, range)

  b_deaths <- b_deaths[,.(value = sum(value)), by = c('age', 'year', 'draw')]
  b_deaths <- dcast(b_deaths, age + year ~ draw, value.var = "value")
  keep <- b_deaths[, 1:2]
  range <- apply(b_deaths[, 3:102], MARGIN = 1, FUN = function(x) quantile(x, 0.975, na.rm = T) - quantile(x, 0.025, na.rm = T))
  b_deaths <- cbind(keep, range)

  deaths_pair <- merge(a_deaths, b_deaths, by = c("age", "year"))
  deaths_pair[, range_diff := (range.x - range.y) * weight]

  return(deaths_pair)
}

extract_incidence <- function(a_inf, b_inf, weight) {
  a_inf <- dcast(a_inf, age + year ~ draw, value.var = "value")
  keep <- a_inf[, 1:2]
  range <- apply(a_inf[, 3:102], MARGIN = 1, FUN = function(x) quantile(x, 0.975, na.rm = T) - quantile(x, 0.025, na.rm = T))
  a_inf <- cbind(keep, range)

  b_inf <- dcast(b_inf, age + year ~ draw, value.var = "value")
  keep <- b_inf[, 1:2]
  range <- apply(b_inf[, 3:102], MARGIN = 1, FUN = function(x) quantile(x, 0.975, na.rm = T) - quantile(x, 0.025, na.rm = T))
  b_inf <- cbind(keep, range)

  inf_pair <- merge(a_inf, b_inf, by = c("age", "year"))
  inf_pair[, range_diff := (range.x - range.y) * weight]

  return(inf_pair)
}

extract_prevalence <- function(a_hivpop, b_hivpop, a_hivnpop, b_hivnpop, weight) {
  a_hivpop <- merge(a_hivpop, a_hivnpop, by = c("age", "year", "draw"))
  a_hivpop <- a_hivpop[, value := value.x / (value.x + value.y)]
  a_hivpop <- a_hivpop[, .(age, year, draw, value)]

  b_hivpop <- merge(b_hivpop, b_hivnpop, by = c("age", "year", "draw"))
  b_hivpop <- b_hivpop[, value := value.x / (value.x + value.y)]
  b_hivpop <- b_hivpop[, .(age, year, draw, value)]

  a_hivpop <- dcast(a_hivpop, age + year ~ draw, value.var = "value")
  keep <- a_hivpop[, 1:2]
  range <- apply(a_hivpop[, 3:102], MARGIN = 1, FUN = function(x) quantile(x, 0.975, na.rm = T) - quantile(x, 0.025, na.rm = T))
  a_hivpop <- cbind(keep, range)

  b_hivpop <- dcast(b_hivpop, age + year ~ draw, value.var = "value")
  keep <- b_hivpop[, 1:2]
  range <- apply(b_hivpop[, 3:102], MARGIN = 1, FUN = function(x) quantile(x, 0.975, na.rm = T) - quantile(x, 0.025, na.rm = T))
  b_hivpop <- cbind(keep, range)

  hivpop_pair <- merge(a_hivpop, b_hivpop, by = c("age", "year"))
  hivpop_pair[, range_diff := (range.x - range.y) * weight]

  return(hivpop_pair)
}

extract_total_effects <- function(var_x, pair_table, iteration,
                                  deaths_dt, hivpop_dt, inf_dt,
                                  deaths_noart_dt, deaths_art_dt,
                                  run_name) {
  pairs <- get_pair(var = var_x, w_pairs = pair_table[variable == var_x, ], iter = iteration, run_name)
  a <- pairs[[1]]
  b <- pairs[[2]]
  a$deaths <- a$deaths[year > 1989]
  a$inf <- a$inf[year > 1989, ]
  a$hivpop <- a$hivpop[year > 1989, ]
  a$hivnpop <- a$hivnpop[year > 1989, ]
  b$deaths <- b$deaths[year > 1989]
  b$inf <- b$inf[year > 1989, ]
  b$hivpop <- b$hivpop[year > 1989, ]
  b$hivnpop <- b$hivnpop[year > 1989, ]

  ### DEATHS
  deaths_pair <- extract_deaths(a$deaths, b$deaths, weight = pair_table[iteration, weight])

  ### INCIDENCE
  inf_pair <- extract_incidence(a$inf, b$inf, weight = pair_table[iteration, weight])

  ### PREVALENCE
  hivpop_pair <- extract_prevalence(a$hivpop, b$hivpop,
    a$hivnpop, b$hivnpop,
    weight = pair_table[iteration, weight]
  )

  ### DEATHS NO ART
  deaths_noart_pair <- extract_deaths_art_spec(a$deaths_noart, b$deaths_noart, weight = pair_table[iteration, weight])

  ### DEATHS ART
  deaths_art_pair <- extract_deaths_art_spec(a$deaths_art, b$deaths_art, weight = pair_table[iteration, weight])


  deaths <- merge(deaths_dt, deaths_pair[, .(age, year, range_diff)], by = c("year", "age"))

  inf <- merge(inf_dt, inf_pair[, .(age, year, range_diff)], by = c("year", "age"))

  hivpop <- merge(hivpop_dt, hivpop_pair[, .(age, year, range_diff)], by = c("year", "age"))

  deaths_noart <- merge(deaths_noart_dt, deaths_noart_pair[, .(age, year, range_diff)], by = c("year", "age"))

  deaths_art <- merge(deaths_art_dt, deaths_art_pair[, .(age, year, range_diff)], by = c("year", "age"))


  out <- rbind(deaths, inf, hivpop,
               deaths_noart, deaths_art, fill = T)

  return(out)
}

decompose_effects <- function(pair_table, var, run_name) {
  y <- pair_table[variable == var]

  effects <- lapply(1:nrow(y), extract_total_effects,
                    run_name = run_name,
    var_x = var, pair_table = pair_table,
    deaths_dt = data.table(expand.grid(list(year = 1970:2022, age = 0:2, range_diff = 0, metric = "deaths"))),
    hivpop_dt = data.table(expand.grid(list(year = 1970:2022, age = 0:2, range_diff = 0, metric = "hivpop"))),
    inf_dt = data.table(expand.grid(list(year = 1970:2022, age = 0:4, range_diff = 0, metric = "inf"))),
    deaths_noart_dt = data.table(expand.grid(list(year = 1970:2022, age = 0:4, range_diff = 0, metric = "deaths_noart"))),
    deaths_art_dt = data.table(expand.grid(list(year = 1970:2022, age = 0:4, range_diff = 0, metric = "deaths_art")))
  )

  effects <- rbindlist(effects)

  effects <- effects[, .(range_diff = sum(range_diff.x, na.rm = T) + sum(range_diff.y, na.rm = T), var = var), by = c("year", "age", "metric")]

  return(effects)
}

save_baseline_results <- function(run_name){
  baseline_in <- qs::qread(paste0(root, '/runs/', run_name, "/scenario_1.qs"))
  baseline <- merge(baseline_in$hivpop[year < 2023], baseline_in$hivnpop[year < 2023], by = c("age", "year", "draw"))
  baseline[, pop := value.x + value.y]
  baseline[, value := value.x / pop]
  baseline <- baseline[, .(age, year, draw, value)]
  baseline <- baseline[, .(lower = quantile(value, 0.025, na.rm = T), upper = quantile(value, 0.975, na.rm = T), value = median(value, na.rm = T)), by = c("age", "year")]
  baseline[, IQR := upper - lower]
  baseline_prev <- baseline[, metric := "hivpop"]

  baseline <- baseline_in$inf[year < 2023]
  baseline <- baseline[year > 1989]
  baseline <- baseline[, .(value = median(value, na.rm = T), lower = quantile(value, 0.025, na.rm = T), upper = quantile(value, 0.975, na.rm = T)), by = c("age", "year")]
  baseline[, IQR := upper - lower]
  baseline_inf <- baseline[, metric := "inf"]

  baseline <- baseline_in$deaths[year < 2023]
  baseline <- baseline[year > 1989]
  baseline <- baseline[, .(value = median(value, na.rm = T), lower = quantile(value, 0.025, na.rm = T), upper = quantile(value, 0.975, na.rm = T)), by = c("age", "year")]
  baseline[, IQR := upper - lower]
  baseline_deaths <- baseline[, metric := "deaths"]

  baseline <- baseline_in$deaths_noart
  baseline <- baseline[year > 1989 & year < 2023]
  baseline <- baseline[, .(value = median(value, na.rm = T), lower = quantile(value, 0.025, na.rm = T), upper = quantile(value, 0.975, na.rm = T)), by = c("age", "year")]
  baseline[, IQR := upper - lower]
  baseline_deaths_noart <- baseline[, metric := "deaths_noart"]

  baseline <- baseline_in$deaths_art
  baseline <- baseline[year > 1989 & year < 2023]
  baseline <- baseline[, .(value = median(value, na.rm = T), lower = quantile(value, 0.025, na.rm = T), upper = quantile(value, 0.975, na.rm = T)), by = c("age", "year")]
  baseline[, IQR := upper - lower]
  baseline_deaths_art <- baseline[, metric := "deaths_art"]

  baseline <- rbind(baseline_prev, baseline_inf, baseline_deaths, baseline_deaths_noart, baseline_deaths_art)
  qs::qsave(baseline, file = paste0(root, '/runs/', run_name, "/baseline_uncertainty.qs"))
}

run_compilation <- function(run_name, uncertainty_params){

  input <- lapply(1, get_weighted_pairs, var_vec = uncertainty_params)
  weight_pairs <- unique(data.table(rbindlist(input[[1]]$weight_pairs)))

  lapply(unique(c(weight_pairs$`standard var index`, weight_pairs$`no var index`)), extract_lmod, run_name = run_name)

  save_baseline_results(run_name)

  out <- lapply(unique(weight_pairs$variable), decompose_effects, pair_table = weight_pairs, run_name = run_name)
  out <- rbindlist(out)
  out[range_diff < 0, range_diff := 0]
  out[, total := sum(range_diff, na.rm = T), by = c("year", "age", "metric")]
  out[, prop := range_diff / total]
  out <- out[, .(year, age, metric, var, prop)]
  qs::qsave(out, file = paste0(root, '/runs/', run_name, "/decomposed_effects.qs"))

}





