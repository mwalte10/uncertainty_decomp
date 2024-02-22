library(data.table)
library(ggplot2)
source(paste0(root, '/central_functions/recreate_uncertainty_functions.R'))
source(paste0(root, '/central_functions/gets_weighted_pairs.R'))

get_lower_and_upper <- function(dt_in, order = c("mat_hiv_births", "mtct", "ctx_effect", "paed_art_val"), band_type) {
  dt_in <- dt_in[, .(year, age, value, var, IQR_prop)]
  dt_in <- dcast(dt_in, year + age + value ~ var, value.var = "IQR_prop")

  for (i in 1:length(order)) {
    if (i == 1) {
      dt_in[, paste0(band_type, "__", order[i], "__a") := value]
      if (band_type == "upperband") {
        dt_in[, paste0(band_type, "__", order[i], "__b") := value + get(order[i])]
      } else {
        dt_in[, paste0(band_type, "__", order[i], "__b") := value - get(order[i])]
      }
    }

    if (i > 1) {
      dt_in[, paste0(band_type, "__", order[i], "__a") := get(paste0(band_type, "__", order[i - 1], "__b"))]
      if (band_type == "upperband") {
        dt_in[, paste0(band_type, "__", order[i], "__b") := get(paste0(band_type, "__", order[i], "__a")) + get(order[i])]
      } else {
        dt_in[, paste0(band_type, "__", order[i], "__b") := get(paste0(band_type, "__", order[i], "__a")) - get(order[i])]
      }
    }
  }

  dt_in <- dt_in[, !which(colnames(dt_in) %in% c("value", order)), with = FALSE]

  dt_in <- melt(dt_in, id.vars = c("age", "year"))
  dt_in[, var := unlist(lapply(strsplit(as.character(dt_in$variable), split = "__"), "[[", 2))]
  dt_in[, maxim := unlist(lapply(strsplit(as.character(dt_in$variable), split = "__"), "[[", 3))]
  dt_in <- dcast(dt_in[, .(year, age, value, var, maxim)], year + age + var ~ maxim, value.var = "value")

  dt_in[, band := gsub(band_type, pattern = "band", replacement = "")]


  return(dt_in)
}

get_iqr_prop <- function(metric_x = "hivpop", dt_decomp, dt_baseline) {
  dt_deaths <- unique(dt_decomp[metric == metric_x, .(metric, year, age, var, prop)])
  dt_baseline <- unique(dt_baseline[metric == metric_x, .(metric, year, value, age, IQR)])
  dt_baseline <- merge(dt_baseline, dt_deaths, by = c("year", "age", "metric"))
  dt_baseline[, IQR_prop := (IQR / 2) * prop]

  return(dt_baseline)
}

prep_for_plotting <- function(metric_x = "hivpop", dt_decomp, dt_baseline, order_vec) {
  dt_iqr <- get_iqr_prop(metric_x, dt_decomp, dt_baseline)

  ## test is now upper
  upper <- get_lower_and_upper(dt_in = dt_iqr, band_type = "upperband", order = order_vec)
  lower <- get_lower_and_upper(dt_in = dt_iqr, band_type = "lowerband", order = order_vec)

  return(list(upper, lower))
}

###########################################
### Plot observed uncertainty
###########################################
# dev.new(width = 9, height = 4.5,noRStudioGD = T)

plot_trend_uncertainty <- function(metric, baseline_dt, decomposed_dt, yaxis_title){
  root <- getwd()
  var_map <- fread(paste0(root, '/central_functions/var_name_map.csv'))
  dt <- rbindlist(prep_for_plotting(metric_x = metric, dt_decomp = decomposed_dt, dt_baseline = baseline_dt,
                                    order_vec = order_vec))
  dt <- merge(dt, var_map, by = "var")
  dt$var_name <- factor(dt$var_name, levels = var_map[var %in% order_vec, var_name])
  dt[band == "lower" & a < 0, a := 0]
  dt[band == "lower" & b < 0, b := 0]

  gg <- ggplot() +
    geom_ribbon(data = dt[band == "upper" & year %in% 1990:2022 & age < 3, ], aes(x = year, ymin = a, ymax = b, fill = as.factor(var_name), color = NULL), alpha = 0.8, show.legend = FALSE) +
    geom_ribbon(data = dt[band == "lower" & year %in% 1990:2022 & age < 3, ], aes(x = year, ymin = b, ymax = a, fill = as.factor(var_name), color = NULL), alpha = 0.8, show.legend = FALSE) +
    geom_line(data = baseline_dt[year %in% 1990:2022 & age < 3 & metric == metric], aes(year, value)) +
    facet_wrap(~age) +
    labs(
      x = NULL, y = yaxis_title,
      fill = "Input uncertainty"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom", strip.text = element_text(size = 13), plot.title = element_text(size = 16),
      axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12),
      panel.background = element_rect(fill = NA, color = "black")
    )

  return(gg)
}

plot_geom_area_uncertainty <- function(metric, baseline_dt, decomposed_dt){
  root <- getwd()
  var_map <- fread(paste0(root, '/central_functions/var_name_map.csv'))
  baseline_iqr <- get_iqr_prop(metric_x = "deaths", dt_decomp = decomposed_dt, dt_baseline = baseline_dt)
  baseline_iqr <- baseline_iqr[, .(year, age, value, IQR_prop, var)]
  baseline_iqr <- merge(baseline_iqr, var_map, by = "var")
  baseline_iqr[, IQR_prop := IQR_prop / value]
  baseline_iqr$var_name <- factor(baseline_iqr$var_name, levels = var_map[var %in% order_vec, var_name])

  dt <- decomposed_dt[metric == "deaths"]
  dt[, total := sum(abs(prop)), by = c("age", "year")]
  dt[, prop := prop / total]
  dt <- merge(dt, var_map, by = "var")
  dt$var_name <- factor(dt$var_name, levels = var_map[var %in% order_vec, var_name])

  dt[is.nan(prop), prop := 0]

  gg <- ggplot() +
    geom_area(data = dt[year %in% 1990:2022 & age < 3], aes(year, prop, fill = as.factor(var_name))) +
    facet_wrap(~age, labeller = label_wrap_gen()) +
    labs(x = NULL, y = "Proportion of total uncertainty", fill = "Uncertainty source") +
    theme_bw() +
    theme_bw() +
    theme(
      legend.position = "bottom", strip.text = element_text(size = 13), plot.title = element_text(size = 16),
      axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 10),
      axis.text.x = element_text(angle = 90, size = 10),
      panel.background = element_rect(fill = NA, color = "black")
    ) +
    labs(x = "", y = "Proportion of uncertainty")

  return(gg)
}

plot_results <- function(baseline_dt, decomposed_dt, order_vec, metric){
  if(metric == 'deaths'){
    yaxis_title = "HIV mortality rate among ALL children"
  }
  if(metric == 'deaths_noart'){
    yaxis_title = "HIV deaths- no ART"
  }
  if(metric == 'deaths_art'){
    yaxis_title = 'HIV deaths- ART'
  }
  if(metric == 'inf'){
    yaxis_title =  "Incidence rate"
  }
  if(metric == 'hivpop'){
    yaxis_title = 'Prevalence'
  }
  gg1 <- plot_trend_uncertainty(metric = metric, baseline_dt, decomposed_dt, yaxis_title = yaxis_title)
  gg2 <- plot_geom_area_uncertainty(metric = metric, baseline_dt, decomposed_dt)

  gg3 <- ggpubr::ggarrange(gg1, gg2, nrow = 2)
  return(gg3)

}

save_plotted_results <- function(baseline_dt, decomposed_dt, order_vec, run_name){
  root <- getwd()
  out_dir <- paste0(root, '/runs/', run_name)

  pdf(paste0(out_dir, 'results.pdf'), width = 9, height = 9)
  plot_results(baseline_dt, decomposed_dt, order_vec, 'inf')
  plot_results(baseline_dt, decomposed_dt, order_vec, 'hivpop')
  plot_results(baseline_dt, decomposed_dt, order_vec, 'deaths')
  plot_results(baseline_dt, decomposed_dt, order_vec, 'deaths_noart')
  plot_results(baseline_dt, decomposed_dt, order_vec, 'deaths_art')
  dev.off()

}

