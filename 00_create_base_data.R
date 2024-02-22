#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 00_create_base_data.R
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

rm(list = ls())
library(eppasm)
devtools::load_all("C:/Users/mwalters/frogger/")
library(leapfrog)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(eppasm)

args <- commandArgs(trailingOnly = TRUE)
run_name <- args[1]
uncertainty_params <- args[2]
uncertainty_params <- unlist(strsplit(uncertainty_params, split = " "))
print(run_name)

root <- getwd()
working_dir <- "Z:/uncertainty_decomp/from_spectrum/"
dir.create(working_dir, recursive = T)

pjnz1 <- paste0(working_dir, "/Botswana2023v4 WPP 02_03_2023 KOS.PJNZ")
specres <- eppasm::read_hivproj_output(pjnz1)
hiv_births <- readxl::read_xlsx(paste0(working_dir, "/uncertainty_analysis_300 draws.XLSX"), sheet = "101. Mothers needing PMTCT")
inc <- readxl::read_xlsx(paste0(working_dir, "/uncertainty_analysis_300 draws.XLSX"), sheet = "52. New HIV infections (15-49)")
pop <- readxl::read_xlsx(paste0(working_dir, "/uncertainty_analysis_300 draws.XLSX"), sheet = "8. Population aged 15-49")
pop <- data.table(pop)
inc <- data.table(inc)
colnames(pop) <- unname(unlist(pop[2, ]))
pop <- pop[-(1:2), ]
colnames(inc) <- unname(unlist(inc[2, ]))
inc <- inc[-(1:2), ]
inc <- melt(inc, id.vars = c("Iteration", "Sex"))
setnames(inc, "value", "inc")
pop <- melt(pop, id.vars = c("Sex", "Iteration"))
setnames(pop, "value", "pop")
inc_rate <- merge(inc, pop, by = c("Iteration", "Sex", "variable"))
inc_rate <- inc_rate[Sex == "Male+Female"]
inc_rate[, inc_rate := inc / pop]

fp1 <- eppasm::prepare_directincid(pjnz1)
fp1$tARTstart <- 61L
out_prop_gte350 <- data.table()
out_prop_lt200 <- data.table()
out_hiv_births <- data.table()

for (i in 1:300) {
  fp1$incidinput <- inc_rate[Iteration == i, inc_rate]
  mod1 <- eppasm::simmod(fp1)
  prop_gte350 <- apply(attr(mod1, "hivpop")[1:2, 1:8, 2, ], MARGIN = 3, FUN = sum) / apply(attr(mod1, "hivpop")[, 1:8, 2, ], MARGIN = 3, FUN = sum)
  prop_lt200 <- apply(attr(mod1, "hivpop")[6:7, 1:8, 2, ], MARGIN = 3, FUN = sum) / apply(attr(mod1, "hivpop")[, 1:8, 2, ], MARGIN = 3, FUN = sum)
  out_prop_gte350 <- rbind(out_prop_gte350, data.table(prop_gte350, year = 1970:2030, iteration = i))
  out_prop_lt200 <- rbind(out_prop_lt200, data.table(prop_lt200, year = 1970:2030, iteration = i))
  print(i / 300)
}

out_prop_gte350[is.nan(prop_gte350), prop_gte350 := 0]
out_prop_lt200[is.nan(prop_lt200), prop_lt200 := 0]


hiv_births <- data.table(hiv_births)
colnames(hiv_births) <- unlist(hiv_births[2, ])
hiv_births <- hiv_births[3:nrow(hiv_births)]
comp <- eppasm::read_hivproj_output(pjnz1)


demp <- leapfrog::prepare_leapfrog_demp(pjnz1)
hivp <- leapfrog::prepare_leapfrog_projp(pjnz1)
hivp <- leapfrog:::prepare_hc_leapfrog_projp(pjnz1, hivp)
hivp$art_mtct[, , 2] <- hivp$art_mtct[, c(3, 2, 1), 2]
hivp$fert_mult_offart[] <- 1
hivp$fert_mult_by_age[] <- 1
#hivp$paed_art_val[which(hivp$artpaeds_isperc)] <- hivp$paed_art_val[which(hivp$artpaeds_isperc)] / 100
## pull in the mean values of the things that are used from the adult model
hivp$mat_hiv_births <- comp$hivpregwomen
hivp$prop_gte350 <- out_prop_gte350[, .(prop_gte350 = mean(prop_gte350)), by = "year"]$prop_gte350
hivp$prop_lt200 <- out_prop_lt200[, .(prop_lt200 = mean(prop_lt200)), by = "year"]$prop_lt200

params_list <- lapply(uncertainty_params, function(x) hivp[[x]])

names(params_list) <- uncertainty_params
# for (i in 1:length(uncertainty_params)) {
#   attr(params_list[[i]], "var") <- uncertainty_params[i]
# }

if(!dir.exists(paste0(root, '/runs/', run_name))){dir.create(paste0(root, '/runs/', run_name))}
saveRDS(uncertainty_params, paste0(root, '/runs/', run_name, "/uncertainty_params.RDS"))
#saveRDS(params_list, paste0(root, '/runs/', run_name, "params_list.RDS"))
saveRDS(demp, paste0("demp.RDS"))
saveRDS(hivp, paste0( "hivp.RDS"))
saveRDS(inc_rate, paste0("inc_rate_draws.RDS"))
saveRDS(hiv_births, paste0("hiv_births_draws.RDS"))
saveRDS(out_prop_lt200, paste0("prop_lt200_draws.RDS"))
saveRDS(out_prop_gte350, paste0("prop_gte350_draws.RDS"))
