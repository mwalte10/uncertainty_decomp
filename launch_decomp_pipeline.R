# library(styler)
working_dir <- "C:/Users/mwalters/Documents/Projects/paediatric_uncertainty/workflow/"
# style_dir(working_dir)
library(hipercow)
root <- getwd()
root <- "Z:/uncertainty_decomp"
hipercow_init(driver = "windows")
hipercow_provision()
resources <- hipercow_resources(cores = 2, )
options(hipercow.max_size_local = Inf)
source_files <- c('02_prep_uncertain_parms.R', '03_setup_decomp_test.R', '04_compile_results.R')
hipercow_environment_create(sources = source_files)
lapply(source_files, source)
source(paste0(root, '/central_functions/gets_weighted_pairs.R'))

run_name = '/04_all_hivp_sources/'

run_table <- data.table::fread('run_table.csv')
uncertainty_params <- run_table[run == run_name,uncertainty_vars]
uncertainty_params <- unlist(strsplit(uncertainty_params, split = ' '))
input <- lapply(1, get_weighted_pairs, var_vec = uncertainty_params)
input <- input[[1]]$total_combos
input[,index := NULL]


# system(paste0("Rscript --vanilla ", working_dir, "00_create_base_data.R"))
# ##Generate baseline mean values
# system(paste0("Rscript --vanilla ", working_dir, "/runs/", run_name, "01_generate_uncertainty.R ", run_name, " ", uncertainty_params))
# system(paste0("Rscript --vanilla ", working_dir, "/runs/", run_name, "01b_plot_uncertainty.R ", run_name))
##Decomposition starts here
##02_prep_uncertain_params.R
# id <- task_create_expr(prep_input_draws(run_name = '/04_all_hivp_sources/',uncertainty_params = uncertainty_params), resources = resources)
# task_status(id)
# task_result(id)

##03_setup_decomp_test.R
args <- 1:nrow(input)
batches <- ceiling(length(args) / 25)
for(batch in 4:batches){
  index <- 1:25 + (batch-1) * 25
  index <- index[index <= length(args)]
  bundle <- task_create_bulk_call(fn = run_total, data = args[index], args = list(run_name = '/04_all_hivp_sources/',
                                                                                  uncertainty_params = uncertainty_params,
                                                                                  input_table = input), resources = resources)
  hipercow_bundle_wait(bundle)
}
args <- 106:nrow(input)
bundle <- task_create_bulk_call(fn = extract_lmod, data = args[1], args = list(run_name = '/04_all_hivp_sources/',
                                                                       uncertainty_params = uncertainty_params,
                                                                      input.table = input), resources = resources)

bundle <- task_create_bulk_call(fn = extract_lmod, data = args[1], args = list(run_name.x = '/04_all_hivp_sources/',
                                                                               uncertainty_params.x = uncertainty_params,
                                                                               input_table.x = input), resources = resources)
hipercow_bundle_wait(bundle)
hipercow_bundle_result(bundle)


##04_compile_results.R
# id <- task_create_expr(run_compilation(run_name = '/04_all_hivp_sources/',uncertainty_params = uncertainty_params))
# task_status(id)
# task_result(id)


































# library(benchmark)
# bench::mark()
# profvis::profvis(simplify = FALSE)
# install.packages('jointprof')
# jointprof::
