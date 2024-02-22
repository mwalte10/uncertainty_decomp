## ----------------------------------------------------------------------------
## get_weighted_pairs.R
## Purpose: Get all combinations of Shapley pairings for decomposition, and their eventual weighting within the decomposition.
## Can be used to generate all of the combinations that need to be run for the decomposition (total_combos in output)
## Also used to get the combinations needed for a given variable to decompose, as well as the weights used

#####Function inputs
############## index: index in the following var_vec that will be used to get the shapley pairs and their corresponding weights
############## var_vec: variables that represent the different modalities of uncertainty

##Helpful resources: https://towardsdatascience.com/shap-explained-the-way-i-wish-someone-explained-it-to-me-ab81cc69ef30


## ----------------------------------------------------------------------------
##make function to pull out pairs from the different combinations and weight them
get_weighted_pairs <- function(index = 1, var_vec = c('mat_prev', 'mtct_rate', 'cotrim', 'numberonart_scalar')){

  ## Create index mat
  input = vector('list', length = length(var_vec))
  names(input) <- var_vec
  for(i in 1:length(var_vec)){
    ##1: standard variation
    ##0: no variation
    input[[i]] <- c(1,0)
  }
  input <- expand.grid(input)
  input <- data.table::data.table(input)
  input[,index := 1:nrow(input)]
  assertthat::assert_that(nrow(input) == 2^(length(var_vec)), msg = 'Expected number of scenarios are created')

  pairs <- lapply(1:length(var_vec), function(index){
    var = var_vec[index]
    stvar <- input[get(var) == 1,]
    no_var <- input[get(var) == 0,]
    stvar[,pair := paste0(1:nrow(stvar))]
    no_var[,pair := paste0(1:nrow(stvar))]

    pairwise <- merge(stvar[,.(`standard var index` = index, pair)], no_var[,.(`no var index` = index, pair)], by = c('pair'))
    pairwise[,pair := as.integer(pair)]

    ##get weights, for ease of calculation pull these from the inputs list where the variable is being altered
    weights = input[get(var) == 0,]
    weights[,index := NULL]
    weights[,c := rowSums(weights)]
    weights[,n := length(var_vec)]
    weights[,weight := (factorial(c) * factorial(n - c - 1)) / factorial(n)]
    assertthat::assert_that(sum(weights$weight) == 1, msg = 'Total weights for this variable are equal to one')

    weights[,pair := 1:nrow(weights)]

    ##Attach weights to pairwise dt
    pairwise <- merge(pairwise, weights[,.(pair, weight)], by = 'pair')
    pairwise[,variable := var]
  })


  return(list(total_combos = input, weight_pairs = pairs))

}

##loop through each of the variables that are being varied here
x = lapply(1:4, get_weighted_pairs, var_vec = c('mat_prev', 'mtct_rate', 'cotrim', 'numberonart_scalar'))


