library(simcausal)
library(rpart)
library(purrr)
library(tidyverse)
library(modeest)
theme_set(theme_minimal())

# function defined to compute contingency table for all combos of all binary variables
# and compute exposure probability and sample size
make_strata_table <- function(dat, A = "A", binary_vars){
  
  # create all combos of binary_vars of length >= 1
  subsets <- unlist(lapply(1:length(binary_vars), function(x) combn(binary_vars, x, simplify = FALSE)),
                    recursive = FALSE)
  
  # summary stats for each combo
  results <- map(subsets, function(vars) {
    dat %>%
      group_by(across(all_of(vars))) %>%
      summarise(n          = n(),
                n_treat    = sum(.data[[A]] == 1),
                n_control  = sum(.data[[A]] == 0),
                proba_exp  = mean(.data[[A]] == 1),
                sample_prop = n()/nrow(dat),
                .groups = "drop") %>%
      mutate(vars = paste(vars, collapse = ","))
  })
  
  # bind all results together
  bind_rows(results) %>%
    arrange(vars, proba_exp)
}
