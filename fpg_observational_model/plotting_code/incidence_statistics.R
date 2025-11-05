################################################################################
# Incidence bins statistical calculation 
################################################################################
IncidencePairwiseCalc <- function(data, y_vars, group_vars,
                                  adjust.method = "BH") {
  
  # build a temporary data frame with a "groups" column
  df <- tidyr::unite(data, col = "groups", all_of(group_vars), sep    = "_", remove = FALSE) 
  valid_groups <- df %>%
    group_by(groups) %>%
    filter(!is.na(y_vars[1])) %>%
    summarise(n = n()) %>%
    filter(n >= 2) %>%
    pull(groups)
  df <- dplyr::filter(df, groups %in% valid_groups)
  
  # for each y-variable, run a pairwise t-test vs groups
  purrr::map_dfr(y_vars, function(y) {
    df %>%
      rstatix::pairwise_t_test(
        formula         = as.formula(sprintf("`%s` ~ groups", y)),
        p.adjust.method = adjust.method
      ) %>%
      mutate(metric = y, .before = 1)
  })
}


IncidencePairwiseEdit <- function(df, bin_order = bin_order_vector, 
                                  group1_col = "group1", group1_split = c("group1_bin", "group1_var"), 
                                  group2_col = "group2", group2_split = c("group2_bin", "group2_var")){
  tmp_df <- df %>%
    # 1) split off bin + sampling for both group1 & group2
    tidyr::separate({{group1_col}}, into  = group1_split, sep = "_", extra = "merge") %>%
    tidyr::separate({{group2_col}}, into  = group2_split, sep = "_", extra = "merge") %>%
    dplyr::mutate(
      pair_bin = case_when(
        group1_bin == group2_bin ~ group1_bin,
        TRUE ~ "Different transmission bins"),
      pair_var = case_when(
        group1_var == group2_var ~ group1_var,
        TRUE  ~ "Different grouping variable")) %>%
    # reorder incidence bins for pairs 
    dplyr::mutate(
      bin_pos1 = match(group1_bin, bin_order),
      bin_pos2 = match(group2_bin, bin_order),
      # use pmin/pmax to pick the smaller/larger index:
      group1_bin_order = bin_order[pmin(bin_pos1, bin_pos2)],
      group2_bin_order = bin_order[pmax(bin_pos1, bin_pos2)]
    ) %>%
    # 4) turn into ordered factors
    dplyr::mutate(
      group1_bin_order = factor(group1_bin_order, levels = bin_order),
      group2_bin_order = factor(group2_bin_order, levels = bin_order)
    ) %>%
    dplyr::select(-bin_pos1, -bin_pos2)
  return(tmp_df)
}
