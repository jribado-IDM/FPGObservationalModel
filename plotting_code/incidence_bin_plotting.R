################################################################################
# Set plotting colors
################################################################################
who_bin_order = c("Very low", "Low", "Moderate", "High")
transmission_colors <- setNames( 
  rev(c("#DE9B71", "#EFBC82", "#FBDFA2", "cornsilk")),
  who_bin_order)


################################################################################
# Trend plotting functions
################################################################################
ProportionBasePlot <- function(df, x_var = "reported_incidence_per_1k", y_var,
                               x_facet = "importation_rate", y_facet = "itn_on",
                               color_variable = "sampling_scheme",
                               shape_variable = NULL,
                               summarize_trend = FALSE){
  
  if(isFALSE(summarize_trend)){
    p <- ggplot(df, aes(x=get(x_var), y=get(y_var), 
                        color=get(color_variable), 
                        # shape=get(shape_variable)
    )) +
      geom_point(alpha=0.1) +
      #facet_grid(get(y_facet)~get(x_facet)) +
      labs(x=x_var, 
           y=y_var,
           color = color_variable) 
  } else{
    p <- df %>%
      ggplot(aes(x=get(x_var), y = get(y_var),
                 color=get(color_variable), 
                 fill=get(color_variable),
                 # shape=get(shape_variable),
      )) +
      geom_smooth(method = "lm", se = 0.68, alpha=0.5)  +
      geom_point(alpha=0.1) +
      labs(x=x_var, y=y_var,
           color = color_variable, 
           fill = color_variable) #+
    #facet_grid(get(y_facet)~get(x_facet)) 
  }
  
  if(color_variable == "sampling_scheme"){
    p <- p + 
      scale_color_manual(values=sampling_scheme_colors) +
      scale_fill_manual(values=sampling_scheme_colors) 
  } else if(color_variable == "season_group") {
    p <- p + 
      scale_color_manual(values=season_peak_colors) +
      scale_fill_manual(values=season_peak_colors)
  } else{
    p <- p + 
      scale_color_viridis_b() +
      scale_fill_viridis_b()
  }
  
  return(p)
} 



PlotBestFit <- function(df, y, x = "reported_incidence_per_1k",
                        span = 0.5,
                        color_variable = "sampling_scheme",
                        plot_points = TRUE) {
  
  predict_df <- df %>%
    mutate(
      .x = .data[[x]],
      .y = .data[[y]]
    )
  
  # fit candidate models
  mod_lin  <- lm(.y ~ .x,                     data = predict_df)
  mod_quad <- lm(.y ~ poly(.x, 2, raw = TRUE), data = predict_df)
  mod_cub  <- lm(.y ~ poly(.x, 3, raw = TRUE), data = predict_df)
  mod_lo   <- loess(.y ~ .x,                   data = predict_df, span = span)
  
  # RMSE helper (use predict_df)
  rmse <- function(mod) {
    pred <- predict(mod, newdata = predict_df)
    sqrt(mean((pred - predict_df$.y)^2, na.rm = TRUE))
  }
  
  stats <- tibble::tibble(
    model = c("linear", "quadratic", "cubic"), # , "loess"),
    RMSE  = c(rmse(mod_lin),
              rmse(mod_quad),
              rmse(mod_cub))
    #           rmse(mod_lo))
  )
  
  # choose the best
  best_name <- stats$model[which.min(stats$RMSE)]
  best_mod  <- switch(best_name,
                      linear    = mod_lin,
                      quadratic = mod_quad,
                      cubic     = mod_cub
                      #                   loess     = mod_lo
  )
  
  # build the plot, mapping color to your grouping variable
  p <- ggplot(predict_df, aes(x = .x, y = .y, 
                              color = .data[[color_variable]], fill = .data[[color_variable]])) +
    {
      if(isTRUE(plot_points)){
        geom_point(alpha = 0.4)    
      } 
    } +
    { 
      if (best_name == "loess") {
        geom_smooth(method = "loess",
                    span   = span,
                    se     = 0.68, alpha=0.2)
      } else {
        form <- switch(best_name,
                       linear    = y ~ x,
                       quadratic = y ~ poly(x, 2, raw = TRUE),
                       cubic     = y ~ poly(x, 3, raw = TRUE))
        geom_smooth(method  = "lm",
                    formula = form,
                    se      = 0.68, alpha=0.2)
      }
    } +
    scale_color_manual(values=sampling_scheme_colors) +
    scale_fill_manual(values=sampling_scheme_colors) +
    {
      if(grepl("prop", y)){
        ylim(c(0,1))
      }
    } #+
    #labs(subtitle = glue::glue("Best smoothing model fit: {best_name}\n(RMSE = {round(stats$RMSE[stats$model==best_name],4)})")) 
  
  list(
    model   = best_mod,
    metrics = stats,
    plot    = p
  )
}


# Create versions of PlotBest fit with limited information
SummarizeByGroups <- function(data,
                              group_cols,     # character vector of columns to group by
                              metric_col,     # single column name to summarise
                              fns = list(     # which statistics to compute
                                count = ~length(.x),
                                mean = ~mean(.x, na.rm=TRUE),
                                sd   = ~sd(  .x, na.rm=TRUE),
                                min  = ~min( .x, na.rm=TRUE),
                                max  = ~max( .x, na.rm=TRUE)
                              )) {
  data %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      across(
        all_of(metric_col),
        .fns   = fns,
        .names = "{fn}_{.col}"
      ),
      .groups = "drop"
    )
}


IncidenceCurveBuild <- function(df,
                                metric = input_metric,
                                incidence = "reported_incidence_per_1k",
                                y_lim = c(0,1),
                                center_lines = center_line_vector){
  
  p <- PlotBestFit(df, 
                   y=metric,
                   plot_points = FALSE)$plot +
    labs(x="Yearly mean reported incidence per 1k", 
         y = y_labs[names(y_labs) == metric],
         color = "Sampling scheme",
         fill = "Sampling scheme") +
    theme(
      #legend.position = "top") + 
      legend.position = c(0.28, 0.80),
      legend.background = element_rect(fill = "#F5F3ED")) +
    scale_x_log10() +
    ylim(y_lim)
  
  # Subset points along genetic metric lines
  points_df <- dplyr::bind_rows(lapply(center_lines, function(c){
    count_df <- dplyr::filter(df, between(.data[[metric]], c - 0.01, c + 0.01)) %>%
      dplyr::mutate(center_point = c)
    return(count_df)
  }))
  
  edges_long <- SummarizeByGroups(
      data         = points_df,
      group_cols   = c("center_point"),
      metric_col   = incidence         # or "incidence" if you prefer a string
    ) %>%
    pivot_longer(
      cols      = c(paste0("min_", incidence), paste0("max_", incidence)),
      names_to  = "which_edge",
      values_to = "incidence_x"
    )
  
  p_points <- p +
    geom_hline(data = tibble(center_line = center_lines),
               aes(yintercept = center_line),
               color = "grey50",
               linetype    = "dashed") +
    geom_point(data = points_df,
               aes(x=.data[[incidence]], y=.data[[metric]]), alpha=0.5) 
  
  
  p_segment <- p_points +
    geom_segment(
      inherit.aes = F,
      data = edges_long,
      aes(x = incidence_x, xend = incidence_x,
          y = center_point, yend = min(y_lim)),
      color = "grey",
      linetype = "dashed",
      size  = 0.75
    ) 

  
  plots_list = list(
    "TrendLineOnly" = p,
    "TrendLinePoints" = p_points,
    "TrendLineSegments" = p_segment
  )
  
  return(plots_list)  
} 


MetricIncidenceMeansPlot <- function(df, 
                                     metric = input_metric, 
                                     incidence = "reported_incidence_per_1k"){
  summary <- df %>%
    dplyr::mutate(rounded = round(.data[[metric]], digits=1)) %>%
    dplyr::group_by(rounded) %>%
    dplyr::summarise(
      total_simulations = n_distinct(sim_id),
      mean = mean(.data[[incidence]]),
      sd = sd(.data[[incidence]]),
      min = min(.data[[incidence]]),
      max = max(.data[[incidence]]),
      .groups = 'drop')
  
  p_inc <- summary %>%
    ggplot(aes(x=rounded, y=mean, 
               size=total_simulations)) +
    geom_errorbar(aes(ymin = mean - sd,
                      ymax = mean + sd),
                  size = 0.5) +
    geom_point(aes(size=total_simulations)) +
    labs(y = "Yearly reported incidence per 1k",
         x = names(y_labs)[y_labs == metric],
         size="Total\nsimulations") +
    theme(
      legend.position = c(0.2, 0.85),
      legend.background = element_rect(fill = "#F5F3ED")
    )
  
  return(p_inc)
}

################################################################################
# Statistical results plotting functions
################################################################################
my_breaks <- c(1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 5e-1)
IncidenceBinSignificancePlots <- function(df, 
                                          x_facet = ".y.",
                                          y_facet = "pair_var"){  
  p <- df %>%
    ggplot(aes(x=group1_bin_order, y=group2_bin_order, fill=p.adj)) +
    geom_tile() +
    facet_grid(as.formula(paste(x_facet, "~", y_facet)),
               labeller = labeller(`.y.` = metric_labels)) +
    scale_fill_gradientn(
      name = "Adjusted p-value",
      colours = rev(c("#67322EFF", "#99610AFF", "#C38F16FF", "#6E948CFF", "#2C6B67FF", "#175449FF", "#122C43FF")),
      values = scales::rescale(my_breaks),  
      breaks = my_breaks,
      limits = c(0, 1),
      labels = scales::number_format(accuracy = 0.01)
    ) +
    guides(fill = guide_colorbar(
      ticks     = TRUE,
      barheight = unit(6, "cm"),
      barwidth  = unit(0.5, "cm")
    )) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_text(data = 
                dplyr::filter(df, p.adj < 0.05),
              aes(label=p.adj.signif), size=2, color="white") + 
    labs(x="", y="",
         title="Statistically significant genetic metric changes between incidence bins",
         subtitle="BH adjusted p-value (**** < 0.00005, *** < 0.0005, ** < 0.005, * < 0.05)")
  return(p)
}