################################################################################
# Plotting colors and shapes
################################################################################
season_peak_colors <- setNames(
  c("#D3AA6AFF", "#95CAA6FF", "#008D98FF"),
  c("Other", "Peak dry season", "Peak wet season")
)

inf_count_colors =  setNames(
  c("#607D8BFF", "#37474FFF"),
  c("Parasite genomes", "Infected individuals"))

year_shapes <- setNames(c(1, 15, 16, 17, 18), 
                        seq(-1, 3))

sampling_scheme_colors <- setNames(
  c("#552000FF", "#8A4D00FF", "#C17D17FF", 
    "#6699FFFF", "#003399FF",
    "#FF6600FF", "#CC0000FF",
    #"#6c1446"
    "#d83397"),
  c("All - Yearly", "Sample - Proportional", "Sample - Even", 
    "Sample - Seasonal (Wet)", "Sample - Peak Seasonal (Wet)",  
    "Sample - Seasonal (Dry)", "Sample - Peak Seasonal (Dry)",
    "Both schemes combined")
)


################################################################################
# Seasonal window objects
################################################################################
highlight_periods <- data.frame(
  xmin = c(8, 20, 32, 44, 56, 68, 80),
  xmax = c(13, 25, 37, 49, 61, 73, 85),
  ymin = -Inf,  # span full y-range
  ymax = Inf
)

################################################################################
# Plotting functions
################################################################################
MonthlyTrendsPlot <- function(df, 
                              x_variable = month_start,
                              y_variable, 
                              bin = "year_who_bin", 
                              sd_ribbons=TRUE,
                              intervention_line = intervention_intercept){
  summary_df <- df %>%
    dplyr::group_by(.data[[bin]], .data[[month_start]], itn_on, importation_rate) %>% 
    summarize(
      count = n(),
      mean = mean(.data[[y_variable]], na.rm = TRUE), 
      sd = sd(.data[[y_variable]], na.rm = TRUE), .groups = 'drop') 
  
  
  p <- summary_df %>%
    dplyr::filter(.data[[month_start]] %in% (min(shifted_full_years):max(shifted_full_years))) %>%
    ggplot(aes(x=.data[[month_start]], 
               y=mean, 
               color=.data[[bin]], 
               linetype=importation_rate, 
               group = paste(.data[[bin]], importation_rate)
    )) +
    geom_rect(data = highlight_subset,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE,
              fill = "blue", alpha = 0.2) +
    {
      if(isTRUE(sd_ribbons)){
        list(
          geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd, 
                          fill = .data[[bin]]), alpha = 0.2),
          geom_line(),
          geom_point()
        )
      } else{
        list(geom_line(aes(x=.data[[x_variable]], y=.data[[y_variable]], group=sim_id), 
                       data = df, alpha=0.25),
             geom_line(),
             geom_point()
        )
      }
    } +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    geom_vline(xintercept = intervention_line, linewidth=1.5) +
    labs(x="Month", 
         y="",
         title=names(y_labs)[y_labs == y_variable],
         color = "Yearly incidence", fill = "Yearly incidence",
         #linetype="Insecticide treated\nnet intervention") +
         linetype="Importation rate") +
    scale_x_continuous(breaks = seq(min(shifted_full_years), max(shifted_full_years), by = 6))
  return(p)
}

MonthlyMetricSummary <- function(df_long, 
                                 month_group = month_start){
  
  tmp <- df_long %>%
    dplyr::group_by(.data[[month_group]], genetic_metric) %>%
    summarise(
      mean_metric = mean(value, na.rm = TRUE),
      sd_metric   = sd  (value, na.rm = TRUE),
      .groups = "drop"
    ) 
  return(tmp)
}


MonthlyCombinedMetrics <- function(df, 
                                   x_variable = month_start,
                                   month_group = month_start,
                                   bin = "start_who_bin"){
  
  base_p <- function(df_long, 
                     x_variable = month_start,
                     y_variable = "value", 
                     color_variable = "genetic_metric", 
                     sd_ribbons=TRUE){
    
    p_df <- df_long
    if(isTRUE(sd_ribbons)){
      p_df <- MonthlyMetricSummary(df_long) 
      y_variable <- "mean_metric"
    }
    
    p <- p_df %>%
      ggplot(aes(x=.data[[x_variable]], 
                 y=.data[[y_variable]], 
                 color=.data[[color_variable]], 
                 fill=.data[[color_variable]],
      )) +
      geom_rect(data = highlight_subset,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                fill = "blue", alpha = 0.2) +  
      {
        if(isTRUE(sd_ribbons)){
          list(
            geom_ribbon(aes(ymin = mean_metric - sd_metric, 
                            ymax = mean_metric + sd_metric, 
                            fill = .data[[color_variable]]), alpha = 0.2),
            geom_line(aes(group=.data[[color_variable]])),
            geom_point()
          )
        } else{
          # TODO: Fix
          list(geom_line(aes(x=.data[[x_variable]], 
                             y=.data[[y_variable]], 
                             group=c(sim_id, .data[[color_variable]])), 
                         alpha=0.25),
               geom_point()
          )
        }
      } +
      geom_vline(xintercept = intervention_intercept, linewidth=1.5) +
      scale_x_continuous(breaks = seq(min(shifted_full_years), max(shifted_full_years), by = 6)) +
      labs(x="Month", fill="", color="") 
    
    return(p)
  }
  
  # Remove incomplete years at start and end of simulation  
  df <- dplyr::filter(df, .data[[x_variable]] %in% (min(shifted_full_years):max(shifted_full_years)))
  # Format to long 
  df_long <- df %>%
    dplyr::select(sim_id, all_of(c(month_group, bin)), n_infections, genome_ids_total_count, all_of(y_vars)) %>%
    tidyr::pivot_longer(
      cols = c("n_infections", "genome_ids_total_count", y_vars),
      names_to = "genetic_metric",
      values_to = "value",
    ) 
  
  incidence_p <- base_p(df_long = select(df, all_of(month_start), reported_incidence_per_1k) %>% 
                          dplyr::mutate(genetic_metric="incidence") %>%
                          dplyr::rename(value=reported_incidence_per_1k)) +
    scale_fill_manual(values="#111111") +
    scale_color_manual(values="#111111") +
    ylim(0, 110) + #max(df$reported_incidence_per_1k)) +
    labs(title="Monthly reported incidence per 1k", y="") +
    theme(legend.position = "none")
  
  count_columns <- c("n_infections", "genome_ids_total_count")
  count_p <- base_p(dplyr::filter(df_long, genetic_metric %in% count_columns) ) +
    scale_color_manual(values = c("#607D8BFF", "#37474FFF"), labels = c("Parasite genomes", "Infected individuals")) +
    scale_fill_manual(values = c("#607D8BFF", "#37474FFF"), labels = c("Parasite genomes", "Infected individuals")) +
    ylim(c(0, max_infections)) +
    labs(title="Infected individuals and genomes", y="", fill="", color="") +
    theme(legend.position = "right")
  
  metric_p <- base_p(dplyr::filter(df_long, !genetic_metric %in% c(count_columns, "superinfection_prop", "effective_coi_mean", "true_coi_mean"))) +
    scale_color_manual(labels=metric_labels, values=metric_colors) +
    scale_fill_manual(labels=metric_labels, values=metric_colors) +
    ylim(c(0,1)) +
    labs(title="Genetic metrics",
         y="Metric proportion", 
         fill="", color="") +
    theme(legend.position = "right")
  
  joint_p <- ggpubr::ggarrange(incidence_p + theme(axis.title.x=element_blank()),
                               count_p + theme(axis.title.x=element_blank()), 
                               metric_p, 
                               nrow=3, align="hv", heights=c(1,1,2))
  return(joint_p)
}


CombinedSamplingSchemeMeansPlot <- function(df){
  pd <- position_dodge(width = .45)                # tweak width to separate points
  p <- df %>%
    arrange(sampling_scheme, year_who_bin, metric) %>%
    ggplot(aes(x = mean, y = year_who_bin, color=sampling_scheme, group=sampling_scheme)) +
    geom_errorbarh(aes(xmin = mean - sd,
                       xmax = mean + sd),
                   height   = .15,
                   linewidth = .35,
                   position = pd) +
    geom_point(size = 3,
               position = pd) +
    xlim(c(0,1)) +
    facet_grid(year~metric,                            
               labeller = as_labeller(metric_labels)) +
    scale_color_manual(values = sampling_scheme_colors, breaks = names(sampling_scheme_colors)) +
    labs(x="Metric proportion\n(Mean -/+ 1 SD)",
         y="WHO transmission\nclassification",
         color = "Sampling scheme") 
  return(p)
}  


# Dot and whsiker plot of means for different temporal sampling schemes
sampling_scheme_ord <- gsub(" - ", "\n", names(sampling_scheme_colors))
IndividualSamplingSchemeMeansPlot <- function(df,
                                              y_variable = "sampling_scheme",
                                              color_variable = "metric",
                                              group_variable = "year"){
  
  df <- df %>%
    dplyr::mutate(sampling_scheme = factor(sampling_scheme, levels = rev(sampling_scheme_ord)),
                  metric = factor(metric, levels = names(metric_labels)))

  
  vline.data <- df %>%
    dplyr::filter(sampling_scheme == sampling_scheme_ord[grepl("Yearly", sampling_scheme_ord)] & year == -1)
  
  pd <- position_dodge2(width = 0.75, 
                        preserve = "total",
                        reverse  = TRUE)
  p <- df %>%
    ggplot(aes(x = mean, 
               y=.data[[y_variable]], 
               color=.data[[color_variable]])) +
    {
      if (length(unique(df$year)) > 1) {
        geom_vline(
          data = vline.data,
          aes(xintercept = mean, 
              linetype = factor(.data[[group_variable]])),
          colour     = "#666666",
          size       = 0.25, 
          show.legend = FALSE
        )
      } else {
        NULL
      }
    } +
    geom_errorbarh(aes(xmin = mean - sd,
                       xmax = mean + sd,
                       group   = factor(.data[[group_variable]])),
                   height   = 0.75, 
                   position = pd) +
    geom_point(aes(shape   = factor(.data[[group_variable]]),
                   group   = factor(.data[[group_variable]])),
               size=1.75,
               position = pd) +
    facet_grid(.~metric,                            
               labeller = as_labeller(metric_labels), 
               scales="free_x") +
    # ggh4x::facetted_pos_scales(
    #   x = list(
    #     metric %in% y_vars[grepl("prop", y_vars)] ~ scale_x_continuous(
    #       limits = c(0, 1),
    #       breaks = seq(0, 1, 0.5),
    #       labels = scales::label_number(drop0trailing=TRUE)),
    #     !metric %in% y_vars[grepl("prop", y_vars)] ~ scale_x_continuous(
    #       limits=c(1,2),
    #       breaks = seq(1,2,0.5))
    #   )
    # ) +
    ggh4x::facetted_pos_scales(
      x = list(
        metric == "poly_coi_prop"          ~ scale_x_continuous(limits = c(0.1, 0.5), breaks = seq(0, 0.4, 0.2)),
        metric == "genome_ids_unique_prop" ~ scale_x_continuous(limits = c(0.2, 1), breaks = seq(0.2, 1,0.2)),
        metric == "cotransmission_prop"    ~ scale_x_continuous(limits= c(0.2, 0.8), breaks = seq(0.2,0.8,0.2)),
        metric == "effective_coi_mean"     ~ scale_x_continuous(limits = c(1,2),
                                                                breaks = seq(1,2,0.5),
                                                                labels = scales::label_number(drop0trailing=TRUE))
        )
      ) +
    {
      if (color_variable == "metric") {
        list(scale_color_manual(values = metric_colors),
             guides(color=FALSE),
             theme(legend.position = "top"))
      } else {
        NULL
      }
    } +
    scale_shape_manual(values=year_shapes) +
    labs(x="Mean -/+ 1 SD",
         #y="Sampling scheme", 
         y="",
         shape="Year", 
         linetype="Year") +
    guides(linetype=FALSE) +
    theme(panel.spacing.x = unit(0.5, "cm"))
  #theme_classic() 
  return(p)
}

