################################################################################
# Purpose: Summary statistic plotting functions from EMOD FPG infection report
# Author: Jessica Ribado
# Date: 04/2025
################################################################################

################################################################################
# set-up 
################################################################################
# load libraries
needed_packages <- c('data.table', 'dplyr', 'tidyr', 'purrr',
                     'ggplot2','ggpubr', 'ggh4x', 'gtools', 'ggtext', 'rstatix')
for(p in needed_packages){
  if(!p %in% installed.packages()[,1]){
    install.packages(p)
    library(p, character.only = TRUE)
  }
  library(p, character.only = TRUE)
}

# load helper functions
source("temporal_plotting.R")
source("incidence_statistics.R")
source("incidence_bin_plotting.R")


# Set directories
dir <- '/mnt/data/malaria/synthetic_genomes/jessica_projects/FPG_ObsModelTesting'
summary_dir <- paste(dir, "infectionFPGReport_sympomaticsOnly_timeCheck", sep="/")
plot_dir <- paste(dir, "plot_testFPG_symptomaticsOnly", sep="/")
if (!dir.exists(plot_dir)){ dir.create(plot_dir, recursive = TRUE) }

# create mapping file to run with observational_extracted.py to give matched sim_id names
# sim_epi_mapping  <- data.table::fread(paste(dir, "sim_data_epi_new_itns_6yr_very_high_biting.csv", sep="/")) %>%
#   dplyr::select(sim_id, outpath) %>% unique() %>%
#   dplyr::mutate(outpath = paste0("\"", gsub("/mnt/idm2/home/", "/mnt/calculon2/", outpath), "/output\"")) %>%
#   dplyr::rename(output_name=sim_id, input_dir=outpath)
# write.table(sim_epi_mapping, file = paste(dir, "sim_mapping_new_itns_6yr_highBiting.csv", sep="/"), quote = F, row.names = F, sep=",")

shift_intervention_time = TRUE
intervention_month_shift = 29
month_start <- ifelse(isTRUE(shift_intervention_time), "intervention_month", "continuous_month")
intervention_intercept <- ifelse(isTRUE(shift_intervention_time), 0, intervention_month_shift)

# shift wet season shading for any intervention shifts 
if(isTRUE(shift_intervention_time)){
  highlight_periods <- highlight_periods - intervention_month_shift
}


################################################################################
# plotting options and functions
################################################################################
# Set plotting options
options(scipen=10000)
theme_set(theme_bw()+
            theme(
              panel.background = element_rect(fill='transparent'), #transparent panel bg
              plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
              panel.grid.major = element_blank(), #remove major gridlines
              panel.grid.minor = element_blank(), #remove minor gridlines
              legend.background = element_rect(fill='transparent', colour = NA), #transparent legend bg
              legend.box.background = element_rect(fill='transparent', colour = NA),
              legend.key = element_rect(fill=NA, colour = NA) #transparent legend panel
            ))


metric_labels <- c(
  true_poly_coi_prop     = "True polygenomic\nproportion",
  effective_poly_coi_prop = "Polygenomic\nproportion",
  variant_poly_coi_prop = "Polygenomic\nproportion",
  all_genomes_unique_prop = "All unique genome\nproportion",
  mono_genomes_unique_prop = "Unique genome\nproportion",
  cotransmission_prop    = "Cotransmission\nproportion",
  effective_coi_mean    = "Mean complexity\nof infection",
  true_coi_mean         = "Mean true complexity of infection",
  rh_poly_inferred_mean = "Mean RH",
  ibs_mean = "Mean identity by state",
  ibd_mean = "Mean identity by descent"
)

metric_colors <- c(
  true_poly_coi_prop        = "mediumpurple4",
  effective_poly_coi_prop = "mediumpurple1",
  variant_poly_coi_prop   = "#7C4B73FF",
  all_genomes_unique_prop  = "#81A88DFF",
  mono_genomes_unique_prop = "#02401BFF",
  cotransmission_prop    = "#D8B70AFF",
  effective_coi_mean     = "#972D15FF",
  true_coi_mean          = "#E68E54FF",
  rh_poly_inferred_mean = "#06A2BCFF",
  ibs_mean = "#5A97C1FF", 
  ibd_mean = "#0A2E57FF"
)


################################################################################
metrics <- c("reported_incidence_per_1k",
             "true_incidence_per_1k",
             "eir",
             "pfpr")

# Function to summarize metrics for intervention shifted summary statistics of summary reports
Monthly2Yearly <- function(df){
  df <- dplyr::select(df, -month)
  group_cols <- setdiff(names(df), metrics)
  
  df_summary <- df %>%
    # 1) group by every column except the four metrics
    group_by(across(all_of(group_cols))) %>%
    # 2) sum each metric, dropping the grouping afterward
    summarise(
      across(all_of(metrics), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  return(df_summary)
}


TimeScaleSummaryStats <- function(df = year_summaries, 
                                  time_scale = "group_year",
                                  bin = "start_who_bin") {
  
  multi_sim_year_summary <- df %>%
    dplyr::group_by(.data[[bin]], sampling_scheme, .data[[time_scale]]) %>%
    summarise(
      total_count = n(),               # <-- row count per group
      across(
        all_of(y_vars),
        .fns = list(
          mean   = ~mean(as.numeric(.x),   na.rm = TRUE),
          median = ~median(as.numeric(.x), na.rm = TRUE),
          sd     = ~sd(.x,     na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
  
  multisim_summary_long <- multi_sim_year_summary %>% 
    pivot_longer(
      # keep the three identifier columns as-is
      cols = -c(eval(bin), sampling_scheme, eval(time_scale), total_count),
      names_to      = c("metric", ".value"),
      names_pattern = "(.+?)_(sum|mean|median|sd)$"
    )
  
  return(multisim_summary_long)    
}

################################################################################
# load data
################################################################################

################################################################################
# Empirical data
################################################################################
senegal_data <- data.table::fread(paste(dir, "../2504_GRSweep/combined_model_predictions_inc24_long.csv", sep="/"), sep=",") %>%
  dplyr::mutate(sampling_scheme = "Senegal data")

################################################################################
# Simulation epi data (provided by analyzer not yet defined)
################################################################################
sim_epi_yearly_path <- ""
sim_epi_monthly_path <- paste(dir, "sim_data_monthly.csv", sep="/")

sim_epi_monthly <- data.table::fread(sim_epi_monthly_path) %>%
  dplyr::rename(simulation_year = year) %>%
  dplyr::mutate(continuous_month = ifelse(simulation_year == 0, month - 1, (simulation_year*12) + month - 1),
                group_year = simulation_year)

if(isTRUE(shift_intervention_time)){
  sim_epi_monthly <- sim_epi_monthly %>% 
    dplyr::mutate(intervention_year = (continuous_month - intervention_month_shift) %/% 12,
                  intervention_month = continuous_month - intervention_month_shift,
                  group_year = intervention_year)
  
  emod_columns <- setdiff(names(sim_epi_monthly), c("simulation_year", "month", "continuous_month", "intervention_month", metrics))
  sim_epi_yearly <- sim_epi_monthly %>%
    dplyr::group_by(across(all_of(c(emod_columns)))) %>%
    dplyr::summarise(
      reported_incidence_per_1k = sum(reported_incidence_per_1k),
      true_incidence_per_1k = sum(true_incidence_per_1k),
      eir = mean(eir),
      pfpr = mean(pfpr), .groups = 'drop'
    )
} else{
  sim_epi_yearly  <- data.table::fread(sim_epi_yearly_path) 
  if("month" %in% names(sim_epi_yearly)) {sim_epi_yearly <- Monthly2Yearly(sim_epi_monthly)}
}

sim_epi_yearly <- sim_epi_yearly %>%
  dplyr::filter(habitat_scale > 6.8) %>%
  dplyr::mutate(output_name = basename(outpath), 
                year_who_bin = case_when(
                  reported_incidence_per_1k < 100 ~ "Very low",
                  reported_incidence_per_1k < 250 ~ "Low",
                  reported_incidence_per_1k < 450 ~ "Moderate",
                  reported_incidence_per_1k > 450 ~ "High"
                ),
                year_who_bin = factor(year_who_bin, levels = c("Very low", "Low", "Moderate", "High")))
sim_epi_yearly$year_incidence_bin <- cut(sim_epi_yearly$reported_incidence_per_1k, 
                                         seq(plyr::round_any(min(sim_epi_yearly$reported_incidence_per_1k), 50, f=floor), 
                                             plyr::round_any(max(sim_epi_yearly$reported_incidence_per_1k), 50, f=ceiling), 50),
                                         include.lowest = TRUE)
sim_start_bins <- dplyr::filter(sim_epi_yearly, group_year==-1) %>% select(sim_id, year_who_bin) %>% unique() %>% rename(start_who_bin=year_who_bin)

sim_epi_yearly <- left_join(sim_epi_yearly, sim_start_bins)
sim_epi_monthly <- inner_join(sim_epi_monthly, sim_start_bins)


################################################################################
# Simulatation genetic data
################################################################################
genetic_summary_files <- list.files(path = summary_dir, pattern = ".*_FPG_ModelSummaries.csv", recursive = T)

# Option 1: Force consistent column types during reading
genetic_stats_all <- dplyr::bind_rows(
  lapply(setNames(genetic_summary_files, genetic_summary_files), function(i) {
    df <- data.table::fread(paste(summary_dir, i, sep="/"))
    
    # Convert ALL columns to character first to eliminate type conflicts
    df[] <- lapply(df, as.character)
    
    return(df)
  }),
  .id = "sim_id"
) %>%
  dplyr::mutate(sim_id = gsub(".*\\/|_FPG_ModelSummaries.csv", "", sim_id))


year_summaries <- dplyr::filter(genetic_stats_all, comparison_type == "group_year" ) %>%
  dplyr::mutate(sampling_scheme = "Sample - Proportional",
                year_group = as.numeric(year_group)) %>%
  dplyr::rename("group_year" = "year_group")
season_summary <- dplyr::filter(genetic_stats_all, comparison_type == "season_bins") %>%
  dplyr::mutate(simulation_year = as.numeric(sub(".*:\\s*(-?\\d+)-.*", "\\1", year_group)),
                season = stringr::str_extract(year_group, "^\\w+"),
                # sampling_scheme = paste0(gsub(".*seasonal|_[0-9].*", "", sampling_scheme), " season (", season, ")")) %>%
                sampling_scheme = paste0("Sample - Seasonal (", season, ")")) %>%
  dplyr::mutate(group_year = case_when(
                  season == "Wet" ~ simulation_year - 2,
                  season == "Dry" & simulation_year == 1 ~ -1,
                  season == "Dry" & simulation_year > 2  ~ simulation_year - 3,
                  .default = NA_real_
                ))
# peak_summary <- dplyr::filter(genetic_stats_all, grepl("Wet|Dry", peak_season)) %>%
#   dplyr::mutate(year = as.numeric(sub(".*:\\s*(-?\\d+)-.*", "\\1", peak_season)),
#                 season = stringr::str_extract(peak_season, "^\\w+"),
#                 sampling_scheme = paste0(sampling_scheme, " (", season, ")"))
year_summaries <- dplyr::bind_rows(year_summaries, 
                                   season_summary, 
                                   #peak_summary
) 
year_summaries <-  dplyr::inner_join(year_summaries, sim_epi_yearly) %>%
  dplyr::mutate(sampling_scheme = factor(sampling_scheme, levels = names(sampling_scheme_colors)))


month_summaries <- dplyr::filter(genetic_stats_all, !is.na(continuous_month)) %>% 
  dplyr::select(-year, -season) %>%
  #dplyr::select(-year, -intervention_month, -month) %>%
  dplyr::inner_join(., sim_epi_monthly) %>%
  dplyr::mutate(
    peak_group = case_when(
      between(month, 10, 12) ~ "Peak wet season",
      between(month, 3, 6) ~ "Peak dry season",
      TRUE ~ "Other"),
    intervention_month = continuous_month - intervention_month_shift,
    year = intervention_month %/% 12
  ) 
month_summaries <- dplyr::inner_join(month_summaries, 
                                     select(sim_epi_yearly, sim_id, year, reported_incidence_per_1k, year_who_bin) %>% 
                                       rename("year_incidence" = "reported_incidence_per_1k"))
max_infections <- max(month_summaries$genome_ids_total_count)
max_incidence <- max(month_summaries$reported_incidence_per_1k)

# set some time limits
#time_vector <- sort(unique(month[[month_start]]))
#shifted_full_years <- time_vector[time_vector%%12 == 0]
shifted_full_years=unique(year_summaries$group_year)
#shifted_full_years <- shifted_full_years[shifted_full_years > -13]
shifted_full_years <- shifted_full_years[shifted_full_years > -2]
highlight_subset <- dplyr::filter(highlight_periods, xmin > min(shifted_full_years) & xmin < max(shifted_full_years)) %>%
  dplyr::mutate(xmax = ifelse(xmax > max(shifted_full_years), max(shifted_full_years), xmax))

sampling_subset <- names(sampling_scheme_colors)[c(2,4,6)]
year_summaries <- dplyr::filter(year_summaries, between(group_year, -1, 2) & sampling_scheme %in% sampling_subset)
month_summaries <- month_summaries %>% dplyr::filter(.[[month_start]] >= min(shifted_full_years) & .[[month_start]] <= max(shifted_full_years)-1)

#write.table(month_summaries, paste(dir, "geneticMetrics_Monthly_itns_6yr_PeakDaysFilter.csv", sep="/"), 
#            sep=",", quote=F, row.names = F)
#write.table(year_summaries, paste(dir, "geneticMetrics_Yearly_itns_6yr_PeakDaysFilter.csv", sep="/"), 
#            sep=",", quote=F, row.names = F)


################################################################################
# Set variables for included summary statistics
y_vars <- c(names(genetic_stats_all)[grepl("prop|_mean", names(genetic_stats_all))])
y_labs <- gsub("\n" ," ", metric_labels)
y_names <- setNames(y_vars, c("AllPoly", "EffPoly", "AllUnique", "MonoUnique", 
                              "Cotx", "EffCOI", "TrueCOI", "RH", "IBS"))  

# Subset genetic metrics of interest
# y_vars <- y_vars[!y_vars %in% c("true_coi_mean", "superinfection_prop")] 
y_vars <- c("effective_poly_coi_prop", "mono_genomes_unique_prop", "cotransmission_prop", "rh_poly_inferred_mean", "effective_coi_mean")


################################################################################
# Senegal data 
################################################################################
district_colors <- c(
  "Bakel" =	"#A6CEE3",
  "Diourbel" = "#B2DF8A",
  "Kedougou" =	"#E31A1C",
  "Kolda"	= "#006600",
  "Makacoulibantang"	= "#993300",
  "Pikine" = "#FFCCFF",
  "Richard Toll"	= "#FF66FF",
  "Sedhiou" =	"#FDBF6F",
  "Tambacounda"	= "#FF7F00",
  "Thies" = "#0000FF",
  "Touba" =	"#CAB2D6",
  "Velingara"	= "#6A3D9A"
)  

axis_size = 12
title_size= 14
sen_p <- senegal_data %>% 
  dplyr::filter(!Metric %in% c("RH", "Monclonality")) %>%
  dplyr::mutate(Metric=gsub(" prop", "\nprop", Metric),
                Metric = gsub("-", "", Metric),
                Metric = gsub("COI", "complexity\nof infection", Metric)) %>%
  dplyr::mutate(Metric = factor(Metric, levels = c("Polygenomic\nproportion", "Unique genome\nproportion*", "Cotransmission\nproportion", "Mean complexity\nof infection"))) %>%
  ggplot(aes(x=Incidence, y=Value, fill=District, size=n))+
  geom_point(shape=21, color="darkgray", alpha=0.8)+
  scale_fill_manual(values = district_colors) +
  labs(x = "Incidence", y = "Metric value", 
       size="Barcoded\nsamples",
       caption = "*calculated using only monogenomic samples") +
  facet_wrap(~ Metric, ncol = 2, scales = "free_y") +
  scale_x_log10() +
  theme(plot.caption = element_text(size = 16))+
  theme_bw() +
  theme(axis.ticks.length.x = unit(0.4, "cm"), axis.ticks.length.y = unit(0.2, "cm")) +
  theme(plot.caption = element_text(size = axis_size),
        axis.text.x = element_text(size = axis_size),  # bigger x-axis labels
        axis.text.y = element_text(size = axis_size),   # bigger y-axis labels
        axis.title.x = element_text(size = title_size),
        axis.title.y = element_text(size = title_size),
        strip.text = element_text(size = axis_size),   # Bigger facet labels
        #legend.text = element_text(size = 20),
        #legend.title = element_text(size = 20)
  )

ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_SenegalMetrics.png"), 
       plot = sen_p, path = plot_dir,
       width = 6, height = 5.25, units = c("in"), dpi = 300)



################################################################################
# Yearly plots
################################################################################
# Best fit line for all years in the simulation
year_summaries <- year_summaries %>% dplyr::mutate(
  sampling_scheme = gsub("Sample - ", "", sampling_scheme), 
  sampling_scheme = factor(sampling_scheme, levels = names(sampling_colors_2)))

all_year_trends <- lapply(y_vars, function(i, x_scale = "log"){ 
  #p <- ProportionBasePlot(year_summaries, y_var = i, summarize_trend = TRUE) +
  p <- PlotBestFit(year_summaries, y = i) 
  p <- p$plot + 
    {
      if(x_scale == "log"){
        scale_x_log10() 
      }
    } + 
    labs(x="Reported incidence per 1k", 
         y=y_labs[names(y_labs) == i],
         color = "Sampling scheme",
         fill = "Sampling scheme") 
  
  ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_Yearly_", names(y_names)[y_names == i], "_log.png"), 
         plot = p, path = plot_dir,
         width = 8, height = 5, units = c("in"), dpi = 300)
  return(p)
})


########
# Scatter plot of yearly v monthly incidence, colored by peak seasons 
# month_summaries <- dplyr::mutate(month_summaries, year=as.character(year))
# all_month_trends <- lapply(y_vars, function(i){
#   p <- ProportionBasePlot(month_summaries, y_var = i, 
#                      color_variable = "season_group", 
#                      # color_variable = "year_incidence"
#                      shape_variable = "year"
#                      ) +
#     labs(x="Monthly reported incidence per 1k", 
#          y=names(y_labs)[y_labs == i],
#          color = "Sampling season", fill = "Sampling season",
#          # color = "Mean yearly\nincidence per 1k", fill = "Mean yearly\nincidence per 1k",
#          shape="Year")
#   
#   ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_MonthlyIncidence_", names(y_names)[y_names == i], ".png"), 
#          plot = p, path = plot_dir,
#          width = 8, height = 5, units = c("in"), dpi = 300)
#   return(p)
# })


########
# Best fit plots for transmission ranges 
senegal_labels <- c(
  effective_poly_coi_prop  = "Polygenomic proportion",
  mono_genomes_unique_prop = "Unique genome proportion*",
  cotransmission_prop    = "Co-transmission proportion",
  rh_poly_inferred_mean = "RH",
  effective_coi_mean    = "Mean COI",
  true_coi_mean         = "Mean COI"
)


for(model_year in c(unique(year_summaries$group_year), "All")){
  
  metric_plots_single <- lapply(y_vars, function(j, x_scale = "log", add_senegal_points=TRUE){
    if(model_year == "All"){
      tmp_df <- year_summaries
    } else{
      tmp_df <- dplyr::filter(year_summaries, group_year==model_year)
    }
    
    tmp_df <- tmp_df %>%
      dplyr::mutate(
        sampling_scheme = gsub("Sample - ", "", sampling_scheme), 
        sampling_scheme = factor(sampling_scheme, levels = names(sampling_colors_2))) 
    regression_run <- PlotBestFit(tmp_df, y = j)
    
    p <- regression_run$plot +
      {
        if(isTRUE(add_senegal_points)){
          geom_point(data = dplyr::filter(senegal_data, 
                                          between(Incidence, min(tmp_df$reported_incidence_per_1k),
                                                  max(tmp_df$reported_incidence_per_1k)) &
                                            Metric == senegal_labels[names(senegal_labels) == j]), 
                     aes(x=Incidence, y=Value), 
                     inherit.aes = FALSE,            # drop the original color/group mapping
                     color       = "black",          # or pick a new color aesthetic
                     shape       = 17,               # a different shape so it stands out
                     size        = 2) 
        }    
      } +
      labs(x="Yearly mean reported incidence per 1k", 
           title=y_labs[names(y_labs) == j],
           color = "Sampling scheme",
           fill = "Sampling scheme") +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) + 
      {
        if(x_scale == "log"){
          scale_x_continuous(breaks = c(100, 250, 450), transform = "log10") 
        }    
      } 
    
    return(p)
  })
  
  p <- ggpubr::ggarrange(
    plotlist=list(metric_plots_single[1][[1]], 
                  metric_plots_single[2][[1]],
                  metric_plots_single[3][[1]],
                  metric_plots_single[4][[1]],
                  metric_plots_single[5][[1]]
    ),
    ncol=3, nrow=2, align="hv", common.legend = T, legend = "top")
  p <- annotate_figure(p,
                       bottom = text_grob("Yearly mean reported incidence per 1k",
                                          size = 12))
  p
  ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_YearlyAllMetrics_AllSchemes_yr", model_year, "_logWithSamples.png"), 
         plot = p, path = plot_dir,
         width = 7.75, height = 5.5, units = c("in"), dpi = 300)
}



##########
# Smoothed line plots of genetic metrics per year
# color_variable = "Intervention"
# line_variable="year"
# single_sim <- lapply(y_vars, function(j, x_scale="log"){
#   p <- dplyr::filter(year_summaries, sampling_scheme == "Sample - Seasonal (Wet)" ) %>%
#     dplyr::mutate(Intervention = ifelse(itn_on == "TRUE", "Insecticide treated nets", "None"),
#                   importation_rate = stringr::str_to_title(importation_rate)) %>%
#     ggplot(aes(x=reported_incidence_per_1k, 
#                y = .data[[j]],
#                color=get(color_variable), 
#                fill=get(color_variable),
#                linetype=as.character(get(line_variable))
#     )) +
#     geom_smooth(method = "lm", formula = y ~ poly(x, 3, raw = TRUE), 
#                 se = 0.68, alpha=0.5)  +
#     geom_point(alpha=0.1) +
#     scale_color_manual(values = c("#D53E4FFF", "#FDAE61FF")) +
#     scale_fill_manual(values = c("#D53E4FFF", "#FDAE61FF")) +
#     labs(x="Yearly mean reported incidence per 1k", 
#          title=names(y_labs)[y_labs == j],
#          # color = "Importation rate",
#          # fill = "Importation rate",
#          # linetype="Intervention"
#          color = "Intervention",
#          fill = "Intervention",
#          linetype="Year") +
#     {
#       if(x_scale == "log"){
#         scale_x_log10() 
#       }    
#     } +
#     theme(axis.title.y=element_blank()) 
#     
#   
#   ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_YearAllWetSeasonPerYear_", names(y_names)[y_names == j], "_log.png"), 
#          plot = p, path = plot_dir,
#          width = 6, height = 5, units = c("in"), dpi = 300)
#   
#   return(p)
# })


###########
# Scatter plot of monthly reported incidence and genetic metrics
# all_month_trends <- lapply(y_vars, function(i){
#   tmp <- dplyr::filter(month_summaries, year==1) %>%
#     dplyr::mutate(Intervention = ifelse(itn_on == "TRUE", "Insecticide treated nets", "None"),
#                   importation_rate = stringr::str_to_title(importation_rate)) 
#   p <- ProportionBasePlot(tmp, y_var = i, 
#                           # color_variable = "season_group", 
#                           color_variable = "year_incidence",
#                           y_facet = "Intervention"
#   ) +
#     labs(x="Monthly reported incidence per 1k", 
#          y=names(y_labs)[y_labs == i],
#          # color = "Sampling season", fill = "Sampling season",
#          color = "Mean yearly\nincidence per 1k", fill = "Mean yearly\nincidence per 1k") 
#   
#   ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_YearAllMonthlyByYearlyIncidence_", names(y_names)[y_names == i], ".png"), 
#          plot = p, path = plot_dir,
#          width = 8, height = 5, units = c("in"), dpi = 300)
#   return(p)
# })


###########
# Line plot of changes in summary statistic on a monthly scale, with incidence to match
lapply(y_vars, function(i){
  month_incidence_p <- MonthlyTrendsPlot(month_summaries, 
                                         y_variable = "reported_incidence_per_1k", 
                                         bin="start_who_bin", 
                                         sd_ribbons = FALSE) + 
    labs(title = "Reported incidence per 1k",
         color = "Starting incidence\nbin", 
         fill = "Starting incidence\nbin") +
    geom_vline(xintercept=intervention_intercept, size=1.5)
  
  month_prevlance_p <- MonthlyTrendsPlot(month_summaries, 
                                         y_variable = "pfpr", 
                                         bin="start_who_bin", 
                                         sd_ribbons = FALSE) + 
    labs(title = "PfPR",
         color = "Starting incidence\nbin", 
         fill = "Starting incidence\nbin") +
    geom_vline(xintercept=intervention_intercept, size=1.5)
  
  month_metric_p <- MonthlyTrendsPlot(month_summaries, y_variable = i, 
                                      bin="start_who_bin", sd_ribbons = FALSE) + 
    labs(
      title=y_labs[names(y_labs) == i],
      color="Starting incidence\nbin", 
      fill = "Starting incidence\nbin") 
  p <- ggpubr::ggarrange(month_incidence_p, 
                         month_prevlance_p,
                         month_metric_p, ncol=1, heights = c(1,1.15,2), 
                         common.legend = T, legend="right")
  
  ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_MonthlyTimescaleWHOBins_", names(y_names)[y_names == i], "_allSims.png"), 
         plot = p, path = plot_dir,
         width = 8, height = 5, units = c("in"), dpi = 300)
})


################################################################################
# Metrics on a single plot
################################################################################  
lapply(unique(month_summaries$start_who_bin), function(i, bin_group = "start_who_bin"){
  tmp_df <- dplyr::filter(month_summaries, .data[[bin_group]] == !!i)
  p <- MonthlyCombinedMetrics(tmp_df) 
  ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_MonthlyAllMetrics_Stratus", gsub(" ", "", i), ".png"), 
         plot = p, path = plot_dir,
         width = 7, height = 7, units = c("in"), dpi = 300)       
})


##### Line plot of differences
multisim_summary_long <- TimeScaleSummaryStats() 
lapply(unique(multisim_summary_long$start_who_bin), function(i){
  tmp_df <- dplyr::filter(multisim_summary_long, 
                          start_who_bin == !!i #& 
                          #importation_rate == "high"
  ) #%>% dplyr::filter(metric != "effective_coi_mean")
  #tmp_df$sampling_scheme <- gsub("Sample - ", "", tmp_df$sampling_scheme)
  if(nrow(tmp_df) > 1){
    p <- IndividualSamplingSchemeMeansPlot(tmp_df) + theme(panel.grid.major = element_blank())
    ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_SamplingSchemeVariation_Stratus", gsub(" ", "", i), "_AllYears.png"), 
           plot = p, path = plot_dir,
           width = 6.5, height = 4, units = c("in"), dpi = 300) 
    
    p_season <- IndividualSamplingSchemeMeansPlot(tmp_df, color_variable = "sampling_scheme", y_variable = "group_year", group_variable = "sampling_scheme")
    ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_SamplingSchemeVariation_Stratus", gsub(" ", "", i), "_ColorBySeason.png"), 
           plot = p_season, path = plot_dir,
           width = 6.5, height = 4, units = c("in"), dpi = 300) 
    
    yearly_change <-TimeSeriesMetricsPlot(dplyr::filter(multisim_summary_long, start_who_bin=="Low"), 
                          color_variable = "sampling_scheme",
                          plot_type = "change",
                          add_ribbon = TRUE)
    ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_SamplingSchemeVariation_Stratus", gsub(" ", "", i), "_LongitudinalDelta.png"), 
           plot = yearly_change, path = plot_dir,
           width = 6.5, height = 4, units = c("in"), dpi = 300) 
    
    yearly_abs <- TimeSeriesMetricsPlot(dplyr::filter(multisim_summary_long, start_who_bin=="Low"), 
                          color_variable = "sampling_scheme",
                          plot_type = "absolute",
                          add_ribbon = TRUE)
    ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_SamplingSchemeVariation_Stratus", gsub(" ", "", i), "_LongitudinalAbsolute.png"), 
           plot = yearly_abs, path = plot_dir,
           width = 6.5, height = 4, units = c("in"), dpi = 300)
    
  }
})


################################################################################
# Incidence plots 
################################################################################ 
# input_metric = y_vars[4]
# center_lines = c(1.4, 2.4)
input_metric = y_vars[2]
# center_lines = c(0.76, 0.90)
center_lines = c(0.73, 0.90)
input_incidence="reported_incidence_per_1k"

inc_plots <- IncidenceCurveBuild(year_summaries %>% dplyr::mutate(
  sampling_scheme = gsub("Sample - ", "", sampling_scheme), 
  sampling_scheme = factor(sampling_scheme, levels = names(sampling_colors_2))),
  incidence = input_incidence,
  metric = input_metric,
  #y_lim = c(1, max(year_summaries[[input_metric]])), 
  y_lim = c(0.3, 1),
  center_lines = center_lines,#[1],
  offset= 0.01
)

lapply(names(inc_plots), function(p){
  ggsave(paste(format(Sys.time(), "%Y%m%d"), p, names(y_names)[y_names == input_metric], "YearsAllPrevalence_log.png", sep="_"),
         plot = inc_plots[[p]],
         path = plot_dir,
         width = 5, height = 4, units = c("in"), dpi = 300)
})


points_df <- dplyr::bind_rows(lapply(center_lines, function(c){
  count_df <- dplyr::filter(year_summaries, between(.data[[input_metric]], c - 0.01, c + 0.01)) %>%
    dplyr::mutate(center_point = c)
  return(count_df)
}))

sampling_inc_summary <- dplyr::bind_rows(
  #SummarizeByGroups(
  #  data         = points_df,
  #  group_cols   = c("center_point"),
  #  metric_col   = incidence) %>% 
  #  dplyr::mutate(sampling_scheme = "Both schemes combined"),
  SummarizeByGroups(
    data         = points_df,
    group_cols   = c("center_point", "sampling_scheme"),
    metric_col   = input_incidence)
) %>%
  dplyr::filter(sampling_scheme != "All - Yearly") %>% 
  dplyr::mutate(
    sampling_scheme = gsub("Sample - ", "", sampling_scheme), 
    sampling_scheme = factor(sampling_scheme, levels = rev(names(sampling_colors_2))))


only_show_A <- TRUE
# Add a 'show' flag and keep all rows (to preserve x-axis)
plot_data <- sampling_inc_summary %>%
  mutate(
    center_point_chr = as.character(center_point),  # fix coercion issue
    center_line_chr = as.character(center_lines[1]),
    
    show = if (only_show_A) {
      center_point_chr == center_line_chr
    } else {
      TRUE
    },
    
    median_plot = ifelse(show, .data[[paste0("median_", input_incidence)]], NA),
    min_plot    = ifelse(show, .data[[paste0("min_", input_incidence)]], NA),
    max_plot    = ifelse(show, .data[[paste0("max_", input_incidence)]], NA)
  )


j <- position_dodge2(width = 0.4, 
                     preserve = "single",
                     reverse  = TRUE)

specific_points_p <- ggplot(plot_data, aes(
  x = as.character(center_point),
  y = .data[[ paste0("median_", input_incidence) ]],
  color = sampling_scheme,
  group = sampling_scheme
)) +
  make_transmission_background(cuts, transmission_colors, axis="y", max_val=650) + 
  
  # Plot error bars only for visible points
  geom_errorbar(aes(
    ymin = min_plot,
    ymax = max_plot
  ),
  size = 1,
  width = 0.4,
  position = j
  ) +
  
  # Plot points only for visible data
  geom_point(
    aes(y = median_plot),
    size = 3,
    position = j
  ) +
  
  scale_color_manual(values = sampling_colors_2) +
  labs(
    y = "Yearly reported incidence per 1k\n(Inferred)",
    x = paste(y_labs[names(y_labs) == input_metric], "(Measured)", sep="\n"),
    color = "Sampling scheme\n(100 samples)"
  ) +
  scale_x_discrete(
    breaks = center_lines,
    labels = LETTERS[1:length(center_lines)]
  ) +
  scale_y_continuous(breaks = c(100, 250, 450), limits = c(0,650)) +
  guides(
    size  = guide_legend(order = 2),
    color = guide_legend(order = 1, reverse = TRUE),
    shape = guide_legend(reverse = TRUE)
  ) +
  theme(
    legend.position = c(0.30, 0.80),
    legend.background = element_rect(fill = NA, color = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.key.size = unit(0.75, "lines")
  )
# specific_points_p <- sampling_inc_summary %>%
#   ggplot(aes(
#     x = as.character(center_point),
#     y = .data[[ paste0("median_", input_incidence) ]],
#     color = sampling_scheme,
#     group = sampling_scheme
#   )) +
#   make_transmission_background(cuts, transmission_colors, axis="y", max_val=550) + 
#   geom_errorbar(aes(
#     #ymin = .data[[paste0("mean_", incidence)]] - .data[[paste0("sd_", incidence)]],
#     #ymax = .data[[paste0("mean_", incidence)]] + .data[[paste0("sd_", incidence)]]), 
#     ymin = .data[[paste0("min_", input_incidence)]],
#     ymax = .data[[paste0("max_", input_incidence)]]),
#     # ymin = .data[[paste0("IQR_25_", incidence)]],
#     # ymax = .data[[paste0("IQR_75_", incidence)]]),
#     size = 1,
#     width = 0.4,
#     position = j) +
#   scale_color_manual(values = sampling_colors_2) +
#   geom_point(#aes(size = .data[[ paste0("count_", incidence) ]]),
#     size=3,         
#     position = j) +
#   labs(y = "Yearly reported incidence per 1k\n(Inferred)",
#        x = paste(y_labs[names(y_labs) == input_metric], "(Measured)", sep="\n"),
#        #size="Total\nsimulations",
#        color="Sampling scheme\n(100 samples)") +
#   scale_x_discrete(breaks=center_lines,
#                    labels=LETTERS[1:length(center_lines)]) +
#   scale_y_continuous(breaks = c(100, 250, 450)) +
#   guides(size  = guide_legend(order = 2),
#          color = guide_legend(order = 1, reverse = TRUE),
#          shape  = guide_legend(reverse = TRUE)
# ) +
#   #scale_y_log10() +
#   theme(
#     legend.position = c(0.30, 0.80),
#     #legend.background = element_rect(fill = "#F5F3ED"),
#     legend.background = element_rect(fill = NA, color = NA),
#     legend.key = element_rect(fill = NA, colour = NA),
#     legend.box       = "horizontal",    # <- stack multiple legends side-by-side
#     legend.box.just  = "left" ,
#     # 1) make the key boxes (and points inside) a bit bigger:
#     legend.key.size  = unit(0.75, "lines"),
#     # 2) add a little extra horizontal space between the two legends:
#   )
specific_points_p

ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_IncidenceMeansPointEstimates_", names(y_names)[y_names == input_metric], "YearsAll_RangeAOnly.png"),
       plot = specific_points_p,
       path = plot_dir,
       width = 3.5, height = 4, units = c("in"), dpi = 300)


################################################################################
# Statistics 
################################################################################  
# bin_order_vector <- levels(year_summaries$year_incidence_bin)
# # All metrics, per year 
# for(input_year in unique(year_summaries$year)){
#   res <- IncidencePairwiseCalc(
#     data       = dplyr::filter(year_summaries, year == input_year & importation_rate == "low"),
#     y_vars     = y_vars,                
#     group_vars = c("year_incidence_bin",     
#                    "sampling_scheme"),
#     adjust.method = "BH"
#   ) %>%
#     IncidencePairwiseEdit(.) %>%
#     dplyr::mutate(pair_var = factor(pair_var, levels = c(names(sampling_scheme_colors), "Different grouping variable")))
#   
#   ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_StatusBins50_Year", input_year, "LowImport.png"), 
#          #plot = IncidenceBinSignificancePlots(dplyr::filter(res, !grepl("Different", pair_var))), 
#          plot = IncidenceBinSignificancePlots(dplyr::filter(res, pair_var %in% sampling_subset)),
#          path = plot_dir,
#          width = 8, height = 6.5, units = c("in"), dpi = 300)
# }
# 
# 
# #### Check means for statistics
# multisim_summary_bybins <- TimeScaleSummaryStats(df = dplyr::filter(year_summaries, importation_rate == "low"), 
#                                                  bin="year_incidence_bin")
# 
# lapply(unique(multisim_summary_bybins$year), function(j){
#   df <- dplyr::filter(multisim_summary_bybins, year == j)
#   
#   count_p <- dplyr::filter(df, sampling_scheme == "All - Yearly") %>% 
#     select(year_incidence_bin, total_count) %>% unique() %>%
#     ggplot(aes(x=year_incidence_bin, y=total_count, fill=year_incidence_bin)) + geom_col() +
#     labs(x="", y="Simulations", fill = "Yearly\nincidence\nbin")
#   means_p <- IndividualSamplingSchemeMeansPlot(df, color_variable = "year_incidence_bin") +
#     labs(color = "Yearly\nincidence\nbin")
#   
#   p <- ggpubr::ggarrange(count_p, means_p, ncol=1, heights = c(.33,1),
#                          common.legend = T, legend="right")
#   
#   ggsave(paste0(format(Sys.time(), "%Y%m%d"), "MetricMeansBySamplingBins_Year", j, "LowImport.png"), 
#          plot = p, path = plot_dir,
#          width = 7, height = 5, units = c("in"), dpi = 300)
# })
# 
# # All metrics divided by year - all sampled infections
# res_by_year <- IncidencePairwiseCalc(
#   data       = year_summaries %>% filter(sampling_scheme == "All - Yearly"),
#   y_vars     = y_vars,                 
#   group_vars = c("year_incidence_bin",     
#                  "year"),
#   adjust.method = "BH"
# ) %>%
#   IncidencePairwiseEdit(.)
# 
# ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_StatusBins50_ByYear_AllInfections.png"), 
#        plot = IncidenceBinSignificancePlots(dplyr::filter(res_by_year, !grepl("Different", pair_var))), 
#        path = plot_dir,
#        width = 12, height = 6.5, units = c("in"), dpi = 300) 
# 
# # Comparing sampling schemes 
# input_year=0
# res <- IncidencePairwiseCalc(
#   data       = dplyr::filter(year_summaries, year == input_year),
#   y_vars     = y_vars,                 
#   group_vars = c("year_incidence_bin",     
#                  "sampling_scheme"),
#   adjust.method = "BH"
# ) %>%
#   IncidencePairwiseEdit(.) %>%
#   dplyr::mutate(pair_var = factor(pair_var, levels = c(names(sampling_scheme_colors), "Different grouping variable")))
