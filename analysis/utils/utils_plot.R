#' Functions needed for plots

#' Visualize average observed treatment effects per trial and per visit, by 
#' using a last observation carried forward (LOCF) principle
#'
#' @param data_locf_effects A data frame average treatment effects
#' @return A plot that visualizes average observed treatment effects per trial 
#'   and per visit
plot_average_effects_locf <- function(data_locf_effects) {
  
  plot <- ggplot(data_locf_effects, aes(x = week, y = estimate, group = trial)) +
    
    geom_errorbar(
      aes(ymin = lower_ci, ymax = upper_ci, colour = as.factor(trial)), 
      width = 0.5, 
      size = 0.3
    ) + 
    geom_point(
      aes(colour = as.factor(trial), fill = as.factor(trial)), 
      shape = 19, 
      size = 1.5, 
      stroke = 1
    ) + 
    geom_line(aes(colour = as.factor(trial)), linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = 'dotted', linewidth = 0.3) +
    scale_color_manual(values=c("#F46D43", "#66C2A5", "#008080")) + 
    labs(x = "Week", y = "Difference in MDS-UPDRS 2+3 change from baseline") +
    ggtitle("Average effect of rasagiline on MDS-UPDRS part 2+3 score\n") +
    coord_cartesian(ylim = c(-8.5, 0)) +
    scale_y_continuous(breaks = c(-8, -6, -4, -2, 0)) +
    scale_x_continuous(breaks = c(0, 6, 10, 12, 14, 20, 24, 26, 36)) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 9, colour = 'black'), 
          axis.text.y = element_text(size = 9, colour = 'black'),
          axis.title.x = element_text(size = 10), 
          axis.title.y = element_text(size = 10),
          plot.title = element_text(size = 10, hjust = 0.5), 
          legend.title = element_blank())
  return(plot)
}

#' Plot median feature importance in a bar plot
#'
#' @param median_feature_importances A data frame with median feature importances
#' @param fill_color A string for the color of the bars
#' @param nb_features_to_show A number indicating how many most important 
#'   features to show 
#' @param ymax A number indicating the limit of the y-axis
#' @return A bar plot of median feature importances
plot_median_feature_importance <- function(median_feature_importance, 
                                           fill_color = "#4B4BCD", 
                                           nb_features_to_show = 10, 
                                           y_max = 0.125) {
  
  data_plot <- median_feature_importance |>
    arrange(desc(median)) |>
    head(nb_features_to_show)
  
  plot <- ggplot(data_plot, aes(x = reorder(feature_text, median), y = median)) +
    geom_bar(
      stat = "identity", 
      fill = fill_color, 
      width=.7, 
      position = position_dodge(width = 0.9)
    ) +
    geom_errorbar(aes(
      ymin = lower_quartile, ymax = upper_quartile), 
      linewidth = 0.2, 
      width = .2, 
      position = position_dodge(.9)
    ) +
    coord_flip() +
    xlab("") +
    ylab("Importance") +
    ylim(c(0, y_max)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size=10), 
      axis.text.y = element_text(size=10),
      axis.title.x = element_text(size=11), 
      axis.title.y = element_text(size=11)
    )
  
  return(plot)
}

#' Plot observed average treatment effect (ATE) per quantile of a given 
#' prioritization rule
#'
#' @param quantile_effects A data frame with quantile effects
#' @param title A string for the plot title
#' @param x_label A string for the x-axis label
#' @param y_min A number indicating the min limit of the y-axis
#' @param y_max A number indicating the max limit of the y-axis
#' @return A plot of quantile effects
plot_quantile_effects <- function(quantile_effects, 
                                  title, 
                                  x_label, 
                                  y_min = -8, 
                                  y_max = 1) {
  
  data_plot <- quantile_effects |>
    mutate(xtick_text = paste(quantile_text, "\n\n\n", 
                              round(estimate, digits=1), "\n", 
                              round(prior_mean, digits=1), "\n",
                              nb_people))
  
  plot <- ggplot(data_plot, aes(x = xtick_text, y = estimate)) +
    geom_hline(yintercept = 0, linetype = 'dotted', linewidth = 0.3) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.12) +
    geom_point(size = 1.5) +
    ggtitle(title) +
    labs(
      x = x_label, 
      y = "\n\n\n\nObserved treatment effect
      (Difference in MDS-UPDRS change from baseline)"
    ) +
    
    annotate("text", x = 0.3, y = -9.5, label = "Mean difference", size = 3.1, hjust = 1) +
    annotate("text", x = 0.3, y = -10, label = "observed", size = 3.1, hjust = 1) +   
    annotate("text", x = 0.3, y = -10.5, label = "estimated", size = 3.1, hjust = 1) +  
    annotate("text", x = 0.3, y = -11, label = "Number of patients", size = 3.1, hjust = 1) +
    
    coord_cartesian(ylim = c(y_min, y_max), clip = "off") + 
    scale_y_continuous(breaks = c(-6, -4, -2, 0)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size=9, colour='black'), 
      axis.text.y = element_text(size=9, colour='black'),
      axis.title.x = element_text(size=10), 
      axis.title.y = element_text(size=10),
      plot.title = element_text(size=11, hjust=0.5)
    )
  
  return(plot)
}

#' Plot observed average treatment effect (ATE) per quantile of a given 
#' single baseline characteristic
#'
#' @param quantile_effects A data frame with quantile effects
#' @param title A string for the plot title
#' @param x_label A string for the x-axis label
#' @param y_min A number indicating the min limit of the y-axis
#' @param y_max A number indicating the max limit of the y-axis
#' @return A plot of quantile effects
plot_quantile_effects_single <- function(quantile_effects, 
                                         title, 
                                         x_label, 
                                         y_min = -12, 
                                         y_max = 2) {
  
  data_plot <- quantile_effects |>
    mutate(xtick_text = paste(quantile_text, "\n\n", 
                              round(estimate, digits=1), "\n", 
                              round(prior_mean, digits=1), "\n",
                              nb_people))
  
  plot <- ggplot(data_plot, aes(x = xtick_text, y = estimate)) +
    geom_hline(yintercept = 0, linetype = 'dotted', linewidth = 0.3) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.12) +
    geom_point(size = 1.5) +
    ggtitle(title) +
    labs(
      x = x_label, 
      y = "\n\n\nObserved treatment effect
      (Difference in MDS-UPDRS change from baseline)"
    ) +
    
    annotate("text", x = 0.3, y = -14.5, label = "Mean difference", size = 3.1, hjust = 1) +
    annotate("text", x = 0.3, y = -15.2, label = "Mean characteristic", size = 3.1, hjust = 1) +   
    annotate("text", x = 0.3, y = -15.9, label = "Number of patients", size = 3.1, hjust = 1) +
    
    coord_cartesian(ylim = c(y_min, y_max), clip = "off") + 
    theme_classic() +
    theme(
      axis.text.x = element_text(size=9, colour='black'), 
      axis.text.y = element_text(size=9, colour='black'),
      axis.title.x = element_text(size=10), 
      axis.title.y = element_text(size=10),
      plot.title = element_text(size=11, hjust=0.5)
    )
  return(plot)
}


#' Plot observed average treatment effect (ATE) per quantile of a given 
#' prioritization rule as a binned calibration curve
#'
#' @param quantile_effects A data frame with quantile effects
#' @param title A string for the plot title
#' @param min_axis A number indicating the min limit of the x- and y-axis
#' @param max_axis A number indicating the max limit of the x- and y-axis
#' @return A binned calibration curve plot of quantile effects
plot_quantile_calibration_curve <- function(quantile_effects, 
                                            title, 
                                            min_axis = -5, 
                                            max_axis = 0) {
  
  # include only quantile effects, not the overall effect
  data_plot <- quantile_effects |>
    filter(!quantile == 0)
  
  plot <- ggplot(data_plot, aes(x = prior_mean, y = estimate)) +
    geom_abline(intercept = 0, linetype = 'dotted', linewidth = 0.3) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.12) +
    geom_point(shape = 21, size = 1, stroke = 1, color = 'black', fill = 'black') +
    geom_line() +
    ggtitle(title) +
    labs(
      x = "Estimated treatment effect
      (Difference in MDS-UPDRS change from baseline)", 
      y = "Observed treatment effect
      (Difference in MDS-UPDRS change from baseline)"
    ) +
    coord_cartesian(
      xlim = c(min_axis, max_axis), 
      ylim = c(min_axis, max_axis), 
      clip = "off"
    ) + 
    theme_classic() +
    theme(
      axis.text.x = element_text(size=9, colour='black'), 
      axis.text.y = element_text(size=9, colour='black'),
      axis.title.x = element_text(size=10), 
      axis.title.y = element_text(size=10),
      plot.title = element_text(size=11, hjust=0.5)
    )
  
  return(plot)
}

#' Plot partial dependence for a number of important features
#'
#' @param partial_dependence A data frame with partial dependences
#' @param feature_importance A data frame with feature importances
#' @param feature_map A data frame that maps feature codes to feature texts
#' @param color A string for the color of the lines
#' @param ymin A number indicating the min limit of the y-axis
#' @param ymax A number indicating the max limit of the y-axis
#' @return A binned calibration curve plot of quantile effects
plot_partial_dependence <- function(partial_dependence, 
                                    feature_importance, 
                                    feature_map, 
                                    color = "#4B4BCD", 
                                    ymin = -3.2,
                                    ymax = -2.2) {
  
  importance_rank <- feature_importance |>
    mutate(rank = rank(-median, ties.method = "first")) |>
    mutate(feature_text_unit = recode_values(
      feature_code, 
      from = feature_map$feature_code, 
      to = feature_map$feature_text_unit))
  
  levels <- importance_rank$feature_text_unit[order(importance_rank$rank)]
  
  data_plot <- partial_dependence |>
    # map variable codes to full text definitions
    mutate(feature_text_unit = recode_values(
      feature_code, 
      from = feature_map$feature_code, 
      to = feature_map$feature_text_unit)) |>
    # add the ranking of the feature
    mutate(rank = recode_values(
      feature_code,
      from = importance_rank$feature_code,
      to = importance_rank$rank)) |>
    mutate(feature_text_unit_fac = as.factor(feature_text_unit)) |>
    mutate(feature_text_unit_fac = factor(
      feature_text_unit_fac, 
      levels = levels))
  
  plot <- ggplot(data = data_plot) +
    geom_line(
      aes(x = grid_value, y = pd_median), 
      color = color, 
      linewidth = 0.5
    ) +
    geom_ribbon(
      aes(x = grid_value, ymin  = pd_lower_quartile, ymax = pd_upper_quartile), 
      alpha = 0.2, 
      fill = color
    ) +
    facet_wrap(
      (vars(feature_text_unit_fac)), 
      ncol = 3, 
      scales = "free", 
      strip.position = "bottom"
    ) +
    coord_cartesian(ylim = c(ymin, ymax), clip = "off") +
    labs(x = "", y = "Mean estimated individualized treatment effect") +
    theme_classic() +
    theme(
      strip.background = element_blank(), 
      strip.placement = "outside",
      strip.text = element_text(margin = margin(b = 20))
    )
  return(plot)
}

