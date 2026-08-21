############################################################################### 
###########################  1. Preparation    ################################ 
############################################################################### 

#-------------------------------------------------------------------------------
#                          Required packages                         
#------------------------------------------------------------------------------- 
library(here)
library(dplyr)
library(ggplot2)
library(ggtext)
library(grid) 

############################################################################### 
##############################  2. Plots    ################################### 
############################################################################### 
here::i_am("scripts/figures/plot_LDP.R")
LDP_csv <- here::here( "results", "figure_data", "Result_LDP.csv" ) 
datLDP <- read.csv(LDP_csv, check.names = FALSE)

n_vec <- seq(4000,60000,by=4000)  


color_order  <- c("CDP", "HDP(K=2)","HDP(K=4)", "LDP")
color_values <- c("CDP" = "tomato",
                  "HDP(K=2)" = "#4B3F72",
                  "HDP(K=4)" = "#2ca02c",
                  "LDP" = "deepskyblue")

## -----------------------------
theme_comb <- function(a){
  if(a == 1){
    theme(
      text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14, color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 1),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.15, "cm"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
    )
  } else {
    theme(
      legend.position = "none",
      text = element_text(size = 10),
      axis.text = element_text(size = 10, color = "black"),
      axis.line.x = element_line(color = "black", linewidth = 1),
      axis.line.y = element_line(color = "black", linewidth = 1),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.15, "cm")
    )
  }
}

legend_title_expr <- expression(bold("Method"))
line_size_main <- 1.2

df_excess <- rbind(
  data.frame(n = n_vec, error = log(as.numeric(datLDP[,1])), type = "CDP"),
  data.frame(n = n_vec, error = log(as.numeric(datLDP[,2])), type = "HDP(K=2)"),
  data.frame(n = n_vec, error = log(as.numeric(datLDP[,3])), type = "HDP(K=4)"),
  data.frame(n = n_vec, error = log(as.numeric(datLDP[,4])), type = "LDP")
) %>%
  mutate(type = factor(type, levels = color_order))

x_breaks <- n_vec[seq(1, length(n_vec), by = 2)]

p_excess <- ggplot(df_excess, aes(x = n, y = error, color = type, linetype = type)) +
  geom_line(linewidth = line_size_main) +
  geom_point(size = 1.6) +
  labs(
    x = "Sample size (n/100)",
    y = "Logarithmic average excess risk"
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_breaks / 100
  ) +
  scale_color_manual(
    name = legend_title_expr,
    values = color_values,
    breaks = color_order,
    limits = color_order
  ) +
  scale_linetype_manual(
    name = legend_title_expr,
    values = c(
      "CDP" = "solid",
      "HDP(K=2)" = "solid",
      "HDP(K=4)" = "solid",
      "LDP" = "solid"
    ),
    breaks = color_order,
    limits = color_order
  ) +
  theme_minimal() +
  theme_comb(1) +
  theme(
    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(0.55, "cm")
  ) +
  guides(
    linetype = guide_legend(override.aes = list(linewidth = 1.1))
  )

#p_excess

figure_file_LDP <- here::here( "results", "figures", "figure_LDP.pdf" )
ggsave(filename = figure_file_LDP, plot = p_excess, width = 5.9, height = 4.2, units = "in", dpi = 300)
