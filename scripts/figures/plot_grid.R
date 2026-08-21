############################################################################### 
###########################  1. Preparation    ################################ 
############################################################################### 

#-------------------------------------------------------------------------------
#                          Required packages                         
#------------------------------------------------------------------------------- 
library(here)
library(ggplot2)
library(patchwork)
library(ggtext)

############################################################################### 
##############################  2. Plots    ################################### 
############################################################################### 
here::i_am("scripts/figures/plot_grid.R")
grid_csv <- here::here( "results", "figure_data", "Result_grid.csv" ) 
data_grid <- read.csv(grid_csv, check.names = FALSE) 

p_vec <-  data_grid$p 
Linear_1 <-  data_grid$linear_SQR
Expand_1 <-  data_grid$Expand_SQR
Bin_1 <-  data_grid$Bin
Knn_1 <-  data_grid$Knn 


plot_df_1 <- data.frame(
  p    = p_vec,
  Linear = Linear_1,
  Expand = Expand_1,
  Bin = Bin_1,
  Knn = Knn_1
)
 

make_long <- function(df_wide) {
  data.frame(
    p = rep(df_wide$p, times = 4),
    Method = factor(
      rep(c("Linear", "Expand", "Bin", "Knn"), each = nrow(df_wide)),
      levels = c("Linear", "Expand", "Bin", "Knn")
    ),
    Cost = c(df_wide$Linear, df_wide$Expand, df_wide$Bin, df_wide$Knn)
  )
}

plot_df_long_1 <- make_long(plot_df_1) 

 

method_colors <- c(
  "Linear" = "tomato",
  "Expand" = "deepskyblue",
  "Bin" = "grey40",
  "Knn" = "#2ca02c"
)

method_labels <- c(
  "Linear" = "\u2113<sub>0</sub>-ERM",
  "Expand" = "Expanded-\u2113<sub>0</sub>-ERM",
  "Bin" = "Binning",
  "Knn" = "KNN"
)

base_theme <- theme_classic(base_size = 14) +
  theme(
    axis.line.x   = element_line(linewidth = 0.8, color = "black"),
    axis.line.y   = element_line(linewidth = 0.8, color = "black"),
    axis.title.x  = element_text(face = "bold", size = 14),
    axis.title.y  = element_text(face = "bold", size = 14),
    axis.text     = element_text(face = "bold", size = 12),
    plot.title    = element_blank(),
    legend.text   = ggtext::element_markdown(size = 12),
    legend.title  = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

p1 <- ggplot(plot_df_long_1, aes(x = p, y = Cost, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = method_colors,
    breaks = c("Linear", "Expand", "Bin", "Knn"),
    labels = method_labels
  ) +
  scale_x_continuous(
    breaks = seq(5, 30, by = 5)
  ) +
  scale_y_continuous(
    breaks =  seq(1.5, 3.0, by = 0.25)
  )+
  labs(
    x = "p",
    y = "Average out-of-sample cost", 
    color = "Method"
  ) +
  base_theme


#p1

  
figure_file_grid <- here::here( "results", "figures", "figure_grid.png" )
ggsave(filename = figure_file_grid, plot = p1, width = 5.5, height = 3.5, units = "in", dpi = 300)
