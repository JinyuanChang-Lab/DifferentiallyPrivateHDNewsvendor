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

 

#######################################################################
## --------------------- 2.1 non-private results --------------------##
#######################################################################
here::i_am("scripts/figures/plot_NN.R")
NN_np_csv <- here::here( "results", "figure_data", "Result_nn_np.csv" ) 
dat_np <- read.csv(NN_np_csv, check.names = FALSE)

keep_cols <- c("n","ER_nn_np","ER_saa","ER_qrf","ER_qrnn")

df_all_1 <- dat_np[dat_np$model == 1,keep_cols,drop = FALSE]

df_all_2 <- dat_np[dat_np$model == 2,keep_cols,drop = FALSE]

method_names <- c("ER_nn_np", "ER_qrf", "ER_qrnn", "ER_saa")

make_long <- function(df_all) {
  df_long <- data.frame(
    n = rep(df_all$n, times = length(method_names)),
    Method = rep(method_names, each = nrow(df_all)),
    Value = c(df_all$ER_nn_np, df_all$ER_qrf, df_all$ER_qrnn, df_all$ER_saa)
  )
  
  df_long$n <- factor(df_long$n, levels = c(500, 1000, 1500, 2000))
  df_long$Method <- factor(df_long$Method, levels = method_names)
  df_long
}

df_long_1 <- make_long(df_all_1)
df_long_2 <- make_long(df_all_2)

method_colors <- c(
  "ER_nn_np"   = "tomato",
  "ER_qrf"  = "#2ca02c",
  "ER_qrnn" = "#4B3F72",
  "ER_saa"  = "deepskyblue"
)

base_theme <- theme(
  panel.grid = element_blank(),
  axis.text.x = element_text(size = 12, face = "bold",color = "black"),
  axis.text.y = element_text(face = "bold", size = 12,color = "black"),
  axis.title.y = element_text(face = "bold", size = 14),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
  panel.border = element_blank(),
  axis.line = element_line(color = "black", linewidth = 0.8),
  axis.ticks = element_line(color = "black"),
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA),
  legend.text = element_text(size = 11, face = "bold"),
  legend.title = element_text(size = 12, face = "bold")
)

pp1 <- ggplot(df_long_1, aes(x = n, y = Value, fill = Method)) +
  geom_boxplot(
    alpha = 0.7,
    outlier.size = 0.2,
    linewidth = 0.3,
    position = position_dodge2(width = 0.75, preserve = "single")
  ) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5),
    linetype = "dashed",
    color = "grey40"
  ) +
  scale_fill_manual(
    values = method_colors,
    labels = c("FNN", "QRF", "QRNN", "SAA")
  ) +
  scale_x_discrete(
    labels = c(
      "500"  = "n = 500",
      "1000" = "n = 1000",
      "1500" = "n = 1500",
      "2000" = "n = 2000"
    )
  ) +
  labs(x = NULL, y = "Boxplots of excess risk", fill = "Method", title = "Model 1") +
  base_theme +
  coord_cartesian(ylim = c(0, 0.15))

pp2 <- ggplot(df_long_2, aes(x = n, y = Value, fill = Method)) +
  geom_boxplot(
    alpha = 0.7,
    outlier.size = 0.2,
    linewidth = 0.3,
    position = position_dodge2(width = 0.75, preserve = "single")
  ) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5),
    linetype = "dashed",
    color = "grey40"
  ) +
  scale_fill_manual(
    values = method_colors,
    labels = c("FNN", "QRF", "QRNN", "SAA")
  ) +
  scale_x_discrete(
    labels = c(
      "500"  = "n = 500",
      "1000" = "n = 1000",
      "1500" = "n = 1500",
      "2000" = "n = 2000"
    )
  ) +
  labs(x = NULL, y = "Boxplots of excess risk", fill = "Method", title = "Model 2") +
  base_theme +
  coord_cartesian(ylim = c(0, 0.35))    

figure_np <- (pp1 | pp2) + plot_layout(guides = "collect") & theme(legend.position = "right")
#figure_np

figure_file_nn <- here::here( "results", "figures", "figure_NN_np.pdf" )
ggsave(filename = figure_file_nn, plot = figure_np, width = 12, height = 4.5, units = "in", dpi = 300)



#######################################################################
## ------------------------ 2.2 private results ---------------------##
#######################################################################

NN_priv_csv <- here::here( "results", "figure_data", "Result_nn_priv.csv" ) 
data_priv  <- read.csv(NN_priv_csv, check.names = FALSE)
colnames(data_priv)  <- c("model","n", "nonp","epsilon=0.5","epsilon=0.7","epsilon=0.9")
data_priv_model_1 <- data_priv[1:10,]
data_priv_model_2 <- data_priv[11:20,] 



n_vec <-  data_priv_model_1$n 
nonp_1 <- data_priv_model_1$nonp
eps_05_1 <- data_priv_model_1$`epsilon=0.5`
eps_07_1 <- data_priv_model_1$`epsilon=0.7`
eps_09_1 <- data_priv_model_1$`epsilon=0.9`
nonp_2 <- data_priv_model_2$nonp
eps_05_2 <- data_priv_model_2$`epsilon=0.5`
eps_07_2 <- data_priv_model_2$`epsilon=0.7`
eps_09_2 <- data_priv_model_2$`epsilon=0.9`


plot_df_1 <- data.frame(
  n    = n_vec,
  nonp = nonp_1,
  eps5 = eps_05_1,
  eps7 = eps_07_1,
  eps9 = eps_09_1
)

plot_df_2 <- data.frame(
  n    = n_vec,
  nonp = nonp_2,
  eps5 = eps_05_2,
  eps7 = eps_07_2,
  eps9 = eps_09_2
)

make_long <- function(df_wide) {
  data.frame(
    n = rep(df_wide$n, times = 4),
    Method = factor(
      rep(c("nonp", "eps5", "eps7", "eps9"), each = nrow(df_wide)),
      levels = c("nonp", "eps5", "eps7", "eps9")
    ),
    Cost = c(df_wide$nonp, df_wide$eps5, df_wide$eps7, df_wide$eps9)
  )
}

plot_df_long_1 <- make_long(plot_df_1)
plot_df_long_2 <- make_long(plot_df_2)


plot_df_long_1$x_show <- plot_df_long_1$n / 1000
plot_df_long_2$x_show <- plot_df_long_2$n / 1000


method_colors <- c(
  "nonp" = "tomato",
  "eps5" = "deepskyblue",
  "eps7" = "grey40",
  "eps9" = "#2ca02c"
)

method_labels <- c(
  "nonp" = "<b>Nonprivate</b>",
  "eps5" = "<b>*epsilon* = 0.5</b>",
  "eps7" = "<b>*epsilon* = 0.7</b>",
  "eps9" = "<b>*epsilon* = 0.9</b>"
)


base_theme <- theme_classic(base_size = 14) +
  theme(
    axis.line.x   = element_line(linewidth = 0.8, color = "black"),
    axis.line.y   = element_line(linewidth = 0.8, color = "black"),
    axis.title.x  = element_text(face = "bold", size = 14),
    axis.title.y  = element_text(face = "bold", size = 14),
    axis.text     = element_text(face = "bold", size = 12),
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.text = ggtext::element_markdown(),
    legend.title  = element_blank(),
    legend.position = "right"
  )

#-------------------------
# Model 1
#-------------------------
p1 <- ggplot(plot_df_long_1, aes(x = x_show, y = Cost, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = method_colors,
    breaks = c("nonp", "eps5", "eps7", "eps9"),
    labels = method_labels
  ) +
  scale_x_continuous(
    breaks = sort(unique(plot_df_long_1$x_show))
  ) +
  labs(
    x = "Sample size (n/1000)",
    y = "Average excess risk",
    title = "Model 1"
  ) +
  base_theme

#-------------------------
# Model 2
#-------------------------
p2 <- ggplot(plot_df_long_2, aes(x = x_show, y = Cost, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = method_colors,
    breaks = c("nonp", "eps5", "eps7", "eps9"),
    labels = method_labels
  ) +
  scale_x_continuous(
    breaks = sort(unique(plot_df_long_2$x_show))
  ) +
  labs(
    x = "Sample size (n/1000)",
    y = "Average excess risk",
    title = "Model 2"
  ) +
  base_theme

figure_priv <- (p1 | p2) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.text = ggtext::element_markdown()
  ) 

#figure_priv
  
figure_file_nn_priv <- here::here( "results", "figures", "figure_NN_priv.pdf" )
ggsave(filename = figure_file_nn_priv, plot = figure_priv, width = 12, height = 4.5, units = "in", dpi = 300)

