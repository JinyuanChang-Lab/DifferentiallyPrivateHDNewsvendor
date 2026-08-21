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
library(tidyr)
library(patchwork)
library(cowplot) 
 
############################################################################### 
##############################  2. Plots    ################################### 
############################################################################### 
here::i_am("scripts/figures/plot_sec5_examples.R")

#--------------------------------------------
#           plot for example 1
#--------------------------------------------

e1_csv <- here::here( "results", "figure_data", "Example1.csv" )
dat1 <- read.csv(e1_csv, check.names = FALSE)


## Plot
plot_df1 <- data.frame(
  p     = dat1[, 1],
  L0ERM = dat1[, 2],
  csERM = dat1[, 3],
  ERM   = dat1[, 5]
)

Oracle <- dat1[1, 4]
p_vec <- dat1$p

long1 <- plot_df1 %>%
  pivot_longer(-p, names_to = "Method", values_to = "Cost")
 
long1$Method <- factor(long1$Method, levels = c("L0ERM","csERM","ERM","Oracle"))

method_colors <- c(
  "L0ERM"   = "#2ca02c",
  "csERM"   = "deepskyblue",
  "ERM"     = "#4B3F72",
  "Oracle"  = "tomato"
)

method_lty <- c(
  "L0ERM"   = "solid",
  "csERM"   = "solid",
  "ERM"     = "solid",
  "Oracle"  = "longdash"
)

method_shapes <- c(
  "L0ERM"   = 16,
  "csERM"   = 16,
  "ERM"     = 16,
  "Oracle"  = NA
)

figure11 <- ggplot(long1, aes(x = p, y = Cost, color = Method)) +
  geom_line(aes(linetype = Method), linewidth = 1) +
  geom_point(aes(shape = Method), size = 2) +
  geom_hline(aes(yintercept = Oracle, color = "Oracle", linetype = "Oracle"),
             linewidth = 0.9) +
  scale_color_manual(
    values = method_colors,
    breaks = c("L0ERM", "csERM", "ERM", "Oracle"),
    labels = c("\u2113<sub>0</sub>-ERM", "cs-ERM", "ERM", "Oracle")
  ) +
  scale_linetype_manual(values = method_lty, breaks = names(method_lty)) +
  scale_shape_manual(values = method_shapes, breaks = names(method_shapes)) + 
  guides(
    color = guide_legend(
      title = "Method",
      override.aes = list(
        linetype = c("solid","solid","solid","longdash"),
        shape    = c(16,16,16,NA),
        linewidth = c(1,1,1,0.9)
      )
    ),
    linetype = "none",
    shape    = "none"
  ) +
  scale_x_continuous(breaks = p_vec, limits = range(p_vec)) +
  labs(x = "p (Dimension)", y = "Average out-of-sample cost") +
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_blank(),
    axis.line.x     = element_line(linewidth = 1.2, color = "black"),
    axis.line.y     = element_line(linewidth = 1.2, color = "black"),
    axis.title      = element_text(face = "bold"),
    axis.text       = element_text(face = "bold"),
    legend.title    = element_text(face = "bold", size = 13),
    legend.text     = element_markdown(size = 12),
    legend.position = "right",
    legend.key.width = grid::unit(1.2, "cm")   
  )

#figure11 
## Save the results  
figure_file_11 <- here::here( "results", "figures", "example_1.png" )
ggsave(filename = figure_file_11, plot = figure11, width = 6, height = 4, units = "in", dpi = 300)

#--------------------------------------------
#           plot for example 2
#--------------------------------------------

e2_csv <- here::here( "results", "figure_data", "Example2.csv" )
dat2 <- read.csv(e2_csv, check.names = FALSE)
epsilon_vec <- seq(0.1, 0.9, by = 0.1)
## Plot
plot_df2 <- data.frame(epsilon=dat2[,1],revenue= dat2[,2],cost=dat2[,3],profit=dat2[,4])
long2 <- plot_df2 %>% pivot_longer(-epsilon, names_to = "Method", values_to = "Cost") %>% mutate(Panel = "(b) Setting 2") 
method_colors <- c("revenue" = "#2ca02c","cost"   = "deepskyblue","profit"= "tomato") 
figure22 <- ggplot(long2, aes(x = epsilon, y = Cost, color = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = method_colors,
    breaks = c(  "revenue", "cost", "profit"),
    labels =c(  "Revenue", "Cost", "Profit")
  ) +
  scale_x_continuous(
    breaks = epsilon_vec,
    limits = range(epsilon_vec)
  ) +
  labs(
    x = "***epsilon*** **(Privacy level)**",
    y = "Average amount",
    color = "Method"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_blank(),
    axis.line.x  = element_line(linewidth = 1.2, color = "black"),
    axis.line.y  = element_line(linewidth = 1.2, color = "black"),
    axis.title   = element_text(face = "bold"), 
    axis.title.x = ggtext::element_markdown(face = "plain"), 
    axis.title.y = element_text(face = "bold"),
    axis.text    = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 12),
    legend.position = "right"
  )
 

#figure22 
## Save the results  
figure_file_22 <- here::here( "results", "figures", "example_2.pdf" ) 
ggsave(filename = figure_file_22, plot = figure22,  width = 6, height = 4, units = "in", dpi = 300)


#--------------------------------------------
#           plot for example 3
#--------------------------------------------

e3_csv <- here::here( "results", "figure_data", "Example3.csv" )
dat3 <- read.csv(e3_csv, check.names = FALSE) 


plot_df33 <- data.frame(
  epsilon       = as.numeric(dat3[, 1]),
  supplier_cost = as.numeric(dat3[, 2]),
  retailer_cost = as.numeric(dat3[, 3])
)


eps_ticks <- c(0.1, 0.5, 0.9)
 
y_breaks <- plot_df33$supplier_cost[match(eps_ticks, plot_df33$epsilon)]
 
if (any(is.na(y_breaks))) {
  y_breaks <- sapply(eps_ticks, function(e0) {
    plot_df33$supplier_cost[which.min(abs(plot_df33$epsilon - e0))]
  })
}
 
fit_quad <- lm(supplier_cost ~ poly(retailer_cost, 2, raw = TRUE), data = plot_df33)
xg <- seq(min(plot_df33$retailer_cost), max(plot_df33$retailer_cost), length.out = 200)
curve_df <- data.frame(
  retailer_cost = xg,
  supplier_cost = as.numeric(predict(fit_quad, newdata = data.frame(retailer_cost = xg)))
)
 
y_min <- min(plot_df33$supplier_cost, na.rm = TRUE)
y_max <- max(plot_df33$supplier_cost, na.rm = TRUE)
pad   <- 0.03 * (y_max - y_min)

figure33_quad <- ggplot(plot_df33, aes(x = retailer_cost, y = supplier_cost)) +
  geom_point(size = 2.2, color = "deepskyblue") +
  geom_line(data = curve_df,
            aes(x = retailer_cost, y = supplier_cost),
            inherit.aes = FALSE,
            linewidth = 1.3, color = "tomato") +
  labs(x = "Retailer's cost", y = "Supplier's cost") +
  scale_y_continuous(
    limits = c(y_min - pad, y_max + pad),
      
    sec.axis = sec_axis(
      transform = ~ .,
      name      = "*epsilon* (Privacy level)",
      breaks    = y_breaks,
      labels    = sprintf("%.1f", eps_ticks)
    )
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    
    axis.title = element_text(face = "bold", size = 18),
    axis.text  = element_text(face = "bold", size = 14),
    
    axis.line  = element_line(linewidth = 1.6, colour = "black"),
    axis.ticks = element_line(linewidth = 1.3, colour = "black"),
    axis.ticks.length = unit(0.22, "cm"),
     
    axis.line.y.right  = element_line(linewidth = 1.6, colour = "black"),
    axis.ticks.y.right = element_line(linewidth = 1.3, colour = "black"),
    axis.text.y.right  = element_text(face = "bold", size = 14, colour = "black"),
    axis.title.y.right = element_markdown(face = "bold", size = 18)
  )

#figure33_quad 
## Save the results  
figure_file_33 <- here::here( "results", "figures", "example_3.pdf" )
ggsave(filename = figure_file_33, plot = figure33_quad, width = 8, height = 6, units = "in", dpi = 300)




#--------------------------------------------
#           plot for example 4
#--------------------------------------------

e4_csv <- here::here( "results", "figure_data", "Example4.csv" )
dat4 <- read.csv(e4_csv, check.names = FALSE) 

n_vec  <- dat4[ ,1]                                   
# non-private 
cont <- dat4[,2]                                   # Excess loss (non-DP) under continuous DGP
mass_at_q <- dat4[,4]                              # Excess loss (non-DP) with mass at the quantile
mass_away <- dat4[,6]                              # Excess loss (non-DP) with mass away from the quantile
# private 
cont_priv <- dat4[,3]                              # Excess loss (DP) under continuous DGP
mass_at_q_priv <- dat4[,5]                         # Excess loss (DP) with mass at the quantile
mass_away_priv <- dat4[,7]                         # Excess loss (DP) with mass away from the quantile


## Plot
## ----------------------------
## styles
## ----------------------------
pt_cols <- c(cont = "deepskyblue", mass_away = "#2ca02c", mass_at_q = "tomato")
pt_pchs <- c(cont = 16,            mass_away = 17,        mass_at_q = 15)

## Baseline slopes colors (keep!)
ln_cols <- c("-1" = "#555555", "-1/2" = "#CCBB44")

## ----------------------------
## panel maker
## ----------------------------
make_panel <- function(cont_vec, mass_away_vec, mass_at_q_vec, title_text,
                       base_dash_key_width = grid::unit(0.7, "cm"),
                       base_dash_lty = "22") {
  
  dat <- tibble(
    n         = n_vec,
    cont      = cont_vec,
    mass_away = mass_away_vec,
    mass_at_q = mass_at_q_vec
  ) |>
    pivot_longer(-n, names_to = "series", values_to = "y") |>
    mutate(
      logn = log(n),
      logy = log(pmax(y, .Machine$double.eps))
    )
  
  base_slope <- c(cont = -1, mass_away = -1, mass_at_q = -0.5)
  
  anchors <- dat |>
    group_by(series) |>
    summarise(x0 = mean(logn), y0 = mean(logy), .groups = "drop") |>
    mutate(
      slope = base_slope[series],
      a     = y0 - slope * x0,
      slope_group = factor(ifelse(slope == -1, "-1", "-1/2"),
                           levels = c("-1", "-1/2"))
    )
  
  xgrid <- tibble(logn = seq(min(dat$logn), max(dat$logn), length.out = 200))
  
  base_lines <- anchors |>
    tidyr::crossing(xgrid) |>
    mutate(logy = a + slope * logn)
   
  dgp_labels <- list(
    cont      = "Continuous benchmark",
    mass_away = expression("Discontinuity away from " * q^{"*"}),
    mass_at_q = expression("Discontinuity at " * q^{"*"})
  )
  
  ggplot() +
    ## ----------------------------
  ## DGPs: points
  ## ----------------------------
  geom_point(
    data = dat,
    aes(x = logn, y = logy, color = series, shape = series),
    size = 2.8, alpha = 0.95
  ) +
    scale_color_manual(
      values = pt_cols,
      name   = "DGPs",
      breaks = c("cont", "mass_away", "mass_at_q"),
      labels = dgp_labels,
      guide = guide_legend(order = 2) 
    ) +
    scale_shape_manual(
      values = pt_pchs,
      name   = "DGPs",
      breaks = c("cont", "mass_away", "mass_at_q"),
      labels = dgp_labels,
      guide = guide_legend(order = 2)
    ) +
    
    ## new color scale for baseline slopes (so it has its own legend)
    ggnewscale::new_scale_color() +
    
    ## ----------------------------
  ## Baseline slopes: colored dashed lines (NO linetype mapping) 
  ## ----------------------------
  geom_line(
    data = base_lines,
    aes(x = logn, y = logy, group = series, color = slope_group),
    linewidth = 1.3,
    linetype  = base_dash_lty
  ) +
    scale_color_manual(
      values = ln_cols,
      name   = "Baseline slopes",
      breaks = c("-1", "-1/2"),
      labels = c("-1", "-1/2")
    ) +
    guides(
      color = guide_legend(
        order = 1,
        keywidth = base_dash_key_width,                          
        override.aes = list(linetype = base_dash_lty, linewidth = 1.3)
      )
    ) +
    
    ## ----------------------------
  ## labels / theme
  ## ----------------------------
  labs(title = title_text, x = "log(n)", y = "log(value)") +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right",
      
      legend.title = ggplot2::element_text(face = "bold", size = 12),
      legend.text  = ggplot2::element_text(size = 12),
      
      legend.spacing.y   = grid::unit(12, "pt"),
      legend.key.height  = grid::unit(18, "pt"),
      legend.key.width   = grid::unit(16, "pt"),
      legend.box.spacing = grid::unit(12, "pt"),
      
      plot.title = ggtext::element_markdown(size = 14, hjust = 0.5),
      
      axis.title.x = ggplot2::element_text(
        face = "bold", size = 13,
        margin = ggplot2::margin(t = 6)
      ),
      axis.title.y = ggplot2::element_text(
        face = "bold", size = 13,
        margin = ggplot2::margin(r = 6)
      ),
      axis.text.x = ggplot2::element_text(face = "bold", size = 12),
      axis.text.y = ggplot2::element_text(face = "bold", size = 12),
      axis.line   = ggplot2::element_line(linewidth = 1.1),
      axis.ticks  = ggplot2::element_line(linewidth = 0.9),
      
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

## ----------------------------
## panels
## ----------------------------
p_public  <- make_panel(
  cont, mass_away, mass_at_q,
  title_text = "&#8467;<sub>0</sub>-ERM",
  base_dash_key_width = grid::unit(0.7, "cm"),    
  base_dash_lty = "22"                      
)

p_private <- make_panel(
  cont_priv, mass_away_priv, mass_at_q_priv,
  title_text = "DP &#8467;<sub>0</sub>-ERM",
  base_dash_key_width = grid::unit(0.7, "cm"),
  base_dash_lty = "22"
)

p_public  <- p_public  + labs(x = NULL, y = NULL)
p_private <- p_private + labs(x = NULL, y = NULL)

combined <- (p_public | p_private) +
  plot_layout(guides = "collect") &
  theme(
    plot.margin = ggplot2::margin(t = 6, r = 12, b = 6, l = 10),
    legend.spacing.y  = unit(10, "pt"),
    legend.key.height = unit(16, "pt"),
    legend.key.width  = unit(14, "pt")
  )

figure44 <- ggdraw() +
  draw_plot(combined, x = 0.07, y = 0.10, width = 0.93, height = 0.88) +
  draw_label("log(n)", x = 0.45, y = 0.06, vjust = 0.5, size = 14, fontface = "bold") +
  draw_label("Logarithmic average excess risk", x = 0.03, y = 0.52, angle = 90,
             vjust = 0.5, size = 14, fontface = "bold") +
  theme(plot.background = element_rect(fill = "white", colour = NA))

#figure44


figure_file_44 <- here::here( "results", "figures", "example_4.png" )
ggsave(filename = figure_file_44, plot = figure44, width = 11, height = 5, units = "in", dpi = 300) 
