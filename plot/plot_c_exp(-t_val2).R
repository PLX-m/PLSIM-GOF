library(ggplot2)
library(patchwork)

df <- data.frame(
  c = rep(c(0, 0.25, 0.5, 0.75, 1), times = 8),
  value = c(
    # n=800, d=s=10：Z1, Z2, CCT, Ln
    0.054, 0.932, 0.954, 0.956, 0.955,
    0.063, 0.922, 0.953, 0.959, 0.962,
    0.070, 0.953, 0.972, 0.970, 0.979,
    0.036, 0.396, 0.646, 0.700, 0.766,
    
    # n=800, d=s=20：Z1, Z2, CCT, Ln
    0.059, 0.799, 0.862, 0.891, 0.892,
    0.057, 0.758, 0.880, 0.889, 0.899,
    0.063, 0.842, 0.933, 0.934, 0.934,
    0.056, 0.262, 0.486, 0.546, 0.622
  ),
  Method = rep(c("Z1", "Z2", "Tc", "Ln"), each = 5, times = 2),
  scenario = rep(c("d10", "d20"), each = 20)
)

# ==============================
# 字体大小参数
# ==============================
base_size        <- 16  # theme_bw 的基础字体大小
title_size       <- 16  # 子图标题字体大小
axis_title_size  <- 16  # 横纵坐标标题字体大小
axis_text_size   <- 14  # 坐标轴刻度字体大小
legend_text_size <- 14  # 图例文字字体大小
legend_title_size <- 16 # 图例标题字体大小


# my_cols <- c(
#   "Z1" = "#00BFC4",
#   "Z2" = "#E69F00",
#   "Tc" = "#F8766D"
# )

my_cols <- c(
  "Z1" = "#3B82B8",
  "Z2" = "#E69F00",
  "Tc" = "#D55E5E",
  "Ln" = "#59A14F"
)



# ==============================
# d=10 图
# ==============================
p1 <- ggplot(subset(df, scenario == "d10"),
             aes(x = c, y = value, color = Method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = my_cols,
    breaks = c("Z1", "Z2", "Tc", "Ln"),
    labels = c(
      expression(Z[n]^{(1)}),
      expression(Z[n]^{(2)}),
      expression(T[C]),
      expression(L[n])
    )
  )+
  # geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "c",   
    y = "Probability of rejection",
    title = "p=q=10"
  ) +
  theme_bw(base_size = base_size) +
  theme(
    plot.title   = element_text(size = title_size),        # 子图标题
    axis.title   = element_text(size = axis_title_size),   # 横纵坐标标题
    axis.text    = element_text(size = axis_text_size),    # 坐标轴刻度
    legend.text  = element_text(size = legend_text_size),  # 图例文字
    legend.title = element_text(size = legend_title_size)  # 图例标题
    # panel.grid = element_blank()   # 去掉网格线
    
  )

# ==============================
# d=20 图
# ==============================
p2 <- ggplot(subset(df, scenario == "d20"),
             aes(x = c, y = value, color = Method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = my_cols,
    breaks = c("Z1", "Z2", "Tc", "Ln"),
    labels = c(
      expression(Z[n]^{(1)}),
      expression(Z[n]^{(2)}),
      expression(T[C]),
      expression(L[n])
    )
  ) +
  # geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "c",
    y = "Probability of rejection",
    title = "p=q=20"
  ) +
  theme_bw(base_size = base_size) +
  theme(
    plot.title   = element_text(size = title_size),        # 子图标题
    axis.title   = element_text(size = axis_title_size),   # 横纵坐标标题
    axis.text    = element_text(size = axis_text_size),    # 坐标轴刻度
    legend.text  = element_text(size = legend_text_size),  # 图例文字
    legend.title = element_text(size = legend_title_size)  # 图例标题
    # panel.grid = element_blank()   # 去掉网格线
    
  )

# ==============================
# 拼接（横向）
# ==============================
p1 + p2


setwd("plot/exp(-t_val2)/")
#
p_final <- (p1 + p2) + plot_layout(guides = "collect")
#
ggsave(
  filename = "power_plot_exp(-t_val2).pdf",
  plot = p_final,
  width = 12,
  height = 5,
  units = "in",
  device = cairo_pdf   # 字体更清晰
)


