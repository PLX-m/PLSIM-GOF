library(ggplot2)
library(patchwork)

df <- data.frame(
  c = rep(c(0, 0.25, 0.5, 0.75, 1), times = 8),
  value = c(
    # n=800, d=s=10：Z1, Z2, CCT, Ln
    0.051, 0.940, 0.987, 0.986, 0.978,
    0.045, 0.954, 0.985, 0.975, 0.979,
    0.047, 0.974, 0.994, 0.992, 0.990,
    0.032, 0.384, 0.710, 0.690, 0.747,
    
    # n=800, d=s=20：Z1, Z2, CCT, Ln
    0.049, 0.767, 0.921, 0.940, 0.950,
    0.039, 0.757, 0.935, 0.948, 0.950,
    0.043, 0.840, 0.960, 0.968, 0.968,
    0.060, 0.300, 0.513, 0.637, 0.670
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


# ==============================
# 颜色
# ==============================
my_cols <- c(
  "Z1" = "#3B82B8",
  "Z2" = "#E69F00",
  "Tc" = "#D55E5E",
  "Ln" = "#59A14F"
)

# ==============================
# 线型
# ==============================
my_linetypes <- c(
  "Z1" = "solid",
  "Z2" = "solid",
  "Tc" = "solid",
  "Ln" = "solid"
)

# ==============================
# 点的形状
# ==============================
my_shapes <- c(
  "Z1" = 16,   # 圆点
  "Z2" = 17,   # 三角形
  "Tc" = 15,   # 方块
  "Ln" = 8     # 星形
)

# 图例标签
my_labels <- c(
  expression(Z[n]^{(1)}),
  expression(Z[n]^{(2)}),
  expression(T[C]),
  expression(L[n])
)



# ==============================
# d=10 图
# ==============================
p1 <- ggplot(
  subset(df, scenario == "d10"),
  aes(
    x = c, y = value,
    color = Method,
    linetype = Method,
    shape = Method
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = my_cols,
    breaks = c("Z1", "Z2", "Tc", "Ln"),
    labels = my_labels
  ) + 
  scale_linetype_manual(
    values = my_linetypes,
    breaks = c("Z1", "Z2", "Tc", "Ln"),
    labels = my_labels
  ) +
  scale_shape_manual(
    values = my_shapes,
    breaks = c("Z1", "Z2", "Tc", "Ln"),
    labels = my_labels
  ) +
  # 5% 显著性水平
  geom_hline(yintercept = 0.05, linetype = "dashed") +
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
p2 <- ggplot(
  subset(df, scenario == "d20"),
  aes(
    x = c, y = value,
    color = Method,
    linetype = Method,
    shape = Method
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = my_cols,
    breaks = c("Z1", "Z2", "Tc", "Ln"),
    labels = my_labels
  ) + 
  scale_linetype_manual(
    values = my_linetypes,
    breaks = c("Z1", "Z2", "Tc", "Ln"),
    labels = my_labels
  ) +
  scale_shape_manual(
    values = my_shapes,
    breaks = c("Z1", "Z2", "Tc", "Ln"),
    labels = my_labels
  ) +
  # 5% 显著性水平
  geom_hline(yintercept = 0.05, linetype = "dashed") +
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

p_final <- (p1 + p2) + plot_layout(guides = "collect")

ggsave(
  filename = "plot/exp(-t_val2)/power_plot_exp(-t_val2).pdf",
  plot = p_final,
  width = 12,
  height = 5,
  units = "in",
  device = cairo_pdf   # 字体更清晰
)


