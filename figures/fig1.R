colors <- paletteer::paletteer_d("tvthemes::AirNomads")
colors[3] = "#FFD700"
colors = colors[c(2, 3, 5, 1)]

####--------------------------------------####
#### part (a) ####
####--------------------------------------####
## expression illustration
library(ggplot2)

set.seed(1)
n_row <- 10; n_col <- 15              # grid size
p_dark <- 0.35                        # fraction of darker tiles

df <- expand.grid(y = seq_len(n_row), x = seq_len(n_col))
df$val <- rbinom(nrow(df), 1, p_dark) # 0/1 for light/dark

p = ggplot(df, aes(x, y, fill = factor(val))) +
  geom_tile(color = "white", size = 0.5, width = 1, height = 1) +  # grid lines
  scale_fill_manual(values = c("0" = "#DEF5E5FF", "1" = "#22A884FF"), guide = "none") +
  coord_fixed(expand = FALSE) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_reverse(expand = c(0,0)) +                               # origin at top-left (like your image)
  theme_void() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

p

ggsave("./figures/output/fig1/pipelinea_expr.pdf", p, width = 7.56, height = 4.5, device = "pdf")

## coordinates 
library(ggplot2)
library(scatterpie)

set.seed(123)

## 1. Generate 20 Visium-like spot coordinates
# start from a grid and then keep points in an ellipse to get an irregular boundary
grid <- expand.grid(x = 1:8, y = 1:5)

mask <- ((grid$x - 5)^2 / 16 + (grid$y - 4)^2 / 9) <= 1.2
pool <- grid[mask, ]

# take 20 spots and jitter them a bit to avoid a perfect grid
n_spots <- 20
spots <- pool[sample(nrow(pool), n_spots), ]
spots$x <- spots$x
spots$y <- spots$y
spots$spot_id <- seq_len(n_spots)

## 2. Plot just the x, y coordinates
p = ggplot(spots, aes(x = x, y = y)) +
  geom_point(size = 8) +
  coord_equal() +
  labs(x = "X", y = "Y") +
  theme_minimal() + 
  theme(
    panel.grid = element_blank(), 
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.line = element_line(color = "black", size = 1)
  )
p

ggsave("./figures/output/fig1/pipelinea_coords.pdf", p, width = 7.56, height = 4.5, device = "pdf")

## 3. Assign 3 random proportions per spot that sum to 1
random_mat <- matrix(rexp(n_spots * 3), ncol = 3)
props <- random_mat / rowSums(random_mat)

spots$p1 <- props[, 1]
spots$p2 <- props[, 2]
spots$p3 <- props[, 3]

## 4. Plot pie charts at each (x, y) coordinate
library(paletteer)
## Pie chart plot with pretty colors
p = ggplot() +
  geom_scatterpie(
    data = spots,
    aes(x = x, y = y, r = 0.47, group = spot_id),
    cols = c("p1","p2","p3"), color = NA, alpha = 0.8
  ) +
  geom_circle(
    data = spots,
    aes(x0 = x, y0 = y, r = 0.47),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  coord_equal() +
  scale_fill_manual(
    values = colors, 
    ,labels = c(
        p1 = "Cell type 1",
        p2 = "Cell type 2",
        p3 = "Cell type 3"
      )
  )+
  # scale_fill_paletteer_d(
  #   "tvthemes::AirNomads"
  #   # "fishualize::Acanthurus_olivaceus"
  #   ,labels = c(
  #     p1 = "Cell type 1",
  #     p2 = "Cell type 2",
  #     p3 = "Cell type 3"
  #   )
  # ) +
  labs(x = "X", y = "Y", fill = "Component") +
  theme_minimal()+ 
  theme(
    panel.grid = element_blank(), 
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.line = element_line(color = "black", size = 1),
    # legend.title = element_blank(),
    legend.position = "none"
    # legend.direction = "horizontal",
    # legend.text = element_blank(),
    # legend.key = element_blank()
  )
p
ggsave("./figures/output/fig1/pipelinea_composition.pdf", p, width = 7.56, height = 4.5, device = "pdf")

####--------------------------------------####
#### part (b) cell type compositions ####
####--------------------------------------####
source("./simulators/idealized/celina_simulator_alt.R")
source("./simulators/idealized/celina_simulator_null.R")
source("./simulators/idealized/stance_simulator_1alt.R")
source("./simulators/idealized/stance_simulator_alt.R")
sim_setup = read.csv("./simulators/idealized/simulator_setup.csv")
seed = 42

# scenario 1
temp_sim <- run_simulator(sim_setup$sim_name[1], seed = seed, phi = sim_setup$phi[1], scene = sim_setup$scene[1], pattern = sim_setup$pattern[1], 
                          control_UMI = FALSE)
indx = sample(1:nrow(temp_sim$cell_metadata), 500)
plot_df = temp_sim$cell_metadata[indx, ]
p = ggplot(plot_df, aes(x = x, y = y, color = type)) +
  geom_point(size = 6, alpha = 0.8)+ 
  coord_equal() +
  scale_color_manual(
    values = colors
  )+
  # scale_color_paletteer_d(
  #   "tvthemes::AirNomads"
  #   # "fishualize::Acanthurus_olivaceus"
  # )+
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.title = element_blank(),
    legend.position = "none",
    panel.border     = element_rect(color = "black", fill = NA, size = 1)
    # legend.direction = "horizontal",
    # legend.text = element_blank(),
    # legend.key = element_blank()
  )
p
ggsave("./figures/output/fig1/pipelineb_scene1.pdf", p, width = 5, height = 5, device = "pdf")

# scenario 2
temp_sim <- run_simulator(sim_setup$sim_name[3], seed = seed, phi = sim_setup$phi[3], scene = sim_setup$scene[3], pattern = sim_setup$pattern[3], 
                          control_UMI = FALSE)
indx = sample(1:nrow(temp_sim$cell_metadata), 500)
plot_df = temp_sim$cell_metadata[indx, ]
types = unique(plot_df$type)
circle_df = data.frame(
  x0 = c(temp_sim$center[1, "x"], temp_sim$center[2, "x"]), 
  y0 = c(temp_sim$center[1, "y"], temp_sim$center[2, "y"]),
  radius = c(temp_sim$radius[1], rep(temp_sim$radius[2]))
)
p = ggplot(plot_df, aes(x = x, y = y, color = type)) +
  geom_point(size = 6, alpha = 0.6)+ 
  coord_equal() +
  geom_circle(data = circle_df, aes(x0=x0, y0=y0, r = radius), color = "black", size = 1, 
              inherit.aes = FALSE)+
  # scale_color_paletteer_d(
  #   # "tvthemes::AirNomads"
  #   "fishualize::Acanthurus_olivaceus"
  # )+
  scale_color_manual(
    values = colors
  )+
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.title = element_blank(),
    legend.position = "none",
    panel.border     = element_rect(color = "black", fill = NA, size = 1)
    # legend.direction = "horizontal",
    # legend.text = element_blank(),
    # legend.key = element_blank()
  )
p
ggsave("./figures/output/fig1/pipelineb_scene2.pdf", p, width = 5, height = 5, device = "pdf")

# scenario 3-6
scene = 8
temp_sim <- run_simulator(sim_setup$sim_name[scene], seed = seed, phi = sim_setup$phi[scene], scene = sim_setup$scene[scene], pattern = sim_setup$pattern[scene], 
                          control_UMI = FALSE)
indx = sample(1:nrow(temp_sim$cell_metadata), 500)
plot_df = temp_sim$cell_metadata[indx, ]
p = ggplot(plot_df, aes(x = x, y = y, color = type)) +
  geom_point(size = 6, alpha = 0.6)+ 
  coord_equal() +
  # scale_color_paletteer_d(
  #   # "tvthemes::AirNomads"
  #   "fishualize::Acanthurus_olivaceus"
  # )+
  scale_color_manual(
    values = colors
  )+
  geom_vline(xintercept = c(0.25, 0.5, 0.75), linetype = "solid", color = "black", size = 1) + 
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.title = element_blank(),
    legend.position = "none",
    panel.border     = element_rect(color = "black", fill = NA, size = 1)
    # legend.direction = "horizontal",
    # legend.text = element_blank(),
    # legend.key = element_blank()
  )
p
ggsave("./figures/output/fig1/pipelineb_scene6.pdf", p, width = 5, height = 5, device = "pdf")

## --- gene expression Hotspots 
source("./simulators/idealized/celina_simulator_alt.R")
source("./simulators/idealized/celina_simulator_null.R")
source("./simulators/idealized/stance_simulator_1alt.R")
source("./simulators/idealized/stance_simulator_alt.R")
scene = 1
seed = 1
temp_sim <- run_simulator(sim_setup$sim_name[scene], seed = seed, phi = sim_setup$phi[scene], scene = sim_setup$scene[scene], pattern = sim_setup$pattern[scene], 
                          control_UMI = FALSE)
circle_df = data.frame(
  x0 = c(temp_sim$center[1, "x"]), 
  y0 = c(temp_sim$center[1, "y"]),
  radius = c(temp_sim$radius[1])
)
plot_df = data.frame(
  # expr = temp_sim$cell_counts["Gene644",],
  expr = temp_sim$cell_counts["Gene701",],
  x = temp_sim$cell_metadata[,"x"],
  y = temp_sim$cell_metadata[,"y"],
  type = temp_sim$cell_metadata[,"type"]
)

selected_df = plot_df %>%
  filter(type %in% c("Type3"))
  # filter(type %in% c("Type1"))

cx <- as.numeric(temp_sim$center$x[1])
cy <- as.numeric(temp_sim$center$y[1])
r  <- as.numeric(temp_sim$radius[1])

p = ggplot() +
  # Grey background in all panels (drop 'type' so it’s drawn in every facet)
  geom_point(
    data = dplyr::select(plot_df, x, y),
    aes(x = x, y = y),
    color = "grey90", size = 3, inherit.aes = FALSE
  ) +
  # Colored points by cell type
  geom_point(
    data = selected_df,
    aes(x = x, y = y, color = expr),
    size = 3
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), ratio = 1, clip = "off") +  # <- one coord
  scale_color_viridis_c(trans =  scales::pseudo_log_trans(sigma = 0.5)) +
  # annotate("segment", x = cx, xend = cx, y = cy - r, yend = cy + r,
  #            colour = "#ED2939", linetype = "dashed", size = 1.5)+
  geom_circle(data = circle_df, aes(x0=x0, y0=y0, r = radius), color = "#ED2939", linetype = "dashed", size = 1.5, 
              inherit.aes = FALSE)+
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.title = element_blank(),
    legend.position = "none",
    panel.border     = element_rect(color = "black", fill = NA, size = 1)
    # legend.direction = "horizontal",
    # legend.text = element_blank(),
    # legend.key = element_blank()
  )
p
ggsave("./figures/output/fig1/pipelineb_hotspot1.pdf", p, width = 5, height = 5, device = "pdf")


## --- gene expression Streaks 
temp_sim <- run_simulator(sim_setup$sim_name[10], seed = seed, phi = sim_setup$phi[10], scene = sim_setup$scene[10], pattern = sim_setup$pattern[10], 
                          control_UMI = FALSE)

plot_df = data.frame(
  expr = temp_sim$cell_counts["Gene644",],
  # expr = scales::rescale(temp_sim$cell_counts["Gene701",]), 
  x = temp_sim$cell_metadata[,"x"],
  y = temp_sim$cell_metadata[,"y"],
  type = temp_sim$cell_metadata[,"type"]
)

selected_df = plot_df %>%
  filter(type %in% c("Type2")) 
  
p = ggplot() +
  # Grey background in all panels (drop 'type' so it’s drawn in every facet)
  geom_point(
    data = dplyr::select(plot_df, x, y),
    aes(x = x, y = y),
    color = "grey90", size = 3, inherit.aes = FALSE
  ) +
  # Colored points by cell type
  geom_point(
    data = selected_df,
    aes(x = x, y = y, color = expr),
    size = 3
  ) +
  scale_color_viridis_c(trans =  scales::pseudo_log_trans(sigma = 1.5))+
  annotate("rect", xmin = 0, xmax = 1, ymin = temp_sim$radius, ymax = 1,
           fill = NA, color = "#ED2939", size = 1.5)+
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.title = element_blank(),
    legend.position = "none",
    panel.border     = element_rect(color = "black", fill = NA, size = 1)
    # legend.direction = "horizontal",
    # legend.text = element_blank(),
    # legend.key = element_blank()
  )
p
ggsave("./figures/output/fig1/pipelineb_streak.pdf", p, width = 5, height = 5, device = "pdf")

## --- gene expression gradient 
source("./simulators/idealized/celina_simulator_alt.R")
source("./simulators/idealized/celina_simulator_null.R")
source("./simulators/idealized/stance_simulator_1alt.R")
source("./simulators/idealized/stance_simulator_alt.R")
temp_sim <- run_simulator(sim_setup$sim_name[16], seed = seed, phi = sim_setup$phi[16], scene = sim_setup$scene[16], pattern = sim_setup$pattern[16], 
                          control_UMI = FALSE)

plot_df = data.frame(
  expr = temp_sim$cell_counts["Gene644",],
  # expr = scales::rescale(temp_sim$cell_counts["Gene701",]), 
  x = temp_sim$cell_metadata[,"x"],
  y = temp_sim$cell_metadata[,"y"],
  type = temp_sim$cell_metadata[,"type"]
)

selected_df = plot_df %>%
  filter(type %in% c("Type3")) 


p = ggplot() +
  # Grey background in all panels (drop 'type' so it’s drawn in every facet)
  geom_point(
    data = dplyr::select(plot_df, x, y),
    aes(x = x, y = y),
    color = "grey90", size = 3, inherit.aes = FALSE
  ) +
  # Colored points by cell type
  geom_point(
    data = selected_df,
    aes(x = x, y = y, color = expr),
    size = 3
  ) +
  scale_color_viridis_c( trans =  scales::pseudo_log_trans(sigma = 1))+
  annotate(
    "segment",
    x = 0.6, xend = 0.6,
    y = 0.1, yend = 0.9,
    colour = "#ED2939", linetype = "dashed", size = 1.5,
    arrow = arrow(length = unit(0.08, "npc"), type = "closed")
  )+
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.title = element_blank(),
    legend.position = "none",
    panel.border     = element_rect(color = "black", fill = NA, size = 1)
    # legend.direction = "horizontal",
    # legend.text = element_blank(),
    # legend.key = element_blank()
  )
p
ggsave("./figures/output/fig1/pipelineb_gradient.pdf", p, width = 5, height = 5, device = "pdf")



####--------------------------------------####
#### part (c) gene expression ####
####--------------------------------------####

## at cell level 
source("./simulators/realistic/pseudo_spot_simulator.R")
load("./data/realistic_packaged/ovarian/data_new_v3_42.RData")

H19_df = data.frame(
  expr = data_new$cell_counts_predrop["H19", ], 
  x = data_new$cell_metadata[, "x"], 
  y = data_new$cell_metadata[, "y"], 
  type = data_new$cell_metadata[, "type"]
)

library(ggh4x)

selected_df = H19_df %>%
  filter(type %in% c("Tumor Cells", "Macrophages", "VEGFA+ Tumor Cells", "Tumor Associated Fibroblasts", "T and NK Cells")) %>%
  mutate(
    expr = expr,
    type = factor(
      type, 
      levels = c("Tumor Cells", "Macrophages", "VEGFA+ Tumor Cells", "Tumor Associated Fibroblasts", "T and NK Cells"),
      labels = c("H19 (Tumor)", "H19 (Macrophage)", "H19 (VEGFA Tumor)", "H19 (Fibroblast)", "H19 (T&NK)")
    )
  ) 

p = ggplot() +
  # Grey background in all panels (drop 'type' so it’s drawn in every facet)
  geom_point(
    data = dplyr::select(H19_df, x, y),
    aes(x = x, y = y),
    color = "grey90", size = 3, inherit.aes = FALSE
  ) +
  # Colored points by cell type
  geom_point(
    data = selected_df,
    aes(x = x, y = y, color = expr),
    size = 3
  ) +
  ggh4x::facet_wrap2(vars(type), nrow = 1) +
  # scale_color_viridis_c(option = "plasma", trans =  scales::pseudo_log_trans(sigma = 1.5))+
  scale_color_viridis_c(trans =  scales::pseudo_log_trans(sigma = 1.5))+
  theme(
    # strip.text = element_text(size = 26, margin = margin(t = 10, b = 10)),
    strip.text = element_blank(),
    strip.background = element_rect(fill = "white", color = NA, linewidth = 1),
    panel.grid = element_blank(), axis.text = element_blank(),
    axis.ticks = element_blank(), axis.title = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = NA, fill = NA, linewidth = 1),
    legend.title = element_blank(), legend.text = element_blank()
  )

p

ggsave("./figures/output/fig1/pipelinec_predrop.pdf",  p, height = 4.22, width = 19.06, device = "pdf")

## at spot level 
H19_df = data.frame(
  expr = scales::rescale(data_new$spot_counts["H19", ]), 
  x = data_new$spot_coords[, "x"], 
  y = data_new$spot_coords[, "y"], 
  label = "H19"
  ) %>%
  filter( x> 10000 & x < 16000 & y > 16000 & y < 22000
    
  )

library(ggh4x)
p = ggplot()+
geom_point(
  data = dplyr::select(H19_df, x, y, expr, label),
  aes(x = x, y = y, color = expr),
  size = 6, inherit.aes = FALSE
) +
scale_color_viridis_c(trans =  scales::pseudo_log_trans(sigma = 0.1))+
theme(
  panel.grid = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.title = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  legend.title = element_blank(), legend.text = element_blank(), 
  legend.position = "none"
)
p
ggsave("./figures/output/fig1/pipelinec_spot.pdf",  p, height = 5, width = 5, device = "pdf")

####--------------------------------------####
#### part (c) pie chart ####
####--------------------------------------####
df <- data.frame(
  category = c("A", "B", "C"),
  value    = c(60, 20, 20)
)
colors <- paletteer::paletteer_d("tvthemes::AirNomads")
pie = ggplot(df, aes(x = "", y = value, fill = category)) +
  geom_col(width = 1, color = "black", linewidth = 1, alpha = 0.5) +                # stacked bar
  coord_polar(theta = "y") +           # wrap into a pie
  theme_void() +                       # remove axes, background
  labs(fill = "Category") + 
  theme(
    legend.position = "none"
  )+ 
  scale_fill_manual(
    values = c(colors[3], colors[4], colors[5])
  )
pie

ggsave("./figures/output/fig1/pipelinc_pie.pdf", pie, width = 4, height = 4, device = "pdf")

####--------------------------------------####
#### part (d) gene expression ####
####--------------------------------------####



####--------------------------------------####
#### part (d) cell type composition ####
####--------------------------------------####
temp_sim <- run_simulator(sim_setup$sim_name[1], seed = 42, phi = sim_setup$phi[1], scene = sim_setup$scene[1], pattern = sim_setup$pattern[1], 
                          control_UMI = FALSE)
colors <- paletteer::paletteer_d("tvthemes::AirNomads")
colors[3] = "#FFD700"
colors = colors[c(2, 3, 5, 1)]

indx = sample(1:nrow(temp_sim$cell_metadata), 500)
plot_df = temp_sim$cell_metadata[indx, ]
p = ggplot(plot_df, aes(x = x, y = y, color = type)) +
  geom_point(size = 6, alpha = 0.8)+ 
  coord_equal() +
  scale_color_manual(
    values = colors
  )+
  # scale_color_paletteer_d(
  #   "tvthemes::AirNomads"
  #   # "fishualize::Acanthurus_olivaceus"
  # )+
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.title = element_text(size = 30),
    legend.position = "right",
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    legend.direction = "vertical",
    # legend.text = element_blank(),
    legend.text = element_blank(), 
    legend.title = element_blank(), 
    legend.key.size = unit(1.5, "lines")
  )
ggsave("./figures/output/fig1/pipelined_tissue1.pdf", p, width = 5, height = 5, device = "pdf")



library(ggplot2)
library(arrow)
library(Seurat)
library(Matrix)
library(Triangulation)
library(readr)
source("./util.R")
source("./simulators/realistic/pseudo_spot_simulator.R")
source("./simulators/decomposed/scDesign_decompose.R")

load(paste0("./data/realistic_packaged/ovarian/ovarian_small_ori.RData"))
load(paste0("./data/realistic_packaged/ovarian/data_new_v3_42.RData"))
colors <- paletteer::paletteer_d("tvthemes::AirNomads")
colors[3] = "#FFD700"
temp_sim = decompose_simulator$new(data_new, data_ori = ovarian_small_ori, dispersion = 0.7, seed = 42, sim_type = c(1), n_markers = 25)

indx = sample(1:nrow(temp_sim$cell_metadata), 500)
plot_df = temp_sim$cell_metadata[indx, ]
colnames(plot_df)[3] = "Type"
top = names(sort(table(plot_df$Type)/length(plot_df$Type), decreasing = TRUE))[1:10]
selected_type = mapply(plot_df$Type, FUN = function(x) ifelse(x %in% top, x, "Other Cell Types"))
plot_df$selected_type = selected_type
plot_df$selected_type = factor(selected_type, levels = top)
p = ggplot(plot_df, aes(x = x, y = y, color = selected_type)) +
  geom_point(size = 6, alpha = 0.6)+ 
  coord_equal() +
  # scale_fill_brewer(
  #   palette = "Set2"
  # )+
  scale_color_manual(
    values = c(colors[1:7],"royalblue", "plum1", "seagreen2", "grey80")
    )+
  theme(
    aspect.ratio   = 1,
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.title = element_text(size = 30),
    legend.position = "right",
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    legend.direction = "vertical",
    # legend.text = element_blank(),
    legend.text = element_blank(), 
    legend.title = element_blank(), 
    legend.key.size = unit(1.5, "lines")
  )
p
ggsave("./figures/output/fig1/pipelined_tissue2.pdf", p, width = 5, height = 5, device = "pdf")


### location 
temp_sim = decompose_simulator$new(data_new, data_ori = ovarian_small_ori, dispersion = 0.7, seed = 42, sim_type = c(2), n_markers = 25)
plot(temp_sim$cell_metadata$x, temp_sim$cell_metadata$y)
# indx_dense = which(temp_sim$cell_metadata$x <=13000 & temp_sim$cell_metadata$x >=10500 & temp_sim$cell_metadata$y <=15000 & temp_sim$cell_metadata$y >=12000)
# indx_loose = which(temp_sim$cell_metadata$x <=16000, temp_sim$cell_metadata$y <=13000)
# indx_loose = setdiff(indx_loose, indx_dense)
indx = sample(1:nrow(temp_sim$cell_metadata), 500)
plot_df = temp_sim$cell_metadata[indx, ]
colors <- paletteer::paletteer_d("tvthemes::AirNomads")
colors[3] = "#FFD700"
colors = colors[c(2, 3, 5, 1)]
# indx = sample(indx_dense, 350)
# indx = c(indx, sample(indx_loose, 100))
# plot_df = temp_sim$cell_metadata[indx, ]
colnames(plot_df)[3] = "Type"
p = ggplot(plot_df, aes(x = x, y = y, color = Type)) +
  geom_point(size = 6, alpha = 0.8)+ 
  coord_fixed(ratio = 1)+
  scale_color_manual(
    values = colors
  )+
  theme(
    aspect.ratio   = 1,
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.title = element_text(size = 30),
    legend.position = "right",
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    legend.direction = "vertical",
    # legend.text = element_blank(),
    legend.text = element_blank(), 
    legend.title = element_blank(), 
    legend.key.size = unit(1.5, "lines")
  )
p
ggsave("./figures/output/fig1/pipelined_tissue3.pdf", p, width = 5, height = 5, device = "pdf")


### cell type distribution 
temp_sim = decompose_simulator$new(data_new, data_ori = ovarian_small_ori, dispersion = 0.7, seed = 42, sim_type = c(1, 2), n_markers = 25)
colors <- paletteer::paletteer_d("tvthemes::AirNomads")
plot_df = temp_sim$cell_metadata[indx, ]
colnames(plot_df)[3] = "Type"
selected_type = mapply(plot_df$Type, FUN = function(x) ifelse(x %in% top, x, "Other Cell Types"))
plot_df$selected_type = selected_type
plot_df$selected_type = factor(selected_type, levels = top)
p = ggplot(plot_df, aes(x = x, y = y, color = selected_type)) +
  geom_point(size = 6, alpha = 0.6)+ 
  coord_fixed(ratio = 1)+
  scale_color_manual(
    values = c(colors[1:7],"royalblue", "plum1", "seagreen2", "grey80")
  )+
  theme(
    aspect.ratio   = 1,
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.title = element_text(size = 30),
    legend.position = "right",
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    legend.direction = "vertical",
    # legend.text = element_blank(),
    legend.text = element_blank(), 
    legend.title = element_blank(), 
    legend.key.size = unit(1.5, "lines")
  )
p
ggsave("./figures/output/fig1/pipelined_tissue4.pdf", p, width = 5, height = 5, device = "pdf")


####--------------------------------------####
#### part (d) gene expression ####
####--------------------------------------####
### 1. pure simulation 
source("./simulators/idealized/celina_simulator_alt.R")
source("./simulators/idealized/celina_simulator_null.R")
source("./simulators/idealized/stance_simulator_1alt.R")
source("./simulators/idealized/stance_simulator_alt.R")
scene = 1
seed = 19 #1
temp_sim <- run_simulator(sim_setup$sim_name[scene], seed = seed, phi = sim_setup$phi[scene], scene = sim_setup$scene[scene], pattern = sim_setup$pattern[scene], 
                          control_UMI = FALSE)

plot_df = data.frame(
  expr = temp_sim$cell_counts["Gene501",],
  x = temp_sim$cell_metadata[,"x"],
  y = temp_sim$cell_metadata[,"y"],
  type = temp_sim$cell_metadata[,"type"]
)

selected_df = plot_df %>%
  filter(type %in% c("Type2"))

cx <- as.numeric(temp_sim$center$x[1])
cy <- as.numeric(temp_sim$center$y[1])
r  <- as.numeric(temp_sim$radius[1])

p = ggplot() +
  # Grey background in all panels (drop 'type' so it’s drawn in every facet)
  geom_point(
    data = dplyr::select(plot_df, x, y),
    aes(x = x, y = y),
    color = "grey90", size = 3, inherit.aes = FALSE
  ) +
  # Colored points by cell type
  geom_point(
    data = selected_df,
    aes(x = x, y = y, color = expr),
    size = 3
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), ratio = 1, clip = "off") +  # <- one coord
  scale_color_viridis_c(
    trans =  scales::pseudo_log_trans(sigma = 0.5) 
    # ,rescaler = function(x, ...) scales::rescale(x, to = c(0.3, 0.6))
  ) +
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.title = element_blank(),
    legend.position = "none",
    panel.border     = element_rect(color = "black", fill = NA, size = 1)
    # legend.direction = "horizontal",
    # legend.text = element_blank(),
    # legend.key = element_blank()
  )
p
ggsave("./figures/output/fig1/pipelined_expression1.2.pdf", p, width = 5, height = 5, device = "pdf")


### realistic ctSVGs

source("./simulators/realistic/pseudo_spot_simulator.R")
load("./data/realistic_packaged/ovarian/data_new_v3_42.RData")

H19_df = data.frame(
  expr = data_new$cell_counts_predrop["H19", ], 
  x = data_new$cell_metadata[, "x"], 
  y = data_new$cell_metadata[, "y"], 
  type = data_new$cell_metadata[, "type"]
)

library(ggh4x)
selected_df = H19_df %>%
  # filter(type %in% c("Tumor Cells"))
  # filter(type %in% c("Macrophages"))
  # filter(type %in% c("Tumor Associated Fibroblasts")) 
  filter(type %in% c("T and NK Cells"))
p = ggplot() +
  # Grey background in all panels (drop 'type' so it’s drawn in every facet)
  geom_point(
    data = dplyr::select(H19_df, x, y),
    aes(x = x, y = y),
    color = "grey90", size = 3, inherit.aes = FALSE
  ) +
  # Colored points by cell type
  geom_point(
    data = selected_df,
    aes(x = x, y = y, color = expr),
    size = 3
  ) +
  scale_color_viridis_c(trans =  scales::pseudo_log_trans(sigma = 0.5))+
  theme_minimal()+
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.title = element_blank(),
    legend.position = "none",
    strip.text = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, size = 1)
    # legend.direction = "horizontal",
    # legend.text = element_blank(),
    # legend.key = element_blank()
  )
p
ggsave("./figures/output/fig1/pipelined_expression4.4.pdf", p, width = 5, height = 5, device = "pdf")


### pseudo spot capture adjustment 
source("./simulators/decomposed/scDesign_decompose.R")
temp_sim = decompose_simulator$new(data_new, data_ori = ovarian_small_ori, dispersion = 0.7, seed = 1, sim_type = c(2, 4), n_markers = 25)

# 2.3-2.4 ctsvg
plot_df = data.frame(
  expr = temp_sim$cell_counts["F13A1", ],
  x = temp_sim$cell_metadata[, "x"],
  y = temp_sim$cell_metadata[, "y"],
  type = temp_sim$cell_metadata[, "type"]
)
selected_df = plot_df %>%
  filter(type %in% c("Type1"))

circle_df = data.frame(
  x0 = c(temp_sim$center[1, "x"]), 
  y0 = c(temp_sim$center[1, "y"]),
  radius = c(temp_sim$radius[1])
)
p = ggplot() +
  # Grey background in all panels (drop 'type' so it’s drawn in every facet)
  geom_point(
    data = dplyr::select(plot_df, x, y),
    aes(x = x, y = y),
    color = "grey90", size = 3, inherit.aes = FALSE
  ) +
  # Colored points by cell type
  geom_point(
    data = selected_df,
    aes(x = x, y = y, color = expr),
    size = 3
  ) +
  scale_color_viridis_c(
    trans =  scales::pseudo_log_trans(sigma = 0.5)
    ,rescaler = function(x, ...) scales::rescale(x, to = c(0.0, 1))
  ) +
  # geom_circle(data = circle_df, aes(x0=x0, y0=y0, r = radius), color = "#ED2939", linetype = "dashed", size = 1.5,
  #             inherit.aes = FALSE)+
  theme_minimal()+
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.title = element_blank(),
    legend.position = "none",
    strip.text = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, size = 1)
    # legend.direction = "horizontal",
    # legend.text = element_blank(),
    # legend.key = element_blank()
  )
p
ggsave("./figures/output/fig1/pipelined_expression3.2.pdf", p, width = 5, height = 5, device = "pdf")



# grid <- tidyr::expand_grid(row = 1:6, col = 1:6) |>
#   mutate(value = 0)
# # 1.1 ctsvg 
# expr = c(
#   1, 1, 1, 1, 1, 1,
#   1, 1, 1, 1, 1, 1,
#   1, 1, 1, 1, 1, 1, 
#   6, 10, 6, 1, 1, 1, 
#   10, 10, 10, 1, 1, 1, 
#   5, 10, 7, 1, 1, 1
# )
# 
# expr = c(
#   1, 1, 1, 1, 1, 1,
#   1, 1, 1, 1, 1, 1,
#   1, 1, 1, 6, 10, 6,  
#   1, 1, 1, 10, 10, 10,
#   1, 1, 1, 5, 10, 7,
#   1, 1, 1, 1, 1, 1
# )
# 
# expr = c(
#   1, 1, 1, 1, 1, 1,
#   1, 1, 1, 1, 1, 1,
#   1, 1, 1, 1, 1, 1, 
#   1, 1, 1, 1, 1, 1,
#   1, 1, 1, 1, 1, 1,
#   1, 1, 1, 1, 1, 1
# )
# 
# expr = c(
#   6, 8, 8, 7, 8, 14,
#   1, 1, 2, 12, 7, 20,
#   1, 1, 1, 6, 10, 6,  
#   1, 1, 1, 1, 10, 10,
#   1, 1, 1, 1, 10, 7,
#   1, 1, 1, 1, 1, 1
# )
# 
# expr = c(
#   1, 10, 10, 8, 1, 1,
#   10, 1, 1, 10, 5, 1,
#   10, 1, 1, 10, 5, 1, 
#   10, 1, 1, 10, 5, 1,
#   8, 1, 1, 8, 4, 1,
#   1, 10, 10, 3, 1, 1
# )
# 
# expr = c(
#   1, 1, 1, 1, 1, 1, 
#   1, 1, 1, 1, 1, 1, 
#   1, 1, 15, 10, 4, 1,  
#   1, 1, 20, 10, 5, 1,  
#   1, 1, 15, 10, 7, 1, 
#   1, 1, 1, 1, 1, 1
# )
# 
# grid$expr = expr + rnorm(36, 0, 1)
# 
# p = ggplot() +
#   geom_tile(
#     data = grid,
#     aes(x = row, y = col, fill = expr), 
#     width = 1, height = 1
#   ) +
#   coord_cartesian(x = c(0.5, 6.5), y = c(0.5, 6.5))+
#   # scale_color_viridis_c(trans =  scales::pseudo_log_trans(sigma = 1))+
#   scale_fill_viridis_c(
#     option = "viridis",
#     trans = scales::pseudo_log_trans(sigma = 1),
#     rescaler = function(x, ...) scales::rescale(x, to = c(0, 1))
#   )+
#   theme_minimal()+
#   theme(
#     panel.grid = element_blank(),
#     panel.background = element_blank(),
#     axis.text = element_blank(),
#     axis.title = element_blank(),
#     axis.line = element_blank(),
#     axis.ticks = element_blank(),
#     # legend.title = element_blank(),
#     legend.position = "none",
#     strip.text = element_blank(),
#     panel.border     = element_rect(color = "black", fill = NA, size = 6)
#     # legend.direction = "horizontal",
#     # legend.text = element_blank(),
#     # legend.key = element_blank()
#   )
# p
# 
# ggsave("./figures/output/fig1/pipelined_ctsvg1.3.pdf", p, width = 5, height = 5, device = "pdf")

# temp_sim = decompose_simulator$new(data_new, data_ori = ovarian_small_ori, dispersion = 0.7, seed = 42, sim_type = c(0), n_markers = 25)
# 
# # 1.1 ctsvg 
# plot_df = data.frame(
#   expr = temp_sim$cell_counts["F13A1", ],
#   x = temp_sim$cell_metadata[, "x"],
#   y = temp_sim$cell_metadata[, "y"],
#   type = temp_sim$cell_metadata[, "type"]
# )
# selected_df = plot_df %>%
#   filter(type %in% c("Type1")) 
# 
# # 1.2 ctsvg 
# plot_df = data.frame(
#   expr = temp_sim$cell_counts["SFRP4", ],
#   x = temp_sim$cell_metadata[, "x"],
#   y = temp_sim$cell_metadata[, "y"],
#   type = temp_sim$cell_metadata[, "type"]
# )
# selected_df = plot_df %>%
#   filter(type %in% c("Type3")) 
# 
# # 1.3 null 
# plot_df = data.frame(
#   expr = temp_sim$cell_counts["F13A1", ],
#   x = temp_sim$cell_metadata[, "x"],
#   y = temp_sim$cell_metadata[, "y"],
#   type = temp_sim$cell_metadata[, "type"]
# )
# selected_df = plot_df %>%
#   filter(type %in% c("Type2")) 
# 
# # 1.4 null 
# plot_df = data.frame(
#   expr = temp_sim$cell_counts["SFRP4", ],
#   x = temp_sim$cell_metadata[, "x"],
#   y = temp_sim$cell_metadata[, "y"],
#   type = temp_sim$cell_metadata[, "type"]
# )
# selected_df = plot_df %>%
#   filter(type %in% c("Type1")) 
# 
# 
# # 1.4 null 
# 
# 
# p = ggplot() +
#   # Grey background in all panels (drop 'type' so it’s drawn in every facet)
#   geom_point(
#     data = dplyr::select(plot_df, x, y),
#     aes(x = x, y = y),
#     color = "grey80", size = 3, inherit.aes = FALSE
#   ) +
#   # Colored points by cell type
#   geom_point(
#     data = selected_df,
#     aes(x = x, y = y, color = expr),
#     size = 3
#   ) +
#   scale_color_viridis_c(trans =  scales::pseudo_log_trans(sigma = 0.5))+
#   theme_minimal()+
#   theme(
#     panel.grid = element_blank(), 
#     panel.background = element_blank(),
#     axis.text = element_blank(), 
#     axis.title = element_blank(), 
#     axis.line = element_blank(),
#     axis.ticks = element_blank(),
#     # legend.title = element_blank(),
#     legend.position = "none",
#     strip.text = element_blank(),
#     panel.border     = element_rect(color = "black", fill = NA, size = 6)
#     # legend.direction = "horizontal",
#     # legend.text = element_blank(),
#     # legend.key = element_blank()
#   )
# p
# ggsave("./figures/output/fig1/pipelined_ctsvg1.4.pdf", p, width = 5, height = 5, device = "pdf")
# 
# # pure simulation
# temp_sim = decompose_simulator$new(data_new, data_ori = ovarian_small_ori, dispersion = 0.7, seed = 42, sim_type = c(3), n_markers = 25)
# 
# # 2.3 null 
# plot_df = data.frame(
#   expr = temp_sim$cell_counts["F13A1", ],
#   x = temp_sim$cell_metadata[, "x"],
#   y = temp_sim$cell_metadata[, "y"],
#   type = temp_sim$cell_metadata[, "type"]
# )
# selected_df = plot_df %>%
#   filter(type %in% c("Type2")) 
# 
# # 2.4 null 
# plot_df = data.frame(
#   expr = temp_sim$cell_counts["AQP8", ],
#   x = temp_sim$cell_metadata[, "x"],
#   y = temp_sim$cell_metadata[, "y"],
#   type = temp_sim$cell_metadata[, "type"]
# )
# selected_df = plot_df %>%
#   filter(type %in% c("Type1")) 
# 
# 
# p = ggplot() +
#   # Grey background in all panels (drop 'type' so it’s drawn in every facet)
#   geom_point(
#     data = dplyr::select(plot_df, x, y),
#     aes(x = x, y = y),
#     color = "grey80", size = 3, inherit.aes = FALSE
#   ) +
#   # Colored points by cell type
#   geom_point(
#     data = selected_df,
#     aes(x = x, y = y, color = expr),
#     size = 3
#   ) +
#   scale_color_viridis_c(trans =  scales::pseudo_log_trans(sigma = 3))+
#   theme_minimal()+
#   theme(
#     panel.grid = element_blank(), 
#     panel.background = element_blank(),
#     axis.text = element_blank(), 
#     axis.title = element_blank(), 
#     axis.line = element_blank(),
#     axis.ticks = element_blank(),
#     # legend.title = element_blank(),
#     legend.position = "none",
#     strip.text = element_blank(),
#     panel.border     = element_rect(color = "black", fill = NA, size = 6)
#     # legend.direction = "horizontal",
#     # legend.text = element_blank(),
#     # legend.key = element_blank()
#   )
# p
# ggsave("./figures/output/fig1/pipelined_ctsvg2.3.pdf", p, width = 5, height = 5, device = "pdf")
