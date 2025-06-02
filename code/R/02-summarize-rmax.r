# Nicholas Ducharme-Barth
# 2025/06/01
# Summarize Rmax prior - Updated for current analysis

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
library(data.table)
library(magrittr)
library(ggplot2)
library(GGally)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
proj_dir = this.path::this.proj()
plot_dir = file.path(proj_dir, "plots", "rmax-priors")
dir.create(plot_dir, recursive = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load data
rmax_dt = readRDS(file.path(proj_dir, "data", "output", "bspm_parameter_priors_parallel.rds"))

# Filter out failed runs for plotting
rmax_clean = rmax_dt[!is.na(rmax) & !is.na(F_est) & !is.na(spr)]
cat("Plotting", nrow(rmax_clean), "successful runs out of", nrow(rmax_dt), "total\n")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# 1. Basic Rmax density plot
p = rmax_clean %>% 
    ggplot() +
    ylim(0, NA) +
    xlab("Rmax") +
    ylab("Density") +
    geom_density(aes(x = rmax), fill = "blue", alpha = 0.25) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    viridis::scale_color_viridis("Model", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    viridis::scale_fill_viridis("Model", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"))

ggsave(filename = "rmax_prior.dens.all.png", plot = p, device = "png", path = plot_dir,
       scale = 1, width = 6, height = 6, units = "in", dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# 2. Input parameter distributions
p = rmax_clean %>% 
    .[, .(sample_id, max_age, M_ref, L1, L2, vbk, l50, sex_ratio, cv_len, 
          selex_l50, selex_slope)] %>%
    melt(., id.vars = "sample_id") %>%
    ggplot() +
    facet_wrap(~variable, scales = "free") +
    ylim(0, NA) +
    xlab("Input") +
    ylab("Count") +
    geom_histogram(aes(x = value, fill = variable), bins = 100, alpha = 0.5) +
    geom_hline(yintercept = 0) +
    viridis::scale_color_viridis("Input\nvariable", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    viridis::scale_fill_viridis("Input\nvariable", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"))

ggsave(filename = "rmax_prior.main.inputs.png", plot = p, device = "png", path = plot_dir,
       scale = 1, width = 12, height = 9, units = "in", dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# 3. Rmax density comparing all vs viable populations
p1 = rmax_clean[, .(rmax)] %>%
    .[, distribution := "all"]
p2 = rmax_clean[rmax > 0, .(rmax)] %>%
    .[, distribution := "survival"]

p = rbind(p1, p2) %>%
    ggplot() +
    ylim(0, NA) +
    xlab("Rmax") +
    ylab("Density") +
    geom_density(aes(x = rmax, fill = distribution), alpha = 0.25) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    viridis::scale_color_viridis("Distribution", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    viridis::scale_fill_viridis("Distribution", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"))

ggsave(filename = "rmax_prior.dens.filter.png", plot = p, device = "png", path = plot_dir,
       scale = 1, width = 6, height = 6, units = "in", dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# 4. Input parameters comparing all vs viable populations
p1 = rmax_clean[, .(sample_id, max_age, M_ref, L1, L2, vbk, l50, sex_ratio, cv_len)] %>%
    .[, distribution := "all"] %>%
    melt(., id.vars = c("sample_id", "distribution")) 

p2 = rmax_clean[rmax > 0] %>%
    .[, .(sample_id, max_age, M_ref, L1, L2, vbk, l50, sex_ratio, cv_len)] %>%
    .[, distribution := "survival"] %>%
    melt(., id.vars = c("sample_id", "distribution"))

p = rbind(p1, p2) %>% 
    ggplot() +
    facet_wrap(~variable, scales = "free") +
    ylim(0, NA) +
    xlab("Input") +
    ylab("Relative count") +
    geom_histogram(aes(x = value, y = after_stat(ncount), fill = distribution), 
                   position = "dodge", bins = 100, alpha = 0.5) +
    geom_hline(yintercept = 0) +
    viridis::scale_color_viridis("Distribution", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    viridis::scale_fill_viridis("Distribution", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"))

ggsave(filename = "rmax_prior.main.inputs.filter.png", plot = p, device = "png", path = plot_dir,
       scale = 1, width = 12, height = 9, units = "in", dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# 5. Pairs plot - biological inputs with outcomes
set.seed(123)
sample_size = min(10000, nrow(rmax_clean))

p = rmax_clean[, .(rmax, max_age, M_ref, L1, L2, vbk, l50, sex_ratio, cv_len)] %>%
    .[sample(1:.N, sample_size)] %>%
    .[, Outcome := "Extinction"] %>%
    .[rmax > 0, Outcome := "Survival"] %>%
    .[, Outcome := as.factor(Outcome)] %>%
    ggpairs(., columns = 1:9, aes(color = Outcome, alpha = 0.4)) +
    viridis::scale_color_viridis("Outcome", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    viridis::scale_fill_viridis("Outcome", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"))

ggsave(filename = "rmax_prior.pairs.continuous_inputs.outcomes.png", plot = p, device = "png", path = plot_dir,
       scale = 1, width = 16, height = 12, units = "in", dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# 6. Pairs plot - weight and selectivity parameters
p = rmax_clean[, .(rmax, weight_a, weight_b, selex_l50, selex_slope, selexNZ_l50, selexNZ_slope)] %>%
    .[sample(1:.N, sample_size)] %>%
    .[, Outcome := "Extinction"] %>%
    .[rmax > 0, Outcome := "Survival"] %>%
    .[, Outcome := as.factor(Outcome)] %>%
    .[, weight_a := log(weight_a)] %>%  # Log transform for better visualization
    ggpairs(., columns = 1:7, aes(color = Outcome, alpha = 0.4)) +
    viridis::scale_color_viridis("Outcome", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    viridis::scale_fill_viridis("Outcome", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"))

ggsave(filename = "rmax_prior.pairs.weight_selex_inputs.outcomes.png", plot = p, device = "png", path = plot_dir,
       scale = 1, width = 16, height = 12, units = "in", dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# 7. Pairs plot - derived quantities with outcomes
p = rmax_clean[, .(rmax, spr, alpha, h, generation_time, F_est)] %>%
    .[sample(1:.N, sample_size)] %>%
    .[, Outcome := "Extinction"] %>%
    .[rmax > 0, Outcome := "Survival"] %>%
    .[, Outcome := as.factor(Outcome)] %>%
    ggpairs(., columns = 1:6, aes(color = Outcome, alpha = 0.4)) +
    viridis::scale_color_viridis("Outcome", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    viridis::scale_fill_viridis("Outcome", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"))

ggsave(filename = "rmax_prior.pairs.derived_quants.outcomes.png", plot = p, device = "png", path = plot_dir,
       scale = 1, width = 16, height = 12, units = "in", dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# 8. Pairs plot - derived quantities for viable populations only
p = rmax_clean[, .(rmax, spr, alpha, h, generation_time, inflection_point)] %>%
    .[sample(1:.N, sample_size)] %>%
    .[rmax > 0 & !is.na(inflection_point)] %>%
    ggpairs(., columns = 1:6, aes(color = "blue", alpha = 0.4)) +
    viridis::scale_color_viridis("Outcome", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    viridis::scale_fill_viridis("Outcome", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"))

ggsave(filename = "rmax_prior.pairs.derived_quants.survival.png", plot = p, device = "png", path = plot_dir,
       scale = 1, width = 16, height = 12, units = "in", dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# 9. Additional plot - Depletion metrics
p = rmax_clean[, .(rmax,inflection_point,dep, dep_sb, spr, F_est)] %>%
    .[sample(1:.N, sample_size)] %>%
    .[rmax > 0 & !is.na(inflection_point)] %>%
    .[,rmax:=NULL] %>%
    .[,inflection_point:=NULL] %>%
    melt(., id.vars = NULL) %>%
    ggplot() +
    facet_wrap(~variable, scales = "free") +
    ylim(0, NA) +
    xlab("Value") +
    ylab("Count") +
    geom_histogram(aes(x = value, fill = variable), bins = 100, alpha = 0.5) +
    geom_hline(yintercept = 0) +
    viridis::scale_color_viridis("Metric", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    viridis::scale_fill_viridis("Metric", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"))

ggsave(filename = "rmax_prior.depletion_metrics.survival.png", plot = p, device = "png", path = plot_dir,
       scale = 1, width = 10, height = 8, units = "in", dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# 10. Weight data comparison plot
p = rmax_clean %>%
    .[, .(rmax,inflection_point,avg_weight_unfished, target_average_weight, sample_id)] %>%
    .[sample(1:.N, sample_size)] %>%
    .[rmax > 0 & !is.na(inflection_point)&avg_weight_unfished<250] %>%
    .[,rmax:=NULL] %>%
    .[,inflection_point:=NULL] %>%
    ggplot() +
    geom_point(aes(x = avg_weight_unfished, y = target_average_weight), alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    xlab("Average Weight (Unfished)") +
    ylab("Target Average Weight") +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"))

ggsave(filename = "rmax_prior.weight_comparison.survival.png", plot = p, device = "png", path = plot_dir,
       scale = 1, width = 8, height = 6, units = "in", dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Total parameter sets:", nrow(rmax_dt), "\n")
cat("Successful runs:", nrow(rmax_clean), "\n")
cat("Proportion successful:", round(nrow(rmax_clean)/nrow(rmax_dt), 3), "\n")
cat("Viable populations (Rmax > 0):", nrow(rmax_clean[rmax > 0]), "\n")
cat("Proportion viable:", round(nrow(rmax_clean[rmax > 0])/nrow(rmax_clean), 3), "\n")

cat("\n=== RMAX SUMMARY ===\n")
print(summary(rmax_clean$rmax))

cat("\n=== SPR SUMMARY ===\n")
print(summary(rmax_clean$spr))

cat("\n=== F ESTIMATE SUMMARY ===\n")
print(summary(rmax_clean$F_est))

cat("\nPlots saved to:", plot_dir, "\n")
