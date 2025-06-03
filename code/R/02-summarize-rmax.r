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
dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")
plot_dir = file.path(proj_dir, "plots", "rmax-priors")
dir.create(plot_dir, recursive = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load data
rmax_dt = readRDS(file.path(proj_dir, "data", "output", "bspm_parameter_priors_parallel.rds"))

# Filter out failed runs for plotting
rmax_clean = rmax_dt[!is.na(rmax)]
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
    .[, .(sample_id, max_age, M_ref, L1, L2, vbk, l50, sex_ratio, cv_len,weight_a,weight_b,h)] %>%
    melt(., id.vars = "sample_id") %>%
    ggplot() +
    facet_wrap(~variable, scales = "free") +
    ylim(0, NA) +
    xlab("Input") +
    ylab("Count") +
    geom_density(aes(x = value, fill = variable), alpha = 0.5) +
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
p2 = rmax_clean[rmax > 0 & rmax<1.5, .(rmax)] %>%
    .[, distribution := "plausible"]

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
p1 = rmax_clean[, .(sample_id, max_age, M_ref, L1, L2, vbk, l50, sex_ratio, cv_len,weight_a,weight_b,h)] %>%
    .[, distribution := "all"] %>%
    melt(., id.vars = c("sample_id", "distribution")) 

p2 = rmax_clean[rmax > 0 & rmax<1.5] %>%
    .[, .(sample_id, max_age, M_ref, L1, L2, vbk, l50, sex_ratio, cv_len,weight_a,weight_b,h)] %>%
    .[, distribution := "plausible"] %>%
    melt(., id.vars = c("sample_id", "distribution"))

p = rbind(p1, p2) %>% 
    ggplot() +
    facet_wrap(~variable, scales = "free") +
    ylim(0, NA) +
    xlab("Input") +
    ylab("Relative count") +
    geom_density(aes(x = value, fill = distribution), alpha = 0.5) +
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

p = rmax_clean[, .(rmax, max_age, M_ref, L1, L2, vbk, l50, sex_ratio, cv_len,weight_a,weight_b,h)] %>%
    .[sample(1:.N, sample_size)] %>%
    .[, Outcome := "All"] %>%
    .[rmax > 0 & rmax<1.5, Outcome := "Plausible"] %>%
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
# 6. Pairs plot - derived quantities with outcomes
p = rmax_clean[, .(rmax, alpha, h, generation_time, inflection_point)] %>%
    .[sample(1:.N, sample_size)] %>%
    .[, Outcome := "All"] %>%
    .[rmax > 0 & rmax<1.5, Outcome := "Plausible"] %>%
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
# 7. Pairs plot - derived quantities for viable populations only
p = rmax_clean[, .(rmax, alpha, h, generation_time, inflection_point)] %>%
    .[sample(1:.N, sample_size)] %>%
    .[rmax > 0 & rmax<1.5 & !is.na(inflection_point)] %>%
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
# 11. NEW: Multi-facet biological relationships plot by outcome

# Function to calculate biological relationships for plotting
calc_biological_relationships = function(param_row) {
    # Extract parameters
    max_age = param_row$max_age
    M_ref = param_row$M_ref
    L1 = param_row$L1
    L2 = param_row$L2
    vbk = param_row$vbk
    age1 = param_row$age1
    age2 = param_row$age2
    cv_len = param_row$cv_len
    maturity_a = param_row$maturity_a
    l50 = param_row$l50
    weight_a = param_row$weight_a
    weight_b = param_row$weight_b
    
    # Age vector
    age_vector = 1:max_age
    
    # Length at age (Francis parameterization)
    length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
    
    # Length bins for calculations
    len_lower = seq(from=0, to=ceiling(max(length_at_age*(1+cv_len))), by=1)
    len_upper = len_lower + 1
    length_vec = len_lower + 0.5
    
    # Natural mortality at length (Lorenzen)
    mortality_at_length = M_ref * (length_vec / L2)^(-1)
    
    # Convert mortality at length to mortality at age using PLA
    pla_LA = pla_function(length(length_vec), length(age_vector), age_vector, len_lower, len_upper, L1, L2, vbk, age1, age2, cv_len)
    mortality_at_age = as.vector(matrix(mortality_at_length, nrow=1, ncol=length(length_vec)) %*% pla_LA)
    
    # Weight at length
    weight_at_length = weight_a * length_vec ^ weight_b
    
    # Maturity at length
    maturity_b = -maturity_a/l50
    maturity_at_length = (exp(maturity_a + maturity_b*length_vec)) / (1 + exp(maturity_a + maturity_b*length_vec))
    
    # Combine results
    results = list(
        # Length at age data
        length_age_data = data.table(
            sample_id = param_row$sample_id,
            age = age_vector,
            length = length_at_age,
            mortality = mortality_at_age,
            variable = "Length at Age",
            x_var = age_vector,
            y_var = length_at_age
        ),
        
        # Weight at length data
        weight_length_data = data.table(
            sample_id = param_row$sample_id,
            length = length_vec,
            weight = weight_at_length,
            variable = "Weight at Length",
            x_var = length_vec,
            y_var = weight_at_length
        )[length <= max(length_at_age) * 1.2], # Limit to reasonable length range
        
        # Maturity at length data
        maturity_length_data = data.table(
            sample_id = param_row$sample_id,
            length = length_vec,
            maturity = maturity_at_length,
            variable = "Maturity at Length",
            x_var = length_vec,
            y_var = maturity_at_length
        )[length <= max(length_at_age) * 1.2],
        
        # Natural mortality at age data
        mortality_age_data = data.table(
            sample_id = param_row$sample_id,
            age = age_vector,
            mortality = mortality_at_age,
            variable = "Natural Mortality at Age",
            x_var = age_vector,
            y_var = mortality_at_age
        )
    )
    
    return(results)
}

# Calculate relationships for a subset of parameter sets
cat("Calculating biological relationships for plotting...\n")
set.seed(456)
plot_sample_size = min(500, nrow(rmax_clean))  # Limit to avoid overplotting

# Stratified sampling to ensure representation across outcomes
rmax_sample = rbind(
    rmax_clean[sample(min(.N, plot_sample_size/3))],  # All
    rmax_clean[rmax > 0 & rmax<1.5][sample(min(.N, plot_sample_size/3))]  # Plausible
)

# Add outcome categories
rmax_sample[, outcome := "All"]
rmax_sample[rmax > 0 & rmax<1.5, outcome := "Plausible"]

# Calculate biological relationships for each parameter set
bio_relationships_list = list()
for(i in 1:nrow(rmax_sample)) {
    if(i %% 100 == 0) cat("Processing parameter set", i, "of", nrow(rmax_sample), "\n")
    
    param_row = rmax_sample[i]
    bio_data = calc_biological_relationships(param_row)
    
    # Combine all data types with outcome information
    combined_data = rbindlist(list(
        bio_data$length_age_data[, .(sample_id, variable, x_var, y_var, outcome = param_row$outcome)],
        bio_data$weight_length_data[, .(sample_id, variable, x_var, y_var, outcome = param_row$outcome)],
        bio_data$maturity_length_data[, .(sample_id, variable, x_var, y_var, outcome = param_row$outcome)],
        bio_data$mortality_age_data[, .(sample_id, variable, x_var, y_var, outcome = param_row$outcome)]
    ), fill = TRUE)
    
    bio_relationships_list[[i]] = combined_data
}

# Combine all biological relationship data
all_bio_data = rbindlist(bio_relationships_list, fill = TRUE)

# Create labels for facets with appropriate units
all_bio_data[, variable_label := factor(variable, 
    levels = c("Length at Age", "Weight at Length", "Maturity at Length", 
               "Natural Mortality at Age"),
    labels = c("Length at Age (cm)", "Weight at Length (kg)", "Maturity at Length", 
               "Natural Mortality at Age (yr⁻¹)"))]

all_bio_data[, outcome := factor(outcome, 
    levels = c("All", "Plausible"),
    labels = c("All", "Plausible"))]

# Create the multi-facet plot
p = all_bio_data %>%
    ggplot(aes(x = x_var, y = y_var, color = outcome)) +
    facet_wrap(~variable_label, scales = "free", ncol = 2) +
    geom_line(aes(group = interaction(sample_id, outcome)), alpha = 0.3, linewidth = 0.5) +
    labs(
        x = "Age (years) / Length (cm)",
        y = "Response Variable",
        color = "Outcome",
        title = "Biological Relationships by Population Outcome"
    ) +
    scale_color_manual("Outcome", values = c("All" = "blue", "Plausible" = "orange")) +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"),
          strip.text = element_text(size = 9),
          axis.text = element_text(size = 8),
          legend.position = "bottom")

ggsave(filename = "rmax_prior.biological_relationships_by_outcome.png", plot = p, device = "png", path = plot_dir,
       scale = 1, width = 12, height = 16, units = "in", dpi = 300, limitsize = TRUE)

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

cat("\nPlots saved to:", plot_dir, "\n")
