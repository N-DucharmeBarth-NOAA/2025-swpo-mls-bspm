# Nicholas Ducharme-Barth
# 2025/06/11
# Sub-sample rmax distribution to match r distribution from posterior samples

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(data.table)
    library(ggplot2)
    library(viridis)
    library(magrittr)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
    proj_dir = this.path::this.proj()
    
#________________________________________________________________________________________________________________________________________________________________________________________________________
# make directories
    dir.create(file.path(proj_dir,"data","output"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(proj_dir,"plots"), showWarnings = FALSE, recursive = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load data
    rmax_exprior_id = readRDS(file.path(proj_dir,"data","output","rmax_exprior_id.rds"))
    prior_dt = fread(file.path(proj_dir,"data","output","bspm_parameter_priors_filtered.csv"))
    posterior_dt = fread(file.path(proj_dir,"data","output","model_runs","0003-2024cpueFPrior_0","hmc_samples.csv"))
    r_posterior = posterior_dt[variable == "r", value]

#________________________________________________________________________________________________________________________________________________________________________________________________________
# sub-sampling functions
    calc_distribution_similarity = function(x, y) {
        ks_test = ks.test(x, y)
        quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95)
        x_quant = quantile(x, quantiles, na.rm = TRUE)
        y_quant = quantile(y, quantiles, na.rm = TRUE)
        
        list(ks_statistic = ks_test$statistic,
             ks_pvalue = ks_test$p.value,
             quantile_diff = abs(x_quant - y_quant),
             mean_diff = abs(mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)),
             sd_diff = abs(sd(x, na.rm = TRUE) - sd(y, na.rm = TRUE)))
    }

    subsample_by_inverse_transform = function(source_values, target_values, n_samples = NULL) {
        if(is.null(n_samples)) n_samples = length(target_values)
        target_ecdf = ecdf(target_values)
        set.seed(123)
        u_vals = runif(n_samples)
        target_quantiles = quantile(target_values, u_vals, na.rm = TRUE)
        
        # Find indices of closest matches
        match_indices = sapply(target_quantiles, function(x) {
            which.min(abs(source_values - x))
        })
        
        return(match_indices)
    }

    subsample_by_importance_sampling = function(source_values, target_values, n_samples = NULL) {
        if(is.null(n_samples)) n_samples = length(target_values)
        source_density = density(source_values, na.rm = TRUE)
        target_density = density(target_values, na.rm = TRUE)
        
        source_dens_at_vals = approx(source_density$x, source_density$y, xout = source_values)$y
        target_dens_at_vals = approx(target_density$x, target_density$y, xout = source_values)$y
        
        source_dens_at_vals[is.na(source_dens_at_vals) | source_dens_at_vals <= 0] = min(source_dens_at_vals[source_dens_at_vals > 0], na.rm = TRUE)
        target_dens_at_vals[is.na(target_dens_at_vals) | target_dens_at_vals <= 0] = min(target_dens_at_vals[target_dens_at_vals > 0], na.rm = TRUE)
        
        weights = target_dens_at_vals / source_dens_at_vals
        weights = weights / sum(weights, na.rm = TRUE)
        
        set.seed(123)
        sample_indices = sample(seq_along(source_values), size = n_samples, replace = TRUE, prob = weights)
        return(sample_indices)
    }

    subsample_by_quantile_matching = function(source_values, target_values, n_strata = 10) {
        quantile_breaks = seq(0, 1, length.out = n_strata + 1)
        target_quantiles = quantile(target_values, quantile_breaks, na.rm = TRUE)
        target_counts = table(cut(target_values, breaks = target_quantiles, include.lowest = TRUE))
        
        sampled_indices = c()
        for(i in 1:n_strata) {
            stratum_mask = source_values >= target_quantiles[i] & source_values <= target_quantiles[i+1]
            indices_in_stratum = which(stratum_mask)
            
            if(length(indices_in_stratum) > 0 && target_counts[i] > 0) {
                n_to_sample = min(target_counts[i], length(indices_in_stratum))
                set.seed(123 + i)
                sampled_stratum = sample(indices_in_stratum, size = n_to_sample, replace = length(indices_in_stratum) < target_counts[i])
                sampled_indices = c(sampled_indices, sampled_stratum)
            }
        }
        return(sampled_indices)
    }

#________________________________________________________________________________________________________________________________________________________________________________________________________
# apply sub-sampling
    prior_original = prior_dt[!is.na(rmax) & is.finite(rmax)]
    prior_clean = prior_dt[!is.na(rmax) & is.finite(rmax) & id %in% rmax_exprior_id]
    rmax_clean = prior_clean$rmax
    id_clean = prior_clean$id
    r_clean = r_posterior[!is.na(r_posterior) & is.finite(r_posterior)]

    idx_subsample_it = subsample_by_inverse_transform(rmax_clean, r_clean)
    idx_subsample_is = subsample_by_importance_sampling(rmax_clean, r_clean)
    idx_subsample_qm = subsample_by_quantile_matching(rmax_clean, r_clean)
    
    rmax_subsample_it = rmax_clean[idx_subsample_it]
    rmax_subsample_is = rmax_clean[idx_subsample_is]
    rmax_subsample_qm = rmax_clean[idx_subsample_qm]

#________________________________________________________________________________________________________________________________________________________________________________________________________
# evaluate performance
    original_similarity = calc_distribution_similarity(rmax_clean, r_clean)
    it_similarity = calc_distribution_similarity(rmax_subsample_it, r_clean)
    is_similarity = calc_distribution_similarity(rmax_subsample_is, r_clean)
    qm_similarity = calc_distribution_similarity(rmax_subsample_qm, r_clean)

    ks_stats = c("Inverse Transform" = it_similarity$ks_statistic,
                 "Importance Sampling" = is_similarity$ks_statistic,
                 "Quantile Matching" = qm_similarity$ks_statistic)

    best_method = names(ks_stats)[which.min(ks_stats)]

    if(best_method == "Inverse Transform") {
        final_indices = idx_subsample_it
    } else if(best_method == "Importance Sampling") {
        final_indices = idx_subsample_is
    } else {
        final_indices = idx_subsample_qm
    }

#________________________________________________________________________________________________________________________________________________________________________________________________________
# create final dataset
    rmax_subsampled_dt = prior_clean[final_indices]
    rmax_subsampled_dt[, subsample_method := best_method]
    rmax_subsampled_dt[, clean_index := final_indices]

#________________________________________________________________________________________________________________________________________________________________________________________________________
# make plots
    plot_dt = rbindlist(list(
        data.table(value = rmax_clean, distribution = "Original rmax", method = "Prior"),
        data.table(value = r_clean, distribution = "Target r posterior", method = "Target"),
        data.table(value = rmax_subsample_it, distribution = "Inverse Transform", method = "Sub-sampled"),
        data.table(value = rmax_subsample_is, distribution = "Importance Sampling", method = "Sub-sampled"),
        data.table(value = rmax_subsample_qm, distribution = "Quantile Matching", method = "Sub-sampled")
    ))

    p1 = ggplot(plot_dt, aes(x = value, fill = distribution)) +
        geom_density(alpha = 0.6) +
        facet_wrap(~method, scales = "free_y") +
        viridis::scale_color_viridis("Distribution",begin = 0.1,end = 0.8,direction = -1,option = "H",discrete=TRUE,drop=FALSE) +
        viridis::scale_fill_viridis("Distribution",begin = 0.1,end = 0.8,direction = -1,option = "H",discrete=TRUE,drop=FALSE) +
        labs(title = "Distribution Comparison: rmax Sub-sampling Methods", x = "Value", y = "Density") +
        theme_minimal() +
        theme(text = element_text(size = 20),panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
              panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
              panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
              strip.background =element_rect(fill="white"),
              legend.key = element_rect(fill = "white"))

    ggsave(filename = "rmax_subsampling_comparison.png", plot = p1, path = file.path(proj_dir,"plots"),
           width = 12, height = 8, units = "in", dpi = 300)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# biological parameter comparison plots
    # Compare input parameters for original vs sub-sampled distributions
    p1_bio = prior_original[, .(id, max_age, M_ref, L1, L2, vbk, l50, sex_ratio, cv_len, weight_a, weight_b, h)] %>%
        .[, distribution := "Prior: all"] %>%
        melt(., id.vars = c("id", "distribution"))

    p2_bio = prior_clean[, .(id, max_age, M_ref, L1, L2, vbk, l50, sex_ratio, cv_len, weight_a, weight_b, h)] %>%
        .[, distribution := "Prior: filtered"] %>%
        melt(., id.vars = c("id", "distribution"))
    
    p3_bio = prior_clean[final_indices] %>%
        .[, .(id, max_age, M_ref, L1, L2, vbk, l50, sex_ratio, cv_len, weight_a, weight_b, h)] %>%
        .[, distribution := "Posterior: sub-sampled"] %>%
        melt(., id.vars = c("id", "distribution"))
    
    p_bio = rbind(p1_bio, p2_bio,p3_bio) %>%
        .[,distribution:=factor(distribution,levels=c("Prior: all","Prior: filtered","Posterior: sub-sampled"))] %>%
        ggplot() +
        facet_wrap(~variable, scales = "free") +
        ylim(0, NA) +
        xlab("Input") +
        ylab("Relative count") +
        geom_density(aes(x = value, fill = distribution), alpha = 0.5) +
        geom_hline(yintercept = 0) +
        viridis::scale_color_viridis("Distribution", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE) +
        viridis::scale_fill_viridis("Distribution", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE) +
        theme(text = element_text(size = 20),
              panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
              panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
              panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
              strip.background = element_rect(fill = "white"),
              legend.key = element_rect(fill = "white"))
    
    ggsave(filename = "rmax_subsampling_biological_parameters.png", plot = p_bio, device = "png", 
           path = file.path(proj_dir,"plots"), scale = 1, width = 16, height = 12, units = "in", dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# save results
    fwrite(rmax_subsampled_dt, file.path(proj_dir,"data","output","bspm_parameter_priors_subsampled_rmax.csv"))
    saveRDS(rmax_subsampled_dt, file.path(proj_dir,"data","output","bspm_parameter_priors_subsampled_rmax.rds"))
    
    summary_stats = list(
        original_rmax_summary = summary(rmax_clean),
        target_r_summary = summary(r_clean),
        subsampled_rmax_summary = summary(rmax_clean[final_indices]),
        method_used = best_method,
        similarity_metrics = list(original = original_similarity,
                                  final = calc_distribution_similarity(rmax_clean[final_indices], r_clean)),
        sample_sizes = list(original_rmax = length(rmax_clean),
                           target_r = length(r_clean),
                           subsampled_rmax = length(final_indices))
    )
    
    saveRDS(summary_stats, file.path(proj_dir,"data","output","rmax_subsampling_summary.rds"))
