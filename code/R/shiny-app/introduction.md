# BSPM Model Analysis - Interactive Explorer

This Shiny application provides comprehensive interactive visualization and comparison of Bayesian Surplus Production Model (BSPM) runs.

## Features

### Model Summary
Compare key statistics and parameter estimates across different model runs in an interactive table format.

### HMC Diagnostics
Examine MCMC convergence and diagnostic plots including:

**Main Diagnostics:**
- **Trace plots**: Monitor chain mixing and convergence
- **Parallel coordinates**: Visualize parameter relationships across iterations
- **Pairs plots**: Examine parameter correlations and potential issues
- **Autocorrelation plots**: Check parameter correlation structure

**Convergence Diagnostics:**
- **R-hat diagnostics**: Assess convergence across chains
- **Effective sample size**: Check sampling efficiency

### Posterior Predictive Checks
Comprehensive model adequacy assessment with separate sections for CPUE and Catch data:

**CPUE Posterior Predictive Checks:**
- **Density overlays**: Compare observed vs predicted distributions
- **Empirical CDF**: Cumulative distribution comparisons
- **PIT ECDF**: Probability integral transform diagnostics
- **Test statistics**: Summary statistic comparisons
- **LOO diagnostics**: Leave-one-out cross-validation plots (PIT, Q-Q, intervals)

**Catch Posterior Predictive Checks:**
- Complete parallel suite for catch data validation
- Same diagnostic approaches applied to catch observations

### Model Fits
Detailed examination of model performance against observed data:

**Index Fits:**
- **Basic fits**: Predicted vs observed abundance indices
- **Fits with PPD**: Include posterior predictive distributions
- **Residual analyses**: Identify patterns in model residuals

**Catch Fits:**
- **Basic fits**: Model performance against catch data
- **Fits with PPD**: Full uncertainty representation
- **Residual analyses**: Diagnostic residual patterns

### Management Plots
Standard fisheries management visualizations:
- **Kobe plots**: Stock status relative to reference points
- **Majuro plots**: Alternative stock status visualization  
- **Forecasts**: Project future stock trajectories under different scenarios
- **Prior-posterior comparisons**: Examine how data updates prior beliefs for both parameters and time series

## Important Notes

- **Single vs Multiple Models**: Some diagnostic plots (HMC diagnostics, LOO analyses) are only available when a single model is selected
