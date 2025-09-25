# 2025 SWPO Striped Marlin Assessment

Bayesian surplus production model for Southwest Pacific Ocean (SWPO) striped marlin

Download SWPO 2025 assessment report:

- **2025 Stock Assessment of Striped Marlin in the Southwest Pacific Ocean: Part II – Bayesian Surplus Production Model**\
  **[WCPFC-SC21-SA-WP-07](https://meetings.wcpfc.int/node/26683)**

Download SWPO striped marlin 2025 diagnostic model and model ensemble:

- Clone the **[2025-swpo-mls-bspm](https://github.com/N-DucharmeBarth-NOAA/2025-swpo-mls-bspm)** repository

## Reference model

The *[diagnostic model](https://github.com/N-DucharmeBarth-NOAA/2025-swpo-mls-bspm/tree/main/data/output/model_runs/0100-dwfn-exeoFSTTGF-cf0.2-nb-qnewK-s1-o52b-o54b_0)* in this repository serves as the reference model that forms the basis of several sections and figures in the stock assessment report.

The diagnostic model is also the entry point when configuring and running the ensemble models that form the basis of scientific advice. When the ensemble includes specific factor levels (for model parameters, likelihood weights, etc.) the diagnostic model has intermediate levels, while other grid members explore higher and lower levels.

Finally, the diagnostic model is also the starting point for the next SWPO striped marlin stock assessment model development. One purpose of this repository is to give the stock assessor a good starting point that is organized and documented.

## Explore data, model settings, and results

The **data** folder includes all the Bayesian surplus production model **[input](data/input)** and **[output](data/output)**.

The **taf** folder extracts the data and results from model format to CSV format that can be examined using Excel, R, or other statistical software. [TAF](https://cran.r-project.org/package=TAF) is a standard reproducible format for stock assessments that is practical for making the model **[data](taf/data)** and **[output](taf/output)** tables available in a format that is easy to examine. The **[report](taf/report)** folder contains example plots.


## Run the assessment model

The SWPO striped marlin 2025 model takes around 30 minutes to run. 

Anyone can run the assessment analysis in TAF format. This project uses [`renv`](https://rstudio.github.io/renv/articles/renv.html) for package management to ensure reproducibility across environments. First use `renv` to install required packages.

#### Using Base R

1. Open an R terminal (recommended [version 4.4.2](https://cloud.r-project.org/) with RTools 4.4 installed)
2. Set your working directory to the cloned repository:
   ```r
   setwd("path/to/repository/")
   ```
3. Source the `.Rprofile` file to bootstrap `renv`:
   ```r
   source(".Rprofile")
   ```
4. If `renv` doesn't bootstrap automatically, run:
   ```r
   renv::restore()
   ```
5. Follow the prompts to install all required packages

#### Using RStudio

1. Open the repository as an [RStudio Project](https://bookdown.org/ndphillips/YaRrr/projects-in-rstudio.html)
2. RStudio should automatically detect the `renv` configuration and prompt for package installation
3. If `renv` doesn't bootstrap automatically, run:
   ```r
   renv::restore()
   ```

#### Using Visual Studio Code

1. Configure VS Code to work with R using the [R extension](https://github.com/REditorSupport/vscode-R)
2. Follow the specific [configuration steps for renv projects](https://github.com/REditorSupport/vscode-R/wiki/Working-with-renv-enabled-projects)
3. Open the repository folder in VS Code
4. Open an R terminal, which should prompt `renv` to bootstrap
5. If `renv` doesn't bootstrap automatically, run:
   ```r
   renv::restore()
   ```

### Run in TAF format

The full assessment model can be run as a TAF analysis. Start R and then run:

- Users seeking to reproduce a single example model should run the R script [`00a_start_here.R`](https://github.com/N-DucharmeBarth-NOAA/2025-swpo-mls-bspm/blob/main/taf/r_code/00a_start_here.R). Note that running the diagnostic case model will take ~30 minutes.
- Users seeking to reproduce the entire model ensemble should run the R script [`00b_start_here.R`](https://github.com/N-DucharmeBarth-NOAA/2025-swpo-mls-bspm/blob/main/taf/r_code/00b_start_here.R). Note that running the full model ensemble will take ~120 minutes.

This runs the assessment model and also extracts the data and output from the model files and makes the results available as CSV files, which can be examined and analyzed further.

### Note

Running the assessment model through the TAF framework will result in minor differences (affecting the 2nd or 3rd decimal place) of [key management quantities](https://github.com/N-DucharmeBarth-NOAA/2025-swpo-mls-bspm/blob/main/taf/output/ens_mgmt_summary.csv) despite all data, model settings and initial conditions being identical. This is due to differences in the random number generator used to set-up the parallel sampling of the posterior distribution. If the exact quantities are desired, as matching the report, please refer to the [main development version of the file](https://github.com/N-DucharmeBarth-NOAA/2025-swpo-mls-bspm/blob/main/data/output/ens_mgmt_summary.csv). 

## License

The code contained in this repository is licensed under the GNU GENERAL PUBLIC LICENSE version 3 ([GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)).

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
