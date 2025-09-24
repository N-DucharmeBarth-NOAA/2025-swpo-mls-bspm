

# NOAA PIFSC
# 2025/09/23
# R code to re-produce the diagnostic model from the 2025 WCPFC SWPO striped marlin stock assessment
# using the Transparent Assessment Framework (TAF) <https://www.ices.dk/data/assessment-tools/Pages/transparent-assessment-framework.aspx>
# Model 0100: DWFN index; baseline priors

# Copyright (c) 2025 NOAA PIFSC
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library("TAF")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# change working directory to ./taf/
    setwd("taf")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# initialize the TAF analysis
# Process the SOFTWARE.bib & DATA.bib  metadata files to set up the files required for the analysis.
# SOFTWARE.bib contains the Stan model code which will be compiled into executables used to fit the model and R helper functions
# DATA.bib contains the input data needed to fit the models
    taf.boot()

# rename the executable
    file.rename(from = file.path("boot","software","bspm_estqsimple_softdep_fullmvprior_x0_sttgamma_flexsigmaC_OPT.stan"),
     to = file.path("boot","software","bspm.stan"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# run the TAF analysis
    source.taf(file.path("code","01_data.R")) # format the data for the Bayesian State-Space Surplus Production model (BSPM)
    # source.taf("code/02a_model.R") # run a representative model from the model ensemble (~4 minutes runtime; single-threaded)
    # source.taf("code/03a_output.R") # summarize output (medians) from the representative model
    # source.taf("code/04a_report.R") # make some plots for the representative model
