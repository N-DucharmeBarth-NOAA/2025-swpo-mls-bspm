

# NOAA PIFSC
# 2025/09/25
# R code to re-produce the model ensemble from the 2025 WCPFC SWPO striped marlin stock assessment
# using the Transparent Assessment Framework (TAF) <https://www.ices.dk/data/assessment-tools/Pages/transparent-assessment-framework.aspx>
# Model 0100: DWFN index; baseline priors
# Model 0102: NZ index; baseline priors
# Model 0105: DWFN index; alternative shape prior
# Model 0107: NZ index; alternative shape prior

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

#________________________________________________________________________________________________________________________________________________________________________________________________________
# run the TAF analysis
    source.taf(file.path("code","01_data.R")) # format the data for the Bayesian State-Space Surplus Production model (BSPM)
    source.taf(file.path("code","02b_model.R")) # run the diagnostic case model (~4 minutes runtime; single-threaded)
    # source.taf("code/03a_output.R") # summarize output (medians) from the representative model
    # source.taf("code/04a_report.R") # make some plots for the representative model
