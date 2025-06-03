# Nicholas Ducharme-Barth
# 2025/06/03
# calculate fishing mortality rates based on change in mean size using
# Gedamke and Hoenig 2006 DOI: 10.1577/T05-153.1
# https://fluke.vims.edu/hoenig/pdfs/Gedamke_and_Hoenig_length_based_Z.pdf
# and assuming population demographics
# Use fishing mortality to calculate SPR (depletion proxy)
# Assumes knife-edge selectivity

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
    library(data.table)
    library(magrittr)

#_____________________________________________________________________________________________________________________________
# define helper functions
    transitional_mean_weight = function(W_inf, Wc, vbk, Z1, Z2, d) {
    # Equation 6
    # modified to be in terms of weight
    # d time in years
    
    # Calculate each component separately for clarity
    term1 = W_inf
    
    numerator_fraction = Z1 * Z2 * (W_inf - Wc) * 
                        (Z1 + vbk + (Z2 - Z1) * exp(-(Z2 + vbk) * d))
    
    denominator_fraction = (Z1 + vbk) * (Z2 + vbk) * 
                            (Z1 + (Z2 - Z1) * exp(-Z2 * d))
    
    mean_weight = term1 - (numerator_fraction / denominator_fraction)
    
    return(mean_weight)
    }

#_____________________________________________________________________________________________________________________________
# define paths
    proj_dir = this.path::this.proj()
    dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")

#_____________________________________________________________________________________________________________________________
# bring in nz rec data
    nz_wt = fread(file.path(proj_dir,"data","input","nz-all-sport-club-weights.csv"))
    colnames(nz_wt) = c("date","club","species","tagged","weight","boat","locality","year")
    nz_wt = nz_wt %>%
                .[!is.na(weight)&!is.na(year)&tagged=="",.(date,year,weight)] %>%
                .[year<1988] %>%
                .[,dd:=sapply(date,function(x)as.numeric(strsplit(x,"/")[[1]][2]))] %>%
                .[,mm:=sapply(date,function(x)as.numeric(strsplit(x,"/")[[1]][1]))] %>%
                .[,yy:=sapply(date,function(x)as.numeric(strsplit(x,"/")[[1]][3]))] %>%
                .[yy>51&yy<100,yy:=yy+1900] %>%
                .[yy<25,yy:=yy+2000] %>%
                .[,year:=yy] %>%
                .[,month:=c(2,2,2,5,5,5,8,8,8,11,11,11)[mm]] %>%
                .[,.(year,month,weight)] %>%
                .[,ts:=paste0(year,"-",month)] %>%
                na.omit(.)
