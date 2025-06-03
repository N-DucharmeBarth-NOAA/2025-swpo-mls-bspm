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
    model_run_dir = file.path(proj_dir, "data", "output")
    dir.create(model_run_dir, recursive = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)

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

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read in biological parameter results
    bio_params_dt = fread(file.path(proj_dir,"data","output","bspm_parameter_priors_filtered.csv"))


negative_log_likelihood = function(par, data, bio_params) {
    # Equation 8 implementation
    # define leading parameters
        Z1 = par[1]
        Z2 = par[2]

    # extract pars from bio_params
        max_age = bio_params$max_age
        L1 = bio_params$L1
        L2 = bio_params$L2
        vbk = bio_params$vbk
        age1 = bio_params$age1
        age2 = bio_params$age2
        weight_a = bio_params$weight_a
        weight_b = bio_params$weight_b
        Lc = bio_params$Lc
    
    # define W_inf and Wc
        age_vector = seq(from=0,to=bio_params$max_age*2,length.out=50)
        length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
        W_inf = weight_a*max(length_at_age)^weight_b
        Wc = weight_a*length_at_age[max(which(length_at_age<Lc))]^weight_b

    # process data
        proc_dat = data[weight>Wc&year>=1952] %>%
                   .[,.(mean_wt=mean(weight),sd_wt=sd(weight),.N),by=year] %>%
                   .[,d:=year-1952]

    # calculate d_vec
        d_vec = proc_dat$d

    # Calculate predicted lengths (using your transitional equation)
        predicted = transitional_mean_weight(W_inf,Wc,vbk,Z1,Z2,d_vec)
    
        log_likelihood = sum(dnorm(proc_dat$mean_wt, 
                             mean = predicted, 
                             sd = proc_dat$sd_wt/sqrt(proc_dat$N), 
                             log = TRUE))
    
    
    return(-log_likelihood)  # Return negative for minimization
}

# Optimization
    tmp_bio_params = bio_params_dt[1]
    result = nlminb(start = c(tmp_bio_params$M_ref+0.1, tmp_bio_params$M_ref+0.2),
                objective = negative_log_likelihood,
                data = nz_wt,
                bio_params = tmp_bio_params,
                lower = rep(tmp_bio_params$M_ref,2),
                upper = rep(3,3))

# calculate relative depletion from Z1 to Z2
    rel_depletion = function(result,bio_params){
        # define parameters
            Z1 = result$par[1]
            Z2 = result$par[2]
            max_age = bio_params$max_age
            L1 = bio_params$L1
            L2 = bio_params$L2
            vbk = bio_params$vbk
            age1 = bio_params$age1
            age2 = bio_params$age2
            weight_a = bio_params$weight_a
            weight_b = bio_params$weight_b
            cv_len = bio_params$cv_len
            maturity_a = bio_params$maturity_a
            l50 = bio_params$l50
            sex_ratio = bio_params$sex_ratio
            reproductive_cycle = bio_params$reproductive_cycle
            
        # age 
            age_vector = 1:max_age
        
        # length
            length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
            len_lower = seq(from=0,to=ceiling(max(length_at_age*(1+cv_len))),by=1)
            len_upper = len_lower + 1
            length_vec = len_lower + 0.5

        # survival
            survival_at_age_Z1 = rep(NA,max_age)
            survival_at_age_Z1[1] = 1
            survival_at_age_Z2 = rep(NA,max_age)
            survival_at_age_Z2[1] = 1
            for(i in 2:length(age_vector)) {
                survival_at_age_Z1[i] = survival_at_age_Z1[i-1] * exp(-Z1)
                survival_at_age_Z2[i] = survival_at_age_Z2[i-1] * exp(-Z2)
            }
            survival_at_age_Z1[max_age] = survival_at_age_Z1[max_age-1] / (1 - exp(-Z1))
            survival_at_age_Z2[max_age] = survival_at_age_Z2[max_age-1] / (1 - exp(-Z2))

        # maturity
            maturity_b = -maturity_a/l50
            maturity_at_length = (exp(maturity_a+maturity_b*length_vec)) / (1+exp(maturity_a+maturity_b*length_vec))

            # calc maturity at age using PLA
                pla_LA = pla_function(length(length_vec), length(age_vector), age_vector, len_lower, len_upper, L1, L2, vbk, age1, age2, cv_len)
                maturity_at_age = as.vector(matrix(maturity_at_length,nrow=1,ncol=length(length_vec)) %*% pla_LA)
                maturity_at_age = maturity_at_age/max(maturity_at_age)
        
        # weight
            weight_at_length = weight_a * length_vec ^ weight_b
            weight_at_age = as.vector(matrix(weight_at_length,nrow=1,ncol=length(length_vec)) %*% pla_LA)

        
        # reproductive output per year (sex-ratio & reproductive cycle)
            reproduction_at_age = (maturity_at_age * weight_at_age * sex_ratio)/reproductive_cycle

        # spr calc
        # lifetime average eggs per recruit in fished and unfished conditions
            epr_Z1 = sum(reproduction_at_age*survival_at_age_Z1)
            epr_Z2 = sum(reproduction_at_age*survival_at_age_Z2)
            rel_dep_n = sum(survival_at_age_Z2)/sum(survival_at_age_Z1)
            rel_dep_ssb = epr_Z2/epr_Z1
            return(data.table(sample_id=bio_params$sample_id,rel_dep_n = rel_dep_n, rel_dep_ssb = rel_dep_ssb))
    }

    rel_depletion(result,tmp_bio_params)

# calculate relative depletion from Z1 to Z2
    mean_wt_resid = function(result,bio_params,data){
        
        # define parameters
            Z1 = result$par[1]
            Z2 = result$par[2]
            max_age = bio_params$max_age
            L1 = bio_params$L1
            L2 = bio_params$L2
            vbk = bio_params$vbk
            age1 = bio_params$age1
            age2 = bio_params$age2
            weight_a = bio_params$weight_a
            weight_b = bio_params$weight_b
            Lc = bio_params$Lc
    
        # define W_inf and Wc
            age_vector = seq(from=0,to=bio_params$max_age*2,length.out=50)
            length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
            W_inf = weight_a*max(length_at_age)^weight_b
            Wc = weight_a*length_at_age[max(which(length_at_age<Lc))]^weight_b

        # process data
            proc_dat = data[weight>Wc&year>=1952] %>%
                    .[,.(mean_wt=mean(weight),sd_wt=sd(weight),.N),by=year] %>%
                    .[,d:=year-1952]

        # calculate d_vec
            d_vec = proc_dat$d

        # Calculate predicted lengths (using your transitional equation)
            predicted = transitional_mean_weight(W_inf,Wc,vbk,Z1,Z2,d_vec)
        
        proc_dat$pred_wt = predicted
        proc_dat$resid = proc_dat$mean_wt - proc_dat$pred_wt
        proc_dat$sample_id=bio_params$sample_id

        return(proc_dat[,.(sample_id,year,d,N,mean_wt,sd_wt,pred_wt,resid)])
    }

    out = mean_wt_resid(result,bio_params,data)       
