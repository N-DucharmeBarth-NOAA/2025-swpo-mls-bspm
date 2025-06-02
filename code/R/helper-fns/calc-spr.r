

# Nicholas Ducharme-Barth
# 2025/05/31
# Calculate the spr
# Assume knife-edge length-based selectivity

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

calc_spr = function(id, max_age, M_ref, L1, L2, vbk, age1, age2, cv_len, maturity_a, l50, weight_a, weight_b,sex_ratio,reproductive_cycle,F,selex_l50,selex_slope){
        require(data.table)
        set.seed(id)

        # age 
            age_vector = 1:max_age
        
        # length
            length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
            len_lower = seq(from=0,to=ceiling(max(length_at_age*(1+cv_len))),by=1)
            len_upper = len_lower + 1
            length_vec = len_lower + 0.5

        # natural morality
            mortality_at_age = M_ref * (max_age / age_vector)^(-1)
        
        # selectivity
            selex_length = 1 / (1 + exp(-selex_slope*(length_vec - selex_l50)))

        # fishing mortality
            f_at_length = F * selex_length
            
            pla_LA = pla_function(length(length_vec), length(age_vector), age_vector, len_lower, len_upper, L1, L2, vbk, age1, age2, cv_len)
            f_at_age = as.vector(matrix(f_at_length,nrow=1,ncol=length(length_vec)) %*% pla_LA)

        # survival
            survival_at_age_unfished = rep(NA,max_age)
            survival_at_age_unfished[1] = exp(-mortality_at_age[1])
            survival_at_age_fished = rep(NA,max_age)
            survival_at_age_fished[1] = exp(-(mortality_at_age[1]+f_at_age[1]))
            for(i in 2:length(age_vector)){
                survival_at_age_unfished[i] = survival_at_age_unfished[i-1]*exp(-mortality_at_age[i])
                survival_at_age_fished[i] = survival_at_age_fished[i-1]*exp(-(mortality_at_age[i]+f_at_age[i]))
            } 

        # maturity
            maturity_b = -maturity_a/l50
            maturity_at_length = (exp(maturity_a+maturity_b*length_vec)) / (1+exp(maturity_a+maturity_b*length_vec))

            # calc maturity at age using PLA
                pla_LA = pla_function(length(length_vec), length(age_vector), age_vector, len_lower, len_upper, L1, L2, vbk, age1, age2, cv_len)
                maturity_at_age = as.vector(matrix(maturity_at_length,nrow=1,ncol=length(length_vec)) %*% pla_LA)
        
        # weight
            weight_at_length = weight_a * length_vec ^ weight_b
            weight_at_age = as.vector(matrix(weight_at_length,nrow=1,ncol=length(length_vec)) %*% pla_LA)

        # fecundity (proxied by weight at age)
            fecundity_at_age = weight_at_age/max(weight_at_age)
        
        # reproductive output per year (sex-ratio & reproductive cycle)
            reproduction_at_age = (maturity_at_age * fecundity_at_age * sex_ratio)/reproductive_cycle

        # spr calc
        # lifetime average eggs per recruit in fished and unfished conditions
            epr_unfished = sum(reproduction_at_age*survival_at_age_unfished)
            epr_fished = sum(reproduction_at_age*survival_at_age_fished)
            spr = epr_fished/epr_unfished

            s_unfished = sum(survival_at_age_unfished*maturity_at_age)
            s_fished = sum(survival_at_age_fished*maturity_at_age)
            dep = s_fished/s_unfished

            sb_unfished = sum(survival_at_age_unfished*maturity_at_age*weight_at_age)
            sb_fished = sum(survival_at_age_fished*maturity_at_age*weight_at_age)
            dep_sb = sb_fished/sb_unfished


        return(data.table(id=id,spr=round(spr,digits=4),dep=round(dep,digits=4),dep_sb=round(dep_sb,digits=4),epr_fished=round(epr_fished,digits=4),epr_unfished=round(epr_unfished,digits=4)))
}
