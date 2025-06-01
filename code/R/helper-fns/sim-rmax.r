

# Nicholas Ducharme-Barth
# 2025/05/29
# Calc Rmax from biological parameters
# Pardo et al. 2016 https://doi.org/10.1139/cjfas-2016-0069
# Hutchings et al. 2012 https://doi.org/10.1890/11-1313.1

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

sim_rmax = function(id, max_age, M_ref, L1, L2, vbk, age1, age2, cv_len, maturity_a, l50, weight_a, weight_b,sex_ratio,reproductive_cycle){
        require(data.table)
        set.seed(id)

        age_vector = 1:max_age

        # natural morality
            mortality_at_age = M_ref * (max_age / age_vector)^(-1)
        
        # survival
            survival_at_age = rep(NA,max_age)
            survival_at_age[1] = exp(-mortality_at_age[1])
            for(i in 2:length(age_vector)){
                survival_at_age[i] = survival_at_age[i-1]*exp(-mortality_at_age[i])
            }

        # length
            length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1))) 

        # maturity
            len_lower = seq(from=0,to=ceiling(max(length_at_age*(1+cv_len))),by=1)
            len_upper = len_lower + 1
            length_vec = len_lower + 0.5
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

            n_vec = 1000*survival_at_age
            alpha = sum(n_vec*((maturity_at_age*fecundity_at_age)/reproductive_cycle))/sum(n_vec) # average recruits per year at equilibrium
            
        # solve for Rmax
        rmax_function = function(rmax){
            euler_lotka_vec = rep(NA,max_age)
            for(i in 1:max_age){
                euler_lotka_vec[i] = survival_at_age[i]*reproduction_at_age[i]*exp(-rmax*i)
            }
            return((sum(euler_lotka_vec)-1)^2)
        }
        
        # t = nlminb(1/(0.4*min(which(maturity_at_age>0.5))), rmax_function, lower = 0, upper = 10)
        rmax = optim(0.2, rmax_function,method="Brent", lower = -1, upper = 1)$par

        if(length(which(maturity_at_age>0.5))>0){
            amat50 = min(which(maturity_at_age>0.5))
        } else {
            amat50 = NA
        }
        
        spr = sum(reproduction_at_age*survival_at_age)
        h = (alpha*spr)/(4+alpha*spr) # sex-ratio already included in calculation of spr
        generation_time = sum(age_vector*reproduction_at_age*survival_at_age)/spr # Grant and Grant 1992

        if(generation_time*rmax<0)
        {
            inflection_point = NA
        } else {
            inflection_point = 0.633 - 0.187*log(generation_time*rmax) # Fowler 1988
        }
        
        return(data.table(id=id,rmax=round(rmax,digits=4),epr_unfished=round(spr,digits=4),alpha=round(alpha,digits=4),h=round(h,digits=4),generation_time=round(generation_time,digits=4),inflection_point=round(inflection_point,digits=4),amat50=amat50,max_age=max_age))
}
