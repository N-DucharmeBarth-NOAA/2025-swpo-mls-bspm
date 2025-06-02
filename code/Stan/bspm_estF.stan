// Original fletcher-schaefer surplus production model code comes from
// the Bayesian biomass dynamic model (https://github.com/cttedwards/bdm/) by Charles Edwards
// Original model code can be specifically found at: https://github.com/cttedwards/bdm/blob/master/R/bdm.R
// This version of the model uses a random effects formulation and estimates:
// logK, r, initial depletion, shape, process error, additional observation error, and fishing mortality error.
data {
    int T; // time dimension
    int I; // number of indices
    real index[T,I]; // matrix of indices
    real obs_removals[T]; // time series of catch
    real sigmao_mat[T,I]; // observation error corresponding to each index (scaled to mean 1)
    real sigmao_input; // observation error
    int lambdas[I]; // lambdas for CPUE
    real sigmac; // observation error for catch

    // priors
    real PriorMean_logK;
    real PriorSD_logK;
    real PriorMean_logr;
    real PriorSD_logr;
    real PriorMean_logsigmap;
    real PriorSD_logsigmap;
    real PriorMean_logshape;
    real PriorSD_logshape;
    real PriorSD_sigmao_add;
    real PriorSD_sigmaf;
}
transformed data{
    int Tm1;
    Tm1 = T-1;
}
parameters {
    // use mean-zero parameters with SD = 1
    // keeps searched parameter space on a similar scale
    real raw_logK;
    real raw_logr;
    real raw_logsigmap;
    real<lower=0> raw_sigmao_add;
    real raw_epsp[T];
    real raw_logshape;
    real<lower=0> raw_sigmaf;
    real<lower=0> raw_F[Tm1];
}
transformed parameters {
    // leading parameters
    real logK;
    real r;
    real dev[T]; // recruitment deviates
    real epsp[T]; // process error (multiplicative)
    real sigmap;
    real shape;
    real sigmaf;
    real F[Tm1];
    

    // process error
    sigmap = exp(raw_logsigmap*PriorSD_logsigmap + PriorMean_logsigmap); // lognormal prior
    real sigmap2;
    sigmap2 = pow(sigmap,2);

    // observation error
    real sigmao[T,I];
    real sigmao_sc;
    real sigmao_add;
    sigmao_add = raw_sigmao_add*PriorSD_sigmao_add;
    sigmao_sc = sigmao_input + sigmao_add;
    for(i in 1:T){
        dev[i] = raw_epsp[i]*sigmap;
        epsp[i] = exp(dev[i]-sigmap2/2);
        for(j in 1:I){
            sigmao[i,j] = sigmao_mat[i,j] * sigmao_sc;
        }
    }

    // non-centered parametrization 
    logK = raw_logK*PriorSD_logK + PriorMean_logK; // normal prior (log-scale)
    r = exp(raw_logr*PriorSD_logr + PriorMean_logr); // lognormal prior
    shape = exp(raw_logshape*PriorSD_logshape + PriorMean_logshape); // lognormal prior
    sigmaf = raw_sigmaf*PriorSD_sigmaf;
    for(t in 1:Tm1){
        F[t] = raw_F[t]*sigmaf; 
    }
    
    // fletcher-schaefer
    // parameters
    real n;
    real dmsy;
    real h;
    real m;
    real g;

    n = shape;
    dmsy = pow((1/n),(1/(n-1)));
    h = 2*dmsy;
    m = r*h/4;
    g = pow(n,(n/(n-1)))/(n-1);

    // population dynamics
        real removals[Tm1]; // catch (relative to logK) in each time step
        real x[T]; // population time series relative to logK
        x[1] = epsp[1];
        for(t in 2:T) {
            if(x[t-1]<=dmsy){
                x[t] = ((x[t-1] + r * x[t-1] * (1 - x[t-1]/h)))*epsp[t]*(exp(-F[t-1]));
                removals[t-1] = ((x[t-1] + r * x[t-1] * (1 - x[t-1]/h)))*epsp[t]*(1-exp(-F[t-1]))*exp(logK);
            } 
            if(x[t-1]> dmsy){
                x[t] = ((x[t-1] + g * m * x[t-1] * (1 - pow(x[t-1],(n-1)))))*epsp[t]*(exp(-F[t-1]));
                removals[t-1] = ((x[t-1] + g * m * x[t-1] * (1 - pow(x[t-1],(n-1)))))*epsp[t]*(1-exp(-F[t-1]))*exp(logK);
            } 
            // // Add safety checks and debugging
            // if(t >= T-5) {  // Debug last few time steps
            //     print("t = ", t, ", x[t-1] = ", x[t-1], ", x[t] = ", x[t]);
            //     print("  epsp[t] = ", epsp[t], ", F[t-1] = ", F[t-1]);
            //     print("  removals[t-1] = ", removals[t-1]);
            // }
        }

    // compute mpd catchability assuming 
    // variable sigmao over time assuming
    // uniform prior on q
        real q[I]; // catchability
        real sigmao2[T,I]; // observation error for each index
        real sum1;
        real sum2;
        real p;
        for(i in 1:I){
            sum1 = 0.0;
            sum2 = 0.0;
            p = 0.0;
            for(t in 1:T){
                sigmao2[t,i] = square(sigmao[t,i]);
                if(index[t,i]>0.0 && x[t]>0.0) {
                    sum1 = sum1 + log(index[t,i]/x[t])/sigmao2[t,i];
                    sum2 = sum2 + 1/sigmao2[t,i];
                    p = p + 1.0;
                }
            }
            if(p>2.0){
                q[i] = exp((0.5 * p + sum1) / sum2);
            } else{
                q[i] = 0.0;
            } 
        }
} 
model {
    // prior densities for
    // estimated parameters
        raw_logK ~ std_normal();
        raw_logr ~ std_normal();
        raw_epsp ~ std_normal();
        raw_logsigmap ~ std_normal();
        raw_sigmao_add ~ std_normal();
        raw_logshape ~ std_normal();
        raw_sigmaf ~ std_normal();
        raw_F ~ std_normal();
    
    // observation model
        real mu_index;
        for(i in 1:I){
            if(lambdas[i]==1){
                for(t in 1:T){
                    if(index[t,i]>0.0 && x[t]>0.0 && q[i]>0.0) {
                        mu_index = log(q[i]*x[t]) - sigmao2[t,i]/2;
                        // print("Location param value (index): ", mu_index,"; t = ",t);
                        target += lognormal_lpdf(index[t,i] | mu_index,sigmao[t,i]);
                    }
                }
            }
        }

        real mu_catch;
        for(t in 1:Tm1){
            mu_catch = log(removals[t]) - 0.5*sigmac^2;
            // print("log(removals[t]): ", log(removals[t]),"; t = ",t);
            // print("Location param value (catch): ", mu_catch,"; t = ",t);
            target += lognormal_lpdf(obs_removals[t]|mu_catch,sigmac);
        }
}
