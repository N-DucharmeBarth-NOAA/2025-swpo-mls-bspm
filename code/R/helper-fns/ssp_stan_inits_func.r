            stan_inits_func = function(Tm1){
                        inits = list(raw_logK=rnorm(1,0,0.25),
                        raw_logr=rnorm(1,0,0.25),
                        raw_epsp=rnorm(Tm1+1,0,0.25),
                        raw_logsigmap=rnorm(1,0,0.25),
                        raw_sigmao_add=abs(rnorm(1,0,0.25)),
                        raw_logshape=rnorm(1,0,0.25),
                        raw_sigmaf=abs(rnorm(1,0,0.25)),
                        raw_F=abs(rnorm(Tm1,0,0.25)))
                    return(inits)
            }
