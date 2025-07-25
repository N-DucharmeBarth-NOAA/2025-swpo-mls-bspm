# Nicholas Ducharme-Barth
# 2025/07/25
# Calculate highest-density intervals (HDI) for management quantities

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(data.table)
    library(magrittr)
    library(HDInterval)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
    proj_dir = this.path::this.proj()

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define helper functions
    calc_hdi = function(name,data,mass=0.95){
        require(data.table)
        require(magrittr)
        require(HDInterval)

        hdi = hdi(data, credMass = mass)
        eti = quantile(data,probs=c(0.5*(1-mass),1-(0.5*(1-mass))))
        names(hdi) = names(eti) = c("lower","upper")

        interval_dt = data.table(metric=rep(name,2)) %>%
                  .[,mass:=rep(mass,2)] %>%
                  .[,type:=c("eti","hdi")] %>%
                  .[,median:=rep(median(data),2)] %>%
                  .[,lower:=c(eti[1],hdi[1])] %>%
                  .[,upper:=c(eti[2],hdi[2])]
        return(interval_dt)
    }
    
    plot_intervals = function(name, data, mass = 0.95, type = c("eti", "hdi", "both"),dens_n=1e3,...) {
        type = match.arg(type)
        
        # Calculate intervals using your function
        interval_dt = calc_hdi(name, data, mass)
        
        # Create density data
        dens = density(data,n=dens_n)
        x = dens$x
        y = dens$y
        
        # Get intervals
        eti = interval_dt[type == "eti"]
        hdi = interval_dt[type == "hdi"]
        
        # Set up plot layout
        if (type == "both") {
            par(mfrow = c(2,1), mar = c(4,4,3,1))
            plot_both = TRUE
        } else {
            par(mfrow = c(1,1), mar = c(5,4,4,2))
            plot_both = FALSE
        }
        
        # Plot ETI panel
        if (type %in% c("eti", "both")) {
            plot(x, y, type = "l", lwd = 2, col = "black",
                main = ifelse(plot_both, "A. Equal-tailed interval", 
                            paste("Equal-tailed interval (", mass*100, "%)", sep="")),
                xlab = name, ylab = "Probability density",...)
            
            # Fill areas
            polygon(c(x[x <= eti$lower], eti$lower), c(y[x <= eti$lower], 0), 
                    col = rgb(1,0,0,0.7), border = NA)
            polygon(c(x[x >= eti$upper], rev(x[x >= eti$upper])), 
                    c(y[x >= eti$upper], rep(0, sum(x >= eti$upper))), 
                    col = rgb(0.5,0,0,0.7), border = NA)
            polygon(c(x[x >= eti$lower & x <= eti$upper], rev(x[x >= eti$lower & x <= eti$upper])), 
                    c(y[x >= eti$lower & x <= eti$upper], rep(0, sum(x >= eti$lower & x <= eti$upper))), 
                    col = rgb(0,0,1,0.7), border = NA)
            
            lines(x, y, lwd = 2, col = "black")
            
            # Add median line
            med_y = approx(x, y, eti$median)$y
            segments(eti$median, 0, eti$median, med_y, col = "darkblue", lwd = 2)
            
            legend("topright", legend = c(paste(mass*100, "% credible", sep=""), "Lower tail", "Upper tail", "Median"), 
                fill = c(rgb(0,0,1,0.7), rgb(1,0,0,0.7), rgb(0.5,0,0,0.7), "darkblue"), bty = "n")
        }
        
        # Plot HDI panel  
        if (type %in% c("hdi", "both")) {
            plot(x, y, type = "l", lwd = 2, col = "black",
                main = ifelse(plot_both, "B. Highest density interval", 
                            paste("Highest density interval (", mass*100, "%)", sep="")),
                xlab = name, ylab = "Probability density",...)
            
            # Fill areas
            polygon(c(x[x <= hdi$lower], hdi$lower), c(y[x <= hdi$lower], 0), 
                    col = rgb(1,0,0,0.7), border = NA)
            polygon(c(x[x >= hdi$upper], rev(x[x >= hdi$upper])), 
                    c(y[x >= hdi$upper], rep(0, sum(x >= hdi$upper))), 
                    col = rgb(0.5,0,0,0.7), border = NA)
            polygon(c(x[x >= hdi$lower & x <= hdi$upper], rev(x[x >= hdi$lower & x <= hdi$upper])), 
                    c(y[x >= hdi$lower & x <= hdi$upper], rep(0, sum(x >= hdi$lower & x <= hdi$upper))), 
                    col = rgb(0,0,1,0.7), border = NA)
            
            # HDI threshold line
            hdi_indices = which(x >= hdi$lower & x <= hdi$upper)
            if(length(hdi_indices) > 0) {
                hdi_y = min(y[hdi_indices])
            } else {
                hdi_y = min(approx(x, y, c(hdi$lower, hdi$upper))$y, na.rm = TRUE)
            }
            abline(h = hdi_y, col = "gray", lwd = 1.5)
            
            lines(x, y, lwd = 2, col = "black")
            
            # Add median line
            med_y = approx(x, y, hdi$median)$y
            segments(hdi$median, 0, hdi$median, med_y, col = "darkblue", lwd = 2)
            
            legend("topright", legend = c(paste(mass*100, "% credible", sep=""), "Lower tail", "Upper tail", "HDI threshold", "Median"), 
                fill = c(rgb(0,0,1,0.7), rgb(1,0,0,0.7), rgb(0.5,0,0,0.7), "gray", "darkblue"), bty = "n")
        }
        
        # Reset plot parameters
        par(mfrow = c(1,1))
        
        # Return the interval data
        invisible(interval_dt)
    }

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read in data
    ens_derived_dt = fread(file=file.path(proj_dir,"data","output","ens_derived_dt.csv"))
    unique(ens_derived_dt$name)

    target_metrics = "recent_D_Dmsy"
    target_dist = na.omit(ens_derived_dt[name==target_metric]$value)

    target_metrics = c("msy","Dmsy","Pmsy","Fmsy","Plrp",
                       "latest_depletion","latest_population","latest_F","latest_removals",
                       "latest_D_Dmsy","latest_F_Fmsy",
                       "recent_D","recent_P","recent_F","recent_removals",
                       "recent_D_Dmsy","recent_F_Fmsy")
    hdi_dt.list = as.list(rep(NA,length(target_metrics)))
    for(i in 1:length(target_metrics)){
        hdi_dt.list[[i]] = calc_hdi(target_metrics[i],na.omit(ens_derived_dt[name==target_metrics[i]]$value),mass=0.95)
    }
    hdi_dt = rbindlist(hdi_dt.list)
    fwrite(hdi_dt,file=file.path(proj_dir,"data","output","ens_hdi_dt.csv"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# plot
    target_metric = "recent_D_Dmsy"
    target_dist = na.omit(ens_derived_dt[name==target_metric]$value)

    png(filename = file.path(proj_dir,"plots","hdi_v_eti.d_dmsy.png"),width = 6, height = 6, units = "in",bg = "white",  res = 300)
        plot_intervals(expression(D[recent]/D[MSY]), target_dist, mass = 0.95, type = c("both"))
    dev.off()
