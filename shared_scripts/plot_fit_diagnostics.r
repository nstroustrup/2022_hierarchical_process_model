library("quantreg")
library("MASS")
library("rms")
#Run a set of regressions to compare the length of vigorous movement span to lifespan and remaining lifespan.
#If "plot_intermediate_states is set, any movement span is also compared to lifespan and remaining lifespan.
#if make_any_plots is set, nice regression diagnostic plots are drawn.
#if calc_boot is T, bootstrapping p values are calculated to see if the slope of the linear model
#relating vigorous movement span to lifespan is less than one.
#set exclude_extreme_quantiles to 1 2 3 to exclude extreme quantiles of death times, vigorous movement times, or both
#from analysis.  In practice, never do this; we include to satisfy reviewer.

plot_fit_diagnostics = function(spans,experiments_to_plot,only_plot_control,plot_intermediate_states,calc_boot,make_any_plots=T,exclude_extreme_quantiles="none",quantiles_to_exclude=NULL){
  par(cex.axis=1, cex.lab=1, cex.main=1, cex.sub=1)
	cr_res = NULL;
	glm_results = NULL
	qr_results = NULL

	for (ii in 1:length(experiments_to_plot)){
	  
  fam = gaussian("identity")
  resid_fun = function(x)x;
		t1 = Sys.time();
		ex = experiments_to_plot[ii]
		sp_ex = subset(spans,spans$Experiment == ex  & spans$expansion != 0 & spans$slow_age != 0);
		sp_ex$excluded_from_fit = F;
		
		groups = unique(sp_ex$group)
		colors = rainbow(length(groups));
		for (g in groups){
			print(paste("Considering ",ex,g));
			flush.console();
			sp = sp_ex[sp_ex$group==g,]
			devices = unique(sp$Device);
			device_colors = rainbow(length(devices));
			
			make_plots = (!only_plot_control || sp$control[1] == T) && make_any_plots
			cex.lab = .8
			point_col=rgb(0,0,0,.5)
			cex = .4;
			sp$Device = factor(sp$Device)
			res = ns_span_device_adjustment(sp);
			sp = res[["spans"]];
			if (exclude_extreme_quantiles != "none"){
			 
			  if (exclude_extreme_quantiles %in% c("death","both")){
			    qt_e = quantile(sp$expansion_v_adj,c(0.05,0.95))
			    if (quantiles_to_exclude == "both"){
  			    sp$excluded_from_fit[sp$expansion_v_adj < qt_e[1] | sp$expansion_v_adj > qt_e[2] ]=T
			    }else if (quantiles_to_exclude == "top"){
			      sp$excluded_from_fit[sp$expansion_v_adj > qt_e[2] ]=T
			    }else{
			      sp$excluded_from_fit[sp$expansion_v_adj < qt_e[1] ]=T
			    }
			  }
			  if (exclude_extreme_quantiles %in% c("vmc","both")){
			    qt_s = quantile(sp$slow_age_v_adj,c(0.05,0.95))
			    if (quantiles_to_exclude == "both"){
    			  sp$excluded_from_fit[sp$slow_age_v_adj < qt_s[1]| sp$slow_age_v_adj > qt_s[2]] = T
			    } else if (quantiles_to_exclude == "top"){
			      sp$excluded_from_fit[sp$slow_age_v_adj > qt_s[2] ]=T
			    }else{
			      sp$excluded_from_fit[sp$slow_age_v_adj < qt_s[1] ]=T
			    }
			 }
			}
			#First, we fit models and make plots comparing vigorous movement and movement to death times
			if (make_plots){
				#PLOT RAW VIG vs LIFE
				m = max(sp$slow_age_v_adj,sp$expansion_v_adj)
				print("Plotting vig adj vs death-time expansion")
				plot(sp$slow_age_v_adj,sp$expansion_v_adj,pch=21,bg=point_col,cex=cex,col=NA,xlim=c(0,m),ylim=c(0,m),xlab="Vigorous Movement Cessation Time (days)",ylab="Death Time (days)",cex.lab=cex.lab)
				
				if (exclude_extreme_quantiles %in% c("death","both")){
				  if (quantiles_to_exclude %in% c("bottom","both"))
				    abline(h=qt_e[1],lty=2,col="gray")
				    if (quantiles_to_exclude %in% c("top","both"))
				  abline(h=qt_e[2],lty=2,col="gray")
				  }
				if (exclude_extreme_quantiles %in% c("vmc","both")){
				  if (quantiles_to_exclude %in% c("bottom","both"))
				  abline(v=qt_s[1],lty=2,col="gray")
				  if (quantiles_to_exclude %in% c("top","both"))
				  abline(v=qt_s[2],lty=2,col="gray")
				}
				title(paste(ex,g))
				#title(ex)
			}
			options(datadist=NULL)
			data_dist <<- datadist(sp);
			options(datadist="data_dist")
			if (length(unique(sp$Device))==1){
				reg_formula = expansion~slow_age
				dev_adj_formula = NULL;
			}else{
				reg_formula = expansion~slow_age+slow_age*Device
				dev_adj_formula = expansion~Device;
			}
      sp = sp[sp$excluded_from_fit ==F,]
			results = quantify_relationship(reg_formula,dev_adj_formula,sp,fam,"vigorous_vs_death",calc_boot);
			gr = results$results
			reg = results$reg;
			gr$replicate_id = which(ex==experiments_to_plot);
			gr$Experiment = ex;
			gr$group=g;
			gr$food = sp$food[1];	
			gr$temperature = sp$temperature[1];	
			glm_results =  rbind(gr,glm_results);
			print(ex)
			print(summary(reg))
			
			if (make_plots)
				abline(a=0,b=1,lwd=2,lty=2,col="gray")

			
			
			if (0){
				#stratify the analysis across each each device (ie scanner)
				if (make_plots){
					rng = quantile(sp$slow_age[sp$group==g],c(.01,.99))
					newx = seq(rng[1],rng[2],.25)
					if (sum(sp$group == g) < 80)
						next
					for ( d in devices){
						pr = predict(reg,newdata=data.frame(slow_age = newx,group=rep(g,length(newx)),Device = rep(d,length(newx))))
						lines(newx,resid_fun(pr),col=device_colors[which(d == devices)],lwd=2)
					}
				}
			}else{
			  #run the analysis looking for the consensus fit across all devices
			  sm = summary(reg);
				if (make_plots)
				abline(a=reg[["coefficients"]][["Intercept"]], b = reg[["coefficients"]][["slow_age"]],col="red",lwd=2)

			}
			if (plot_intermediate_states){
			
				if (make_plots){
					#PLOT RAW MOVEMENT vs LIFE
					m = max(sp$movement_cessation_m_adj,sp$expansion_m_adj)
					plot(sp$movement_cessation_m_adj,sp$expansion_m_adj,pch=21,bg=point_col,cex=cex,col=NA,xlim=c(0,m),ylim=c(0,m),xlab="Weak Movement Cessation Time (days)",ylab="Death Time (days)",cex.lab=cex.lab)
				 	title(paste(ex,g))
					abline(a=0,b=1,lwd=2,lty=2,col="gray")
				}

				options(datadist=NULL)
				data_dist <<- datadist(sp);
				options(datadist="data_dist")
				if (length(unique(sp$Device))==1){
					reg_formula = expansion~movement_cessation
					dev_adj_formula = NULL;
				}else{
					reg_formula = expansion~movement_cessation+movement_cessation*Device
					dev_adj_formula = expansion~Device;
				}
				
				results = quantify_relationship(reg_formula,dev_adj_formula,sp,fam,"movement_vs_death",calc_boot);
				
				gr_m = results$results
				reg_m = results$reg;
				if (!is.null(reg_m)){
  				gr_m$replicate_id = which(ex==experiments_to_plot);
  				gr_m$Experiment = ex;
  				gr_m$group=g;
  				gr_m$food = sp$food[1];	
  				gr_m$temperature = sp$temperature[1];		
  				glm_results =  rbind(gr_m,glm_results);
  				
  				if (0){
  					#run the analysis looking for the consensus fit across all devices
  					if (make_plots){
  						rng = quantile(sp$slow_age[sp$group==g],c(.01,.99))
  						newx = seq(rng[1],rng[2],.25)
  						if (sum(sp$group == g) < 80)
  							next
  						for ( d in devices){
  							pr = predict(reg_m,newdata=data.frame(slow_age = newx,group=rep(g,length(newx)),Device = rep(d,length(newx))))
  							lines(newx,resid_fun(pr),col=device_colors[which(d == devices)],lwd=2)
  						}
  					}
  				}else{
  				  #run the analysis looking for the consensus fit across all devices
  					if (make_plots)
  					abline(a=reg_m[["coefficients"]][["Intercept"]], b = reg_m[["coefficients"]][["movement_cessation"]],col="red",lwd=2)
  				}
				}
			}
			
			if (0){
				#plot quantile regression
				taus = c(.25,.5,.75)
				tau_color = rainbow(length(taus));
				for (tau in taus){
					qres = rq(expansion~slow_age_v_adj,sp,tau=tau);
					cf= summary(qres,se="boot")$coefficients
					cf2 = summary(qres,se="rank")$coefficients;
					if (make_plots){
						col = color[which(tau == taus)]
						if (cf["slow_age_v_adj","Pr(>|t|)"] < .05){
							abline(a=cf2["(Intercept)",1],b=cf2["slow_age_v_adj",1],col=col)
						}
					}
				}
			}
			
			#Second, we fit models and make plots comparing vigorous movement and movement to the remaining lifespan
			#sp = sp_ex[sp_ex$group==g,
			#sp$Device = factor(sp$Device)
			#res = ns_span_device_adjustment(sp);]
			spp = sp#[sp$group==g,]
			if (length(unique(spp$Device))==1){
				ee = spp$expansion;
				ss = spp$slow_age;
			}else{
				
				vr = ns_single_bj_group(Surv(spp$slow_age,1-spp$Censored)~spp$Device,);
				cf = vr$coefficients
				ee = spp$expansion;
				ss = spp$expansion
				for (i in 1:(dim(cf)[1])){
					cur = spp$Device==as.character(cf$bj_group[i])
					ee[cur] = ee[cur]*exp(-cf$coefficient[i])
					ss[cur] = ss[cur]*exp(-cf$coefficient[i])
				}	
			}

			cr = corr(data.frame(ee=ee,ss=ss))
			dif = median(sp$expansion[sp$group==g]-sp$slow_age[sp$group==g])
			mean_vig = median(sp$expansion[sp$group==g]);
			cr_res = rbind(data.frame(lbl=sp$lbl[1],experiment=ex,group=g,cr=cr,dif=dif,mean_vig=mean_vig,N=(dim(spp)[1])),cr_res)

			res = resid_fun(residuals(reg))
			pp = resid_fun(predict(reg))
			
			if (make_plots){
				#PLOT RESIDUAL of VIG vs RESIDUAL LIFE
				plot(pp,res,pch=21,bg=point_col,cex=cex,col=NA,log="x",xlab="Predicted Death time (days)",ylab="Residual Death Time (days)")
			  	title(paste(ex,g))
				if (quantile(sp[sp$group == g,"slow_age"],.25) == 0)
					next
				rr = res[sp$group == g]
				ppp = pp[sp$group == g]
				if (length(rr) < 80)
					next
				los = loess(rr~ppp,span=1);
				rng = quantile(ppp,c(.05,.95));
				newx = seq(rng[1],rng[2],.25)
				pr = predict(los,newdata=data.frame(ppp = newx))
				lines(newx,pr,col=colors[which(g == groups)],lwd=2)
				legend("topright",title=paste(ex,g),bty="n",col=colors[g==groups],legend=g,lty=1,cex=.9)
			}


		
				
				sp$d_col = apply(sp,1,FUN=function(x,devices,device_colors){return(device_colors[which(x[["Device"]]==devices)])},devices,device_colors)
				#only plot control groups (usually wild-type animals at 25C)
				
				options(datadist=NULL)
				data_dist <<- datadist(sp);
				options(datadist="data_dist")
				results = quantify_relationship(ev_d~slow_age+Device*slow_age,ev_d~Device,sp,fam,"d_vs_slow_age",F);
				gr_rl = results$results
				gr_rl$replicate_id = which(ex==experiments_to_plot);
				gr_rl$Experiment = ex;
				gr_rl$group=g;
				gr_rl$food = sp$food[1];	
				gr_rl$temperature = sp$temperature[1];		
				glm_results =  rbind(gr_rl,glm_results);
				#Try some fancy quantile regression (not included in manuscript)
				if (make_plots){
					#PLOT DEVICE-ADJUSTED VIG vs (LIFE-VIG) with QUANTILE REGRESSION
					plot(sp$slow_age_v_adj,sp$ev_d,pch=21,bg=point_col,cex=cex,col=NA, xlab="Vigorous Movement Cessation Time (days)",ylab="Remaining Lifespan (days)",cex.lab=cex.lab)
				 	title(paste(ex,g))
				}
				
			
				los = loess(ev_d~slow_age_v_adj,data=sp,span=1);
				rng = quantile(sp$slow_age_v_adj,c(.001,.999));
				newx = seq(rng[1],rng[2],.25)
				pr = predict(los,newdata=data.frame(slow_age_v_adj = newx))
				if (make_plots)
				lines(newx,pr,lwd=1,col="green")
				reg_coef = coef(results$reg)
				abline(a=reg_coef[["Intercept"]],b=reg_coef[["slow_age"]],col="red")
				#browser()
				
				taus = seq(.1,.9,.2)
				tau_color = rainbow(length(taus));
				for (tau in taus){
					qres = rq(ev_d~slow_age_v_adj,sp,tau=tau);
					cf= summary(qres,se="boot")$coefficients
					cf2 = summary(qres,se="rank")$coefficients;
					if (0 && make_plots){
						col = tau_color[which(tau == taus)]
						if (cf["slow_age_v_adj","Pr(>|t|)"] < .05){
							abline(a=cf2["(Intercept)",1],b=cf2["slow_age_v_adj",1],col=col)
						
						}
					}
					if (tau == .5){
						qr_results = rbind(data.frame(stage="vigorous_vs_death",
									      Experiment = ex, group=g,
									      food = sp$food[1],
									      tau = tau,
									      Intercept = cf2["(Intercept)",1],
									      Intercept_li=cf2["(Intercept)",2],
									      Intercept_ui=cf2["(Intercept)",3],
									      val = cf2["slow_age_v_adj",1],
									      val_li=cf2["slow_age_v_adj",2],
									      val_ui=cf2["slow_age_v_adj",3],
									      p=cf["slow_age_v_adj","Pr(>|t|)"],
									      N = dim(sp)[1])
									      ,qr_results);
					}
				}
				if (0 && make_plots){
					legend("topright",legend=taus,col=color,lty=1,bty="n",cex=.8)
				}
				
			
				
				#Try some fancy quantile regression (not included in manuscript)

				if (plot_intermediate_states){
				  options(datadist=NULL)
				  data_dist <<- datadist(sp);
				  options(datadist="data_dist")	
				  results = quantify_relationship(em_d~movement_cessation+Device*movement_cessation,em_d~Device,sp,fam,"dm_vs_movement",F);
				  if ( make_plots){
					#PLOT DEVICE-ADJUSTED MOVEMENT vs (LIFE-VIG) with QUANTILE REGRESSION
					plot(sp$movement_cessation_m_adj,sp$em_d,pch=21,bg=point_col,cex=cex,col=NA, xlab="Weak Movement Cessation Time (days)",ylab="Remaining Lifespan (days)",cex.lab=cex.lab)
				  title(paste(ex,g))
					reg = rlm(em_d~movement_cessation_m_adj,data=sp);
					sm = summary(reg);
					abline(a=sm$coefficients["(Intercept)",1], b = sm$coefficients["movement_cessation_m_adj",1],col="red",lwd=2)
					los = loess(em_d~movement_cessation_m_adj,data=sp,span=1);
					rng = quantile(sp$movement_cessation_m_adj,c(.001,.999));
					newx = seq(rng[1],rng[2],.25)
					pr = predict(los,newdata=data.frame(movement_cessation_m_adj = newx))
					lines(newx,pr,lwd=1,col="green")
					
					if (0){
						taus = seq(.1,.9,.2)
						tau_color = rainbow(length(taus));
						for (tau in taus){
							qres = rq(em_d~movement_cessation_m_adj,sp,tau=tau);
							cf= summary(qres,se="boot")$coefficients
							cf2 = summary(qres,se="rank")$coefficients;
							if (0 && make_plots){
								col = tau_color[which(tau == taus)]
								if (cf["movement_cessation_m_adj","Pr(>|t|)"] < .05){
									abline(a=cf2["(Intercept)",1],b=cf2["movement_cessation_m_adj",1],col=col)
								}
							}
							if (tau == .5){
								qr_results = rbind(data.frame(stage="movement_vs_death",
											      Experiment = ex, group=g,
											      food = sp$food[1],
											      tau = tau,
											      Intercept = cf2["(Intercept)",1],
											      Intercept_li=cf2["(Intercept)",2],
											      Intercept_ui=cf2["(Intercept)",3],
											      val = cf2["movement_cessation_m_adj",1],
											      val_li=cf2["movement_cessation_m_adj",2],
											      val_ui=cf2["movement_cessation_m_adj",3],
											      p=cf["movement_cessation_m_adj","Pr(>|t|)"],
											      N = dim(sp)[1])
											      ,qr_results);
							}
						}
					}
				  }
				}
				if (make_plots){
					if (1){
						#PLOT DEVICE-ADJUSTED VIG vs (LIFE-VIG)  with nice colors
						plot(sp$slow_age_v_adj,sp$ev_d,pch=21,bg=sp$d_col,cex=cex,col=NA,xlab="Vigorous Movement Cessation Time (days)",ylab="Remaining Lifespan (days)",cex.lab=cex.lab)
					 	title(paste(ex,g))
					}	

					#PLOT CDF TO COMPARE EFFECT OF DEVICES ON VIG
					plot(c(),c(),xlim=range(sp$slow_age_v_adj),ylim=c(0,1),type="n",xlab="Vigorous Movement Cessation Time (days)",ylab="CDF")
					title(paste(ex,g))
					for (d in devices){
						cd = ecdf(sp$slow_age_v_adj[sp$Device==d]);
						lines(cd,col=device_colors[which(d==devices)]);
					}
					#PLOT CDF TO COMPARE EFFECT OF DEVICES on DEATH
					plot(c(),c(),xlim=range(sp$expansion_v_adj),ylim=c(0,1),type="n",xlab="Death Time (days)",ylab="CDF")
					title(paste(ex,g))
					for (d in devices){
						cd = ecdf(sp$expansion_v_adj[sp$Device==d]);
						lines(cd,col=device_colors[which(d==devices)]);
					}
				}
				if (0){
					qt = quantile(sp$slow_age_v_adj,c(1,2)/3)
					sp$split = "upper";
					sp$split[sp$slow_age_v_adj < qt[2]] = "middle";
					sp$split[sp$slow_age_v_adj < qt[1]] = "lower";
					splts = c("upper","middle","lower")
					sp$split = factor(sp$split,levels=splts)
					spl_colors=c("black","turquoise","red")
					sv = survfit(Surv(sp$expansion,1-sp$Censored)~split,data=sp)
					plot(sv,col=spl_colors,xlab="Age (days)",ylab="Fraction Surviving")
					legend("topright",bty="n",legend=splts,col=spl_colors,lty=1,title=paste(ex,g))

					sv = survfit(Surv(sp$ev_d,1-sp$Censored)~split,data=sp)
					plot(sv,col=spl_colors,xlab="Time(days)",ylab="Fraction Surviving")
					legend("topright",bty="n",legend=splts,col=spl_colors,lty=1,title=paste(ex,g))
					
					
					sp$d_raw = sp$expansion-sp$slow_age
					sv = survfit(Surv(sp$d_raw,1-sp$Censored)~split,data=sp)
					plot(sv,col=spl_colors,xlab="Time (days)",ylab="Fraction Surviving")
					legend("topright",bty="n",legend=splts,col=spl_colors,lty=1,title=paste(ex,g))
				}
			}
			
		}
		t2 = Sys.time();
		print(difftime(t2,t1,units="min"))
	return(list("cr_res" =cr_res,
      	"glm_results" = glm_results,
      	"qr_results" = qr_results));
}