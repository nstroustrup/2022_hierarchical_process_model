library(grDevices)

ns_plot_pca_set_standard_plot_size = function(number_of_pcas_to_display){
	ns_plot_pca_figure_width = 4;
	ns_plot_pca_figure_aspect_ratio = 1.41;
	n = number_of_pcas_to_display
	windows(width = n*ns_plot_pca_figure_width,height=4*ns_plot_pca_figure_width/ns_plot_pca_figure_aspect_ratio,rescale="fixed");
	layout(matrix(1:(3*n),nrow=3,ncol=n))
	par(mar = c(4, 4, 2, 2))
}

get_individual_spans = function(data,observation_start_time){
	d = data[ !is.na(data$By.Hand.Movement.Cessation.Relative.Time) &
			!is.na(data$Best.Guess.Death.Associated.Expansion.Time) &
			!is.na(data$Best.Guess.Death.Associated.Expansion.Time),]
			
	m = data.frame(group=d$group,Plate.Name=d$Plate.Name,Group.ID=d$Group.ID,
			movement_cessation = d$First.Observation.Relative.Time-d$By.Hand.Movement.Cessation.Relative.Time-observation_start_time,
			expansion = d$First.Observation.Relative.Time-d$Best.Guess.Death.Associated.Expansion.Time-observation_start_time,
			slow_age = d$First.Observation.Relative.Time-d$Best.Guess.Fast.Movement.Cessation.Relative.Time-observation_start_time)
	
	if (any(m$movement_cessation < 0) || any(m$expansion < 0) || any(m$slow_age < 0))
		warning("observation start time is cropping data!");
	#crop any events that happened before observation started.
	m$movement_cessation[m$movement_cessation < 0] = 0
	m$expansion[m$expansion < 0] = 0
	m$slow_age[m$slow_age < 0] = 0
	
	if ("Experiment.Name" %in% names(d))
		m$Experiment=d$Experiment.Name;
	m$wrm = paste(m$Plate.Name,m$Group.ID)
	rows_to_grab = c();
	for (i in unique(m$wrm)){
		rows_to_grab = c(which(m$wrm==i)[1],rows_to_grab);
	}
	m2 = m[rows_to_grab,]
	mmm = m2$movement_cessation<m2$slow_age
	m2$slow_age[mmm] = m2$movement_cessation[mmm];	#bug with by-hand annotations of slow age time
	
	
	m2$duration_alive_but_nonmoving = m2$expansion-m2$movement_cessation;
	m2$duration_alive_but_nonmoving_crop = m2$duration_alive_but_nonmoving;
	m2$duration_alive_but_nonmoving_crop[m2$duration_alive_but_nonmoving_crop<0]=0;
	m2$duration_changing_posture = m2$movement_cessation-m2$slow_age
	
	m2$status = 1;
	return(m2);
}

#generates the pretty survival curves with state info overlaid
plot_pretty_survival_spans = plot_spans  = function(spans,legend_title){
	x = 1:length(spans$slow_age);
	spans$total=spans$slow_age+spans$duration_changing_posture+spans$duration_alive_but_nonmoving_crop;
	d = spans[order(spans$total,spans$slow_age,spans$duration_changing_posture,spans$duration_alive_but_nonmoving_crop),]
	max_y = max(d$slow_age + d$duration_changing_posture+d$duration_alive_but_nonmoving);
	plot(c(),c(),type="n",ylim=c(0,max(x)+1),xlim=c(0,max_y),xaxt="n",xlab="Time (days)",ylab="Survival")
	axis(1,at=seq(0,max_y,2),xlab="Time (days)")
	for (i in x){
		a = d[i,];
		xx =c(i,i+1,i+1,i,i);
		polygon(c(0,0,a$slow_age,a$slow_age,0),max(x)-xx,col="purple",border=NA);
		polygon(a$slow_age+c(0,0,a$duration_changing_posture,a$duration_changing_posture,0),max(x)-xx,col="orange",border=NA);
		polygon(a$slow_age+a$duration_changing_posture+c(0,0,a$duration_alive_but_nonmoving_crop,a$duration_alive_but_nonmoving_crop,0),max(x)-xx,col="dark red",border=NA);
	}
	lines(c(0,d$total),c(max(x),max(x)-x),col="black",lty=1)
	legend("topright",title=legend_title,legend=" ",bty="n")
	d = data_hmm[order(data_hmm$slow_age/data_hmm$total,data_hmm$duration_changing_posture/data_hmm$total),]
	plot(c(),c(),type="n",xlim=c(0,max(x)),ylim=c(0,1),xlab="Animal ID",ylab="Fraction of lifespan")
	for (i in x){
		a = d[i,];
		tot = a$slow_age+a$duration_changing_posture+a$duration_alive_but_nonmoving_crop;
		xx =c(i,i+1,i+1,i,i);
		polygon(xx,c(0,0,a$slow_age,a$slow_age,0)/tot,col="purple",border=NA);
		polygon(xx,(a$slow_age+c(0,0,a$duration_changing_posture,a$duration_changing_posture,0))/tot,col="orange",border=NA);
		polygon(xx,(a$slow_age+a$duration_changing_posture+c(0,0,a$duration_alive_but_nonmoving_crop,a$duration_alive_but_nonmoving_crop,0))/tot,,col="dark red",border=NA);
	}
	slow = mean(d$slow_age)
	posture = mean(d$duration_changing_posture);
	duration_alive_but_nonmoving = mean(d$duration_alive_but_nonmoving);
	max_y =duration_alive_but_nonmoving+posture+slow;
	plot(c(),c(),type="n",ylim=c(0,1),xlim=c(0,max_y),ylab="Survival",axt="n",xlab="Time (days)")
	axis(1,at=seq(0,max_y,1),xlab="Time (days)")
	xx= c(0,1,1,0,0);
	polygon(c(0,0,slow,slow,0),xx,col="purple",border=NA);
	polygon(c(0,0,posture,posture,0)+slow,xx,col="orange",border=NA);
	polygon(c(0,0,duration_alive_but_nonmoving,duration_alive_but_nonmoving,0)+slow+posture,xx,col="dark red",border=NA);
}

#compares span lengths as they covary across individuals
pretty_plot_vmc_lifespan_spans_comparison =compare_individual_spans = function(spans,observation_start_time,proportional_to_group_size=F){
	if ( is.logical(observation_start_time))
		stop("compare_individual_spans()::You must specifiy an observation start time");
	max_y =1;
	groups = unique(spans$group)
	groups = groups[order(as.character(groups))]
	spans$slow_age= spans$slow_age-observation_start_time;
	spans$slow_age= spans$slow_age-observation_start_time;
	spans$slow_age[spans$slow_age<0]=0;
	if (proportional_to_group_size){
		xlim=c(0,dim(spans)[1]+length(groups))
		x_group_step = aggregate(slow_age~group,data=spans,FUN=length);
		x_group_step = x_group_step[order(as.character(x_group_step$group)),]
		x_group_pos = cumsum(x_group_step$slow_age+1)-x_group_step$slow_age/2
	}else{
		group_width = 5;
		group_border = .5;
		xlim=c(0,(group_width+group_border)*length(groups))
		x_group_pos = (0:(length(groups)-1)+.5)*(group_width+group_border);
	}
	plot(c(),c(),type="n",xlim=xlim,ylim=c(0,max_y),ylab="Fraction of Life",xaxt="n",xlab="Time (days)",yaxt="n")
	spans$duration_changing_posture_crop = spans$duration_changing_posture;
	spans$duration_changing_posture_crop[spans$duration_changing_posture_crop < 0] = 0;
	
	axis(1,at= x_group_pos,labels=groups)
	axis(2,at=seq(0,1,.1))
	group_cutoff = 0;
	for (gg in groups){
		ss = spans[spans$group == gg,]
		ss$tot = ss$slow_age+ss$duration_changing_posture_crop+ss$duration_alive_but_nonmoving_crop;

		ss = ss[order(ss$slow/ss$tot),]
		x = 1:length(ss$slow_age);
		if (length(ss$slow_age) < 5)
			next;
		for (i in x){
				a = ss[i,];
				if (proportional_to_group_size){
					xx =c(i,i+1,i+1,i,i)+group_cutoff;
				}else{
					xx = c(i,i+1,i+1,i,i)/dim(ss)[1]*group_width+ (which(gg == groups)-1)*(group_width+group_border);
				}
				polygon(xx,c(0,0,a$slow_age,a$slow_age,0)/a$tot,col="purple",border=NA);
				polygon(xx,(a$slow_age+c(0,0,a$duration_changing_posture_crop,a$duration_changing_posture_crop,0))/a$tot,col="orange",border=NA);
				polygon(xx,(a$slow_age+a$duration_changing_posture_crop+c(0,0,a$duration_alive_but_nonmoving_crop,a$duration_alive_but_nonmoving_crop,0))/a$tot,col="dark red",border=NA);
		}
		group_cutoff = group_cutoff+length(ss$group)+1;
		#print(length(ss$group))
	}
}

#generates a plot comparing the average span length across individuals
compare_average_spans = function(spans){
	slow = aggregate(slow_age~group,data=spans,FUN=mean)
	posture = aggregate(duration_changing_posture~group,data=spans,FUN=mean) 
	duration_alive_but_nonmoving = aggregate(duration_alive_but_nonmoving_crop~group,data=spans,FUN=mean)
	stats = data.frame(group=slow$group,slow=slow$slow_age,posture=posture$duration_changing_posture,duration_alive_but_nonmoving=duration_alive_but_nonmoving$duration_alive_but_nonmoving_crop)
	stats$total = stats$duration_alive_but_nonmoving+stats$posture+stats$slow;
	max_y =1;#max(stats$duration_alive_but_nonmoving+stats$posture+stats$slow);
	g = unique(spans$group);
	plot(c(),c(),type="n",ylim=c(0,length(g)),xlim=c(0,max_y),ylab="Survival",xaxt="n",xlab="Time (days)",yaxt="n")
	axis(1,at=seq(0,max_y,1))
	axis(2,at=c(.5,1.5),label=stats$group)
	for (gg in stats$group){
		ss = stats[stats$group == gg,]
		pos = which(gg==stats$group)-1
		yy= c(pos,pos+.95,pos+.95,pos,pos);

		polygon(c(0,0,ss$slow,ss$slow,0)/ss$total,yy,col="purple",border=NA);
		polygon((c(0,0,ss$posture,ss$posture,0)+ss$slow)/ss$total,yy,col="orange",border=NA);
		polygon((c(0,0,ss$duration_alive_but_nonmoving,ss$duration_alive_but_nonmoving,0)+ss$slow+ss$posture)/ss$total,yy,col="dark red",border=NA);
	}
}



#plots a pca.res generated by the fda library, using custom coloring and shaded intervals
ns_plot_pca_fd = function (x, nx = 128, pointplot = F, harm = 0, expand = 0, shading= 0, cycle = FALSE, colors=rainbow(length(expand)),ltys=rep(1,length(expand)),lines_not_shape=F,...) 
{
    ns_plot_pca_figure_aspect_ratio = 1.41;
    shading_cols = col2rgb(colors)/255;
    pcafd <- x
    if (!(inherits(pcafd, "pca.fd"))) 
        stop("Argument 'x' is not a pca.fd object.")
    harmfd <- pcafd[[1]]
    basisfd <- harmfd$basis
    rangex <- basisfd$rangeval
    if (length(nx) > 1) {
        argvals <- nx
        nx <- length(x)
    }
    else {
        argvals <- seq(rangex[1], rangex[2], length = nx)
    }
    fdmat <- eval.fd(argvals, harmfd)
    meanmat <- eval.fd(argvals, pcafd$meanfd)
    dimfd <- dim(fdmat)
    nharm <- dimfd[2]
    plotsPerPg <- sum(par("mfrow"))
    harm <- as.vector(harm)
    if (harm[1] == 0) 
        harm <- (1:nharm)
        
    for (jharm in harm) {
	    #if (jharm == 2) {
		#op <- par(ask = TRUE)
		#on.exit(par(op))
	   # }
	    iharm <- jharm
	    vecharm <- fdmat[, iharm]
	    pcmat  <- matrix(rep(meanmat,length(expand)),ncol=length(expand))+ 
		      matrix(rep(vecharm,length(expand)),ncol=length(expand))*matrix(rep(expand,length(meanmat)),byrow=T,ncol=length(expand))

	    pcmat_shading  <- matrix(rep(meanmat,length(shading)),ncol=length(shading))+ 
		      matrix(rep(vecharm,length(shading)),ncol=length(shading))*matrix(rep(shading,length(meanmat)),byrow=T,ncol=length(shading))

	    if (pointplot) 
		plottype <- "p"
	    else plottype <- "l"
	    percentvar <- round(100 * pcafd$varprop[iharm], 1)
	    rng2 = range(c(range(pcmat),range(pcmat_shading)))

	    plot(c(), c(), type = "n", ylim = rng2, xlim=range(argvals),
	    	xlab="Time (days)",ylab = "Body size (relative)", #asp=ns_plot_pca_figure_aspect_ratio,
		main = paste("PCA function", iharm, "(Percentage of variability", 
		  percentvar, ")"))
	    abline(v=0,col="gray")
	    abline(h=1,col="gray")
	    if (pointplot) {
		for (i in seq(1,dim(pcmat_shading)[2],2)){
			cc=rgb(shading_cols[1,ceiling(i/2)],shading_cols[2,ceiling(i/2)],shading_cols[3,ceiling(i/2)],alpha=.1);
			if (lines_not_shape){
				cc="#000000";
				lines(argvals, pcmat_shading[, i], col=cc)
				lines(argvals, pcmat_shading[,i+1], col=cc,lty=2)
			}
			else{
				polygon(c(argvals,rev(argvals)), c(pcmat_shading[, i],rev(pcmat_shading[,i+1])), col=cc,border=NA,)
				lines(argvals, pcmat_shading[,i+1], col=cc,,lty=2,lwd=2)
			}
		}
		for (i in 1:dim(pcmat)[2]){
			points(argvals, pcmat[, i], lty = ltys[i],col=colors[i],lwd=2)
		}
	    }
	    else {
		for (i in seq(1,dim(pcmat_shading)[2],2)){
			cc=rgb(shading_cols[1,ceiling(i/2)],shading_cols[2,ceiling(i/2)],shading_cols[3,ceiling(i/2)],alpha=.1);
			if (lines_not_shape){
				cc="#000000";
				lines(argvals, pcmat_shading[, i], col=cc)
				lines(argvals, pcmat_shading[,i+1], col=cc,lty=2)
			}
			else{
				lines(argvals, pcmat_shading[,i+1], col=cc,lty=2,lwd=2)
				polygon(c(argvals,rev(argvals)), c(pcmat_shading[, i],rev(pcmat_shading[,i+1])), col= cc,border=NA)
			}
		}
		for (i in 1:dim(pcmat)[2]){
			lines(argvals, pcmat[, i], lty = ltys[i],col=colors[i],lwd=2)
		}
	    }
	}
  
   # invisible(NULL)
}

#plots a few different animals from mutant and wildtype group, selected as the 50, 75, and 90th quantile largest changing animals.
ns_plot_representative_animals = function(animal_scores,animals_to_use,data2,harm,control_group,test_group,xlim=NULL){

	ns_plot_pca_figure_aspect_ratio = 1.41;
	
	pca_column = paste0("pca_",harm);
	wt_animals_to_use = subset(animal_scores,animal_scores$group==control_group & animal_scores$animal %in% animals_to_use);
	mut_animals_to_use = subset(animal_scores,animal_scores$group==test_group& animal_scores$animal %in% animals_to_use);
	
	#wt_animals_to_use$distance_from_median = abs(wt_animals_to_use$pca_1 - medians[1]);
	#mut_animals_to_use$distance_from_median = abs(mut_animals_to_use$pca_1 - medians[2]);

	wt_animals_to_use = wt_animals_to_use[order(wt_animals_to_use[,pca_column]),]
	mut_animals_to_use = mut_animals_to_use[order(mut_animals_to_use[,pca_column]),]

	animals_to_plot = c(as.character(wt_animals_to_use$animal[c( round(dim(wt_animals_to_use)[1]/2),
									round(3*dim(wt_animals_to_use)[1]/4),
									round(9*dim(wt_animals_to_use)[1]/10)
								    )
								  ]),
			   as.character(mut_animals_to_use$animal[c( round(dim(mut_animals_to_use)[1]/2),
									round(3*dim(mut_animals_to_use)[1]/4),
									round(9*dim(mut_animals_to_use)[1]/10)
								    )
								  ]))
	animal_rank = c(3,2,1,3,2,1);
	data_to_plot= subset(data2,data2$animal %in% animals_to_plot);
	if (is.null(xlim))
		xlim = range(data_to_plot$time)
	rng=NULL;
	for (i in animals_to_plot){
		dd = subset(data_to_plot,data_to_plot$animal==i)
		sm = ksmooth(dd$time,dd$nI)
		#print(c(i,range(sm$y[!is.na(sm$y)])))	
		rng = range(c(sm$y[!is.na(sm$y)],rng))	#data can have gaps in measurement, which show up as NA values in sm$y when the gap width is larger than the smoothing kernel
	}	
	plot(c(),c(),xlim=xlim,ylim=rng,type="n",xlab="Time (days)",ylab="Relative animal Size")#,asp=ns_plot_pca_figure_aspect_ratio)
	    abline(v=0,col="gray")
	    abline(h=1,col="gray")
	for (i in animals_to_plot){
		dd = subset(data_to_plot,data_to_plot$animal==i)
		#dd = subset(dd,dd$By.Hand.Post.mortem.Contraction.Stop.Time < 0 & dd$Best.Guess.Fast.Movement.Cessation.Relative.Time> 0)
		cur_rank = animal_rank[which(i == animals_to_plot)];
		if (cur_rank == 1){ lty=1;lwd=1;}
		if (cur_rank == 2){ lty=2;lwd=1;}
		if (cur_rank == 3){ lty=3;lwd=1;}
		col="red";
		if (dd$group[1] == control_group)
		col="black";
		sm = ksmooth(dd$time,dd$nI)
		lines(sm$x,sm$y,col=col,lty=lty,lwd=lwd)
		text(xlim[2]-.5,sm$y[length(sm$y)],i);
	}
	legend("topleft",cex=.5,border="n",legend=c(control_group,test_group),col=c("black","red"),lty=1);
}

#plot three figures covering different aspects of variation between animals
ns_plot_pca_diagnostics = function(pca.res,animal_scores,nharm,legend_title,control_group,test_group,analysis_range){
    compare_groups = !is.na(test_group)
    ns_plot_pca_figure_aspect_ratio = 1.41;
	for (i in 1:nharm){
		pca_col = paste0("pca_",i);
		rng = quantile(animal_scores[,pca_col],c(.01,.99))
		wt_scores = animal_scores[animal_scores$group == control_group ,pca_col]
		if (length(wt_scores)==0)
			stop("No control animals were found in animal_scores");
		if(compare_groups){
			mut_scores =  animal_scores[animal_scores$group == test_group,pca_col]
		}else mut_scores = c();
		if (compare_groups && length(mut_scores)==0)
			stop("No test animals were found in animal_scores");
		medians = c(median(wt_scores),median(mut_scores))
		if (compare_groups){
			confint = c(quantile(wt_scores,c(.1,.90)),quantile(mut_scores,c(.1,.9)))
			ns_plot_pca_fd(pca.res,expand=medians,shading=confint,harm=i,colors=c("black","red",rep("black",2),rep("red",2)),ltys=c(1,1,2,2,2,2),shading_col=c(rgb(1,0,0,.05),rgb(0,0,0,.05)))
		}else{
			confint = quantile(wt_scores,c(.1,.90))
			ns_plot_pca_fd(pca.res,expand=medians,shading=confint,harm=i,colors=rep("black",3),ltys=c(1,2,2),shading_col=rgb(0,0,0,.05))
		}
		wt_scores_ecdf = ecdf(wt_scores)
		if(compare_groups){
			mut_scores_ecdf = ecdf(mut_scores)			
		}
		plot(wt_scores_ecdf,xlim=rng,verticals=T,pch=NA,asp=ns_plot_pca_figure_aspect_ratio,xlab="Expansion Magnitude",ylab="cdf")
		if(compare_groups){
			ks = ks.test(wt_scores,mut_scores,"two.sided",exact=T)
			lines(mut_scores_ecdf,col="red",verticals=T,pch=NA)
			legend("bottomright",legend=c(control_group,test_group),title=paste0(legend_title," p=",round(ks[["p.value"]],4)),col=c("black","red"),bty="n",lty=1)
		}
		print("Unusual Animals: ")
		print(animals_to_use[which(pca.res$scores[,i] < rng[1] )])#| pca.res$scores[,i] > rng[2])])

		ns_plot_representative_animals(animal_scores,animals_to_use,data2,harm=i,control_group=control_group,test_group=test_group,xlim=analysis_range);
	}
}

#selects the subset of data that is useful
#returns the data and prints the summary.
ns_summarize_and_select_data_to_use = function(all_data,analysis_range){
	all_data$By.Hand.Movement.Cessation.Relative.Time = as.numeric(all_data$By.Hand.Movement.Cessation.Relative.Time);
	all_data$Best.Guess.Death.Associated.Expansion.Time = as.numeric(all_data$Best.Guess.Death.Associated.Expansion.Time);
	all_data$Best.Guess.Fast.Movement.Cessation.Relative.Time = as.numeric(all_data$Best.Guess.Fast.Movement.Cessation.Relative.Time);
	
#	browser()
	data = subset(all_data,all_data$Excluded %in% c("0","FALSE") & all_data$Number.of.Animals.In.Clump %in% c(0,1) & !is.na(all_data$By.Hand.Movement.Cessation.Relative.Time))
	#find how many animals expanded and how many did not.
	data$expanded = !is.na(data$Best.Guess.Death.Associated.Expansion.Time)
	data$not_expanded = is.na(data$Best.Guess.Death.Associated.Expansion.Time)
	data$animal = paste(data$Plate.Name,data$Group.ID)
	data$By.Hand.Post.mortem.Contraction.Stop.Time = 0

	rr = aggregate(expanded~animal+group,data=data,FUN=sum)
	rr$expanded = rr$expanded > 0
	rr$not_expanded = rr$expanded == 0
	expansion_info = aggregate(expanded~group,data=rr,FUN=sum)
	expansion_info$not_expanded =  aggregate(not_expanded~group,data=rr,FUN=sum)$not_expanded
	expansion_info$frac_expanded = expansion_info$expanded/(expansion_info$expanded + expansion_info$not_expanded)

	data = subset(data, !is.na(data$By.Hand.Post.mortem.Contraction.Stop.Time) & !is.na(data$Best.Guess.Fast.Movement.Cessation.Relative.Time))
	
	rr2 = aggregate(expanded~animal+group,data=data,FUN=sum)
	rr2$expanded = rr2$expanded > 0
	rr2$not_expanded = rr2$expanded == 0
	excluded_info = aggregate(expanded~group,data=rr2,FUN=sum)
	expansion_info$number_after_exclusion=excluded_info$expanded;
	#print(sp)
	print(expansion_info)
	#now we focus on animals whose data fits within the specified ranges
	data = subset(data,!is.na(data$Best.Guess.Death.Associated.Expansion.Time))
	data$time = data$Best.Guess.Death.Associated.Expansion.Time;
	data2 = subset(data,data$time > analysis_range[1] & data$time < analysis_range[2]);

	animals = unique(data2$animal);
	animals_to_use = c();
	data2$nI = data2$Normalized_Stabilitzed_Intensity;
	for (a in animals){
		dd = subset(data2,data2$animal==a)
		#we don't want to inspect times before the animals have stopped fast moving, or after they have stopped contracting
		rng = range(dd$time[dd$Total.Stabilized.Intensity > 0*1000 & dd$By.Hand.Post.mortem.Contraction.Stop.Time <= 0 & dd$Best.Guess.Fast.Movement.Cessation.Relative.Time >= 0])

		if (rng[1] > analysis_range[1]+.25 || rng[2] < analysis_range[2]-.25)
			next;
		#normalize to start
		mm = mean(dd$Total.Stabilized.Intensity[dd$time>-0.25 & dd$time<0]);
		if (!is.finite(mm))
			next;
		data2$nI[data2$animal==a] = data2$Total.Stabilized.Intensity[data2$animal==a]/mm

		animals_to_use = c(animals_to_use,a);
	}
	if (is.null(animals_to_use))
		stop("No animals met the criterea");
	return(list(data=data2,animals_to_use=animals_to_use));
}

#does the spline fitting
ns_fit_splines_and_pca = function(data2,animals_to_use,nbasis_c,smoothing_param,nharm){

		#we build a spline basis to fit all animals
		splbas = create.bspline.basis(rangeval=range(data2$time),norder=3,nbasis=nbasis_c ,breaks=NULL)
		cvec0 <- matrix(0,nbasis_c,1)
		Wfd0<- fd(cvec0, splbas)
		growfdPar2.3 <- fdPar(Wfd0, 1, smoothing_param)

		#we fit all the animals
		animal_fits = list();
		for (i in 1:length(animals_to_use)){
			dd = subset(data2,data2$animal==animals_to_use[i])
			animal_fits[[i]] = smooth.basis(argvals=dd$time,y=dd$nI,returnMatrix=TRUE,fdParobj=growfdPar2.3);
		}


		combined_fits = as.fd(smooth.basis(argvals=
			matrix(c(animal_fits[[1]]$argvals,animal_fits[[2]]$argvals),ncol=2),
			y=matrix(c(animal_fits[[1]]$y,animal_fits[[2]]$y),ncol=2),
			returnMatrix=TRUE,fdParobj=growfdPar2.3));


		for (i in 3:length(animals_to_use)){
			combined_fits_prev = combined_fits;
			combined_fits$coefs = cbind(combined_fits$coefs,as.fd(animal_fits[[i]])$coefs)
		}
		#do the functional PCA
		pca.res = pca.fd(combined_fits,nharm=nharm)
		

		#tabulate the results and write to disk
		data2 = subset(data2,data2$animal %in% animals_to_use)
		data2$pca_1 = 0;
		data2$pca_2 = 0;
		animal_scores = NULL;
		pca.res$group = "";
		for (i in 1:length(animals_to_use)){
			data2$pca_1[data2$animal == animals_to_use[[i]] ] = pca.res$scores[i,1]
			cur_group = as.character(unique(data2$group[data2$animal == animals_to_use[[i]] ]));
			cd = data.frame(animal=animals_to_use[[i]],group=cur_group);
			for (j in 1:nharm){
				cd[,paste0("pca_",j)] = pca.res$scores[i,j];
			}

			animal_scores = rbind(cd,animal_scores);
		}
		return(list(pca.res = pca.res,animal_scores=animal_scores, combined_fits = combined_fits));
}
#depreciated!
if (0){
ns_run_span_regression_models = function(spans2,control_group,test_group){
	N_test = length(spans2$group[spans2$group==test_group])
	N_control = length(spans2$group[spans2$group==control_group])

	spans2 = subset(spans2,spans2$group%in% c(control_group,test_group))
	spans2$duration_changing_posture[spans2$duration_changing_posture<0] = 0;
	spans2$duration_alive_but_nonmoving[spans2$duration_alive_but_nonmoving<0] = 0;
	spans2$group = factor(spans2$group,levels= c(control_group,test_group))

	lm_res = glm(slow_age~group,data=spans2,family=quasipoisson,contrasts=list(w=contr.treatment(2)));
	
	
	lm_res2 = glm(expansion~slow_age*group,data=spans2,family=quasipoisson,contrasts=list(w=contr.treatment(2)));

	print(summary(lm_res))
	mean_slow_age = mean(spans2$slow_age);
	#cross term is the effect of group on the relationship between slow age and posture
	effect_on_frac_slow = exp(summary(lm_res)$coefficients[grep(":",rownames(summary(lm_res)$coefficients)),"Estimate"])
	
	ci = exp(confint(lm_res,level=.95)[grep(":",rownames(confint(lm_res))),]);
	effect_on_frac_slow_u = ci[1]
	effect_on_frac_slow_l = ci[2]
	effect_on_frac_slow_pvalue=summary(lm_res)$coefficients[grep(":",rownames(summary(lm_res)$coefficients)),"Pr(>|t|)"]
	effect_on_abs_slow = exp(contrast(lm_res, list(group = test_group,slow_age=mean_slow_age), list(group = control_group,slow_age=mean_slow_age))$Contrast)
	effect_on_abs_slow_u =  exp(contrast(lm_res, list(group = test_group,slow_age=mean_slow_age), list(group = control_group,slow_age=mean_slow_age))$Upper)
	effect_on_abs_slow_l = exp(contrast(lm_res, list(group = test_group,slow_age=mean_slow_age), list(group = control_group,slow_age=mean_slow_age))$Lower)
	effect_on_abs_slow_pvalue=summary(lm_res)$coefficients[which(substring(rownames(summary(lm_res)$coefficients),1,5)=="group"),"Pr(>|t|)"]


	lm_res = glm(duration_alive_but_nonmoving~group,data=spans2,family=quasipoisson);
	print(summary(lm_res))
	effect_on_alnm = exp(contrast(lm_res, list(group = test_group), list(group = control_group))$Contrast)
	effect_on_alnm_u = exp(contrast(lm_res, list(group = test_group), list(group = control_group))$Upper)
	effect_on_alnm_l = exp(contrast(lm_res, list(group = test_group), list(group = control_group))$Lower)
	effect_on_alnm_pvalue=summary(aov(lm_res))[[1]][["Pr(>F)"]][1]

	res = data.frame(test_group=test_group,
			   control_group=control_group,
			   median_abs_slow_test = median(spans2$slow_age[spans2$group == test_group]),
			   median_abs_slow_control = median(spans2$slow_age[spans2$group == control_group]),
			   median_frac_slow_test = median(spans2$slow_age[spans2$group == test_group]/spans2$expansion[spans2$group == test_group]),
			   median_frac_slow_control = median(spans2$slow_age[spans2$group == control_group]/spans2$expansion[spans2$group == control_group]),
			   median_almn_test = median(spans2$duration_alive_but_nonmoving[spans2$group == test_group]),
			   median_almn_control = median(spans2$duration_alive_but_nonmoving[spans2$group == control_group]),
			   effect_on_frac_slow=effect_on_frac_slow,effect_on_frac_slow_u=effect_on_frac_slow_u,effect_on_frac_slow_l=effect_on_frac_slow_l,
			   effect_on_frac_slow_pvalue=effect_on_frac_slow_pvalue,
			   effect_on_abs_slow=effect_on_abs_slow,effect_on_abs_slow_u=effect_on_abs_slow_u,effect_on_abs_slow_l=effect_on_abs_slow_l,
			   effect_on_abs_slow_pvalue=effect_on_abs_slow_pvalue,
			   effect_on_alnm=effect_on_alnm,effect_on_alnm_u=effect_on_alnm_u,effect_on_alnm_l=effect_on_alnm_l,
			   effect_on_alnm_pvalue=effect_on_alnm_pvalue,
		  		 N_test=N_test,N_control=N_control)
	rownames(res) = NULL;
	return(res);
}
}
ns_calc_simple_pca_stats = function(animal_scores,pca_fit,p){
	grp1 = which(animal_scores$group==p[1])
	if (is.na(p[2])){
		browser()
		return(data.frame(test=NA,
				control=p[1],
				N_test = NA,
				N_control = length(grp1),
				pca_frac_var_1=pca_fit$varprop[1],
				pca_score_diff_1=NA,
				pca_pval=NA	
		))
	}
	
	grp2 = which(animal_scores$group==p[2])
	diff_1 = mean(animal_scores$pca_1[grp2])-mean(animal_scores$pca_1[grp1])
	#diff_2 = abs((mean(animal_scores$pca_2[grp2])-mean(animal_scores$pca_2[grp1]))/mean(animal_scores$pca_2[grp1]))
	ks_res_1 = ks.test(animal_scores$pca_1[grp1],animal_scores$pca_1[grp2],"two.sided",exact=T)
	#ks_res_2 = ks.test(animal_scores$pca_2[grp1],animal_scores$pca_2[grp2],"two.sided",exact=T)
	return(data.frame(test=p[2],
		control=p[1],
		N_test = length(grp2),
		N_control = length(grp1),
		pca_frac_var_1=pca_fit$varprop[1],
		pca_score_diff_1=diff_1,
		pca_pval=ks_res_1$p.value		
	))
}

#depreciated
if (0){
	ns_calc_simple_span_stats = function(spans,p){
		grp1 = which(spans$group==p[1])
		grp2 = which(spans$group==p[2])
		return(data.frame(test=p[2],
			control=p[1],
			N_test = length(grp2),
			N_control = length(grp1),
			slow_test=mean(spans$slow_age[grp2]),
			slow_test_sd=sd(spans$slow_age[grp2]),
			slow_control=mean(spans$slow_age[grp1]),
			slow_control_sd=sd(spans$slow_age[grp1]),
			weak_test=mean(spans$duration_changing_posture[grp2]),
			weak_test_sd=sd(spans$duration_changing_posture[grp2]),
			weak_control=mean(spans$duration_changing_posture[grp1]),
			weak_control_sd=sd(spans$duration_changing_posture[grp1]),
			nonmoving_test=mean(spans$duration_alive_but_nonmoving[grp2]),
			nonmoving_test_sd=sd(spans$duration_alive_but_nonmoving[grp2]),
			nonmoving_control=mean(spans$duration_alive_but_nonmoving[grp1]),
			nonmoving_control_sd=sd(spans$duration_alive_but_nonmoving[grp1]),
			frac_slow_test=mean(spans$slow_age[grp2]/spans$expansion[grp2]),
			frac_slow_test_sd=sd(spans$slow_age[grp2]/spans$expansion[grp2]),
			frac_slow_control=mean(spans$slow_age[grp1]/spans$expansion[grp1]),
			frac_slow_control_sd=sd(spans$slow_age[grp1]/spans$expansion[grp1]),
			frac_nonmoving_test=mean(spans$duration_alive_but_nonmoving[grp2]/spans$expansion[grp2]),
			frac_nonmoving_test_sd=sd(spans$duration_alive_but_nonmoving[grp2]/spans$expansion[grp2]),
			frac_nonmoving_control=mean(spans$duration_alive_but_nonmoving[grp1]/spans$expansion[grp1]),
			frac_nonmoving_control_sd=sd(spans$duration_alive_but_nonmoving[grp1]/spans$expansion[grp1])
		))
	}
}

