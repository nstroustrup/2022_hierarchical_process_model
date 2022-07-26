#estimates dS and dR of the regression model
#log(yi) = Beta_T*T + Beta_X*X + e
#where yi is the event time from the set Y = {V,D} 
#T specifies the type of event (eg does yi describe an event from V or D)
#I specifies the intervention, eg whether y is observed in a wild-type population or a mutant population
#In this case, Beta_t is our best esimate of R_f and Beta_x is our best estimate of Delta_Drift_s.
ns_rms_contrasts_to_data_frame_2 = function(rms_data,group_names){
  resp = cbind(rms_data$Contrast,rms_data$SE,rms_data$Lower,rms_data$Upper,rms_data$Z,rms_data$Pvalue)
  colnames(resp) = c("coefficient", "se", "lower", "upper", "Z", "p")
  resp = as.data.frame(resp)
  resp$group = group_names;
  return(resp);
}

ns_estimate_dS_and_dR = function(death_times,group_column="group",event_columns=c("slow_age","expansion"),device_column=NULL,ref_group=NULL){
  if (length(event_columns) != 2)
    stop("Provide exactly two event columns")
  df = melt(death_times,id.vars=c(group_column,device_column),measure.vars=event_columns)
  df$c = 1;
  names(df)[names(df)==device_column] = "device";
  names(df)[names(df)==group_column] = "group";
  
  data_dist <<- datadist(df);
  options(datadist="data_dist")
  
  if (is.null(device_column)){
    res = bj(Surv(value,c)~variable*group, time.inc=.0001, data=df,control=list(iter.max=60), 
             x=TRUE, y=TRUE)
  }else{
    res = bj(Surv(value,c)~variable*group+device, time.inc=.0001, data=df,control=list(iter.max=60), 
             x=TRUE, y=TRUE)
  }
  if (is.null(ref_group))
    ref_group = unique(df$group)[1];
  
  var_levels = unique(df$variable);
  group_levels = unique(df$group)
  group_levels = group_levels[group_levels != ref_group];
  effect_of_mut_on_both= NULL;
  for (g in rev(unique(group_levels))){
    effect_of_mut_on_both = rbind(ns_rms_contrasts_to_data_frame_2(rms::contrast(res,
                                                                                 list(group=ref_group,variable=var_levels),
                                                                                 list(group=g,variable=var_levels),type="average",usebootcoef=TRUE),g),effect_of_mut_on_both)
  }
  effect_of_mut_on_slow = ns_rms_contrasts_to_data_frame_2(rms::contrast(res,
                                                                         list(group=rep(ref_group,length(group_levels)),variable=event_columns[1]),
                                                                         list(group=group_levels,variable=event_columns[1]),usebootcoef=TRUE),group_levels)
  effect_of_mut_on_exp= ns_rms_contrasts_to_data_frame_2(rms::contrast(res,
                                                                       list(group=rep(ref_group,length(group_levels)),variable=event_columns[2]),
                                                                       list(group=group_levels,variable=event_columns[2]),usebootcoef=TRUE),group_levels)
  
  
  cols = c("coefficient","lower","upper","group")
  est_dS = effect_of_mut_on_both[cols]
  est_dS$var = "dS"
  est_dR_C = effect_of_mut_on_slow$coefficient-effect_of_mut_on_exp$coefficient;
  est_dR_L = effect_of_mut_on_slow$lower-effect_of_mut_on_exp$upper;
  est_dR_U = effect_of_mut_on_slow$upper-effect_of_mut_on_exp$lower
  est_dR = data.frame(coefficient=est_dR_C,lower=est_dR_L,upper=est_dR_U,group=effect_of_mut_on_slow$group)
  est_dR$var = "dR"
  cmb = rbind(est_dS,est_dR);
  return(list(coefficients=cmb,model=res))
}


if(0){
  #run example
  #one intervention
  dS1 = 2;
  dR1 = .9
  
  
  #another intervention
  dS2 = 1.75;
  dR2 = 1.1
  
  dD1 = dS1/sqrt(dR1)
  dV1 = dS1*sqrt(dR1)
  dD2 = dS2/sqrt(dR2)
  dV2 = dS2*sqrt(dR2)
  device_effect = .8;
  N=50;
  death_times = data.frame(slow_age=c(rinvgauss(N,8,1000),rinvgauss(N,8/dV1,1000),rinvgauss(N,8/dV2,1000)),
                           expansion= c(rinvgauss(N,16,1000),rinvgauss(N,16/dD1,1000),rinvgauss(N,16/dD2,1000)),
                           group = c(rep("wt",N),rep("m1",N),rep("m2",N)),
                           device = sample(c("a","b"),3*N,replace=T)
  )
  death_times$slow_age[death_times$device=="a"] = death_times$slow_age[death_times$device=="a"]*device_effect;
  
  plot(survfit(Surv(slow_age)~group,death_times)	,xlim=c(0,25))
  lines(survfit(Surv(expansion)~group,death_times) , col="red")
  result = ns_calc_dS_and_dR(death_times,device_column="device",ref_group = "wt");
  print(cbind(exp(result$coefficients[,c("coefficient","lower","upper")]),result$coefficient[,c("group","var")]))
}