scan_survival_kulldorff = function(Lon,Lat,time,
                                   cens,simulations = 499,
                                   scan_prop = 0.3,type = "both"){
  
  if(class(Lon)!="numeric"){stop("Lon must be numeric",call. = F)}
  if(class(Lat)!="numeric"){stop("Lat must be numeric",call. = F)}
  if((class(time)!="numeric")|(sum(time<0)!=0)){stop("time must be numeric and greater than 0",call. = F)}
  if(length(which(cens==1|cens==0))!=length(cens)){stop("cens must be 0 or 1",call. = F)}
  if(!all(sapply(list(length(Lon),length(Lat)),function(x){identical(x,length(time))}))){
    stop("Lon, Lat and time must be the same length",call. = F)
  }
  if(class(simulations)!="numeric"|class(scan_prop)!="numeric"){
    stop("simulations and scan_prop must be numeric",call. = F)}
  if((scan_prop<=0|scan_prop>0.5)){stop("scan must be greater than 0 and less than 0.5",call. = F)}
  if(simulations<1|(simulations%%1)!=0){
    stop("simulations must be greater than 0 and not decimal",call. = F)}
  if(!(type%in%c("low","both","high"))){stop("type must be low, both or high",call. = F)}
  
  time = (time-min(time))/(max(time)-min(time))
  
  database          = data.frame(Lon,Lat,time,cens)
  distance          = as.matrix(dist(database[,c("Lon","Lat")], method = "euclidian",diag = T,upper = T))
  dist_list         = split(distance, seq(nrow(distance)))
  dist_order_index  = lapply(dist_list,  function(x){sort(x,index.return=T)$ix})
  dist_order_value  = lapply(dist_list,  function(x){sort(x,index.return=T)$x})
  dist_order_value_unique = lapply(dist_order_value, unique)
  
  obs_max = round(scan_prop*dim(distance)[1])
  
  message("finding the zones")
  
  zones = list()
  for(i in 1:dim(distance)[1]){
    j = 1
    positions = dist_order_index[[i]][which(dist_order_value[[i]]%in%dist_order_value_unique[[i]][1:(j+1)])]
    while(length(positions)<=obs_max){
      vetor            = rep(0,dim(distance)[1])
      vetor[positions] = 1 
      zones = append(zones,list(vetor))
      j = j+1
      positions = dist_order_index[[i]][which(dist_order_value[[i]]%in%dist_order_value_unique[[i]][1:(j+1)])]
    }
  }
  zones = unique(zones)
  
  r_in  = lapply(zones,function(x){sum(x*database$cens)})          
  r_out = lapply(zones,function(x){sum(((x-1)*-1)*database$cens)})       
  t_in  = lapply(zones,function(x){sum(x*database$time)}) 
  t_out = lapply(zones,function(x){sum(((x-1)*-1)*database$time)})  
  
  indicator = rep(1,length(zones))
  if(type!="both"){
    if(type=="low"){
      div_in = mapply(function(x,y){x/y},x=t_in,y=r_in)
      div_out = mapply(function(x,y){x/y},x=t_out,y=r_out)
      indicator[which(div_in>div_out)]=0
      zones = zones[c(which(indicator==1))]
      r_in  = lapply(zones,function(x){sum(x*database$cens)})          
      r_out = lapply(zones,function(x){sum(((x-1)*-1)*database$cens)})       
      t_in  = lapply(zones,function(x){sum(x*database$time)}) 
      t_out = lapply(zones,function(x){sum(((x-1)*-1)*database$time)})
    }else{
      div_in = mapply(function(x,y){x/y},x=t_in,y=r_in)
      div_out = mapply(function(x,y){x/y},x=t_out,y=r_out)
      indicator[which(div_in<div_out)]=0
      zones = zones[c(which(indicator==1))]
      r_in  = lapply(zones,function(x){sum(x*database$cens)})          
      r_out = lapply(zones,function(x){sum(((x-1)*-1)*database$cens)})       
      t_in  = lapply(zones,function(x){sum(x*database$time)}) 
      t_out = lapply(zones,function(x){sum(((x-1)*-1)*database$time)})
    }
  }
  
  R = sum(database$cens)
  t_total = sum(database$time)
  # calculating lambda
  lambda = mapply(function(r_in,r_out,t_in,t_out,R,t_total){
    (((r_in/(t_in))^r_in)*((r_out/(t_out))^r_out))/((R/(t_total))^R)
  },r_in = r_in, r_out = r_out, t_in = t_in, t_out = t_out, R = R, t_total = t_total)
  
  cluster_very_prob = zones[[which.max(lambda)]]
  value_cluster     = max(lambda)
  positions_cluster_very_prob = which(cluster_very_prob==1)
  
  value_cluster_random = c()
  for(i in 1:simulations){
    linhas = sample(1:dim(database)[1])
    cens_sample = database$cens[linhas]
    time_sample = database$time[linhas]
    
    r_in  = lapply(zones,function(x){sum(x*cens_sample)})          
    r_out = lapply(zones,function(x){sum(((x-1)*-1)*cens_sample)})       
    t_in  = lapply(zones,function(x){sum(x*time_sample)}) 
    t_out = lapply(zones,function(x){sum(((x-1)*-1)*time_sample)})
    
    lambda = mapply(function(r_in,r_out,t_in,t_out,R,t_total){
      (((r_in/(t_in))^r_in)*((r_out/(t_out))^r_out))/((R/(t_total))^R)
    },r_in = r_in, r_out = r_out, t_in = t_in, t_out = t_out, R = R, t_total = t_total)
    
    value_cluster_random  = c(value_cluster_random,max(lambda))
    message("simulation number ",i)
  }
  
  w = (sum(value_cluster_random<value_cluster))
  p_value = 1-(w/(simulations+1))
  output = list(positions = positions_cluster_very_prob,
                LLR       = value_cluster,
                p_value   = p_value)
  output
}
