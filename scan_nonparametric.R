scan_nonparametric = function(Lon,Lat,variable,simulations,
                              scan_prop = 0.3,type = "high",exact=F){
  
  if(class(Lon)!="numeric"){stop("Lon must be numeric",call. = F)}
  if(class(Lat)!="numeric"){stop("Lat must be numeric",call. = F)}
  if(class(variable)!="numeric"){stop("variable must be numeric",call. = F)}
  if(!all(sapply(list(length(Lon),length(Lat)),function(x){identical(x,length(variable))}))){
    stop("Lon, Lat and variable must be the same length",call. = F)
  }
  if(class(simulations)!="numeric"|class(scan_prop)!="numeric"){
    stop("simulations and scan_prop must be numeric",call. = F)}
  if((scan_prop<=0|scan_prop>0.5)){stop("scan must be greater than 0 and less than 0.5",call. = F)}
  if(simulations<1|(simulations%%1)!=0){
    stop("simulations must be greater than 0 and not decimal",call. = F)}
  if(!(type%in%c("low","both","high"))){stop("type must be low, both or high",call. = F)}
  if(class(exact)!="logical"){stop("exact must be logical",call. = F)}
  
  database = data.frame(Lon,Lat,variable)
  distance = as.matrix(dist(database[,c("Lon","Lat")], method = "euclidian",diag = T,upper = T))
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
  
  x_z  = lapply(zones,function(x){sum(x*database$variable)})          
  n_z  = lapply(zones,function(x){length(which(x==1))})       
  mu_z = mapply(function(x,y){x/y},x=x_z,y=n_z) 
  X    = sum(database$variable)
  N    = dim(database)[1] 
  lambda_z = mapply(function(x_z,n_z){(X-x_z)/(N-n_z)},x_z=x_z,n_z=n_z) 
  
  indicator = rep(1,length(zones))
  if(type!="both"){
    if(type=="low"){
      indicator[which(mu_z>lambda_z)]=0
      zones = zones[c(which(indicator==1))]
    }else{
      indicator[which(mu_z<lambda_z)]=0
      zones = zones[c(which(indicator==1))]
    }
  }
  
  # calculating the p-value
  x_z  = lapply(zones,function(x){x*database$variable})          
  x_z  = lapply(x_z,function(x){x[which(x!=0)]})
  x_c  = lapply(zones,function(x){((x-1)*-1)*database$variable})
  x_c  = lapply(x_c,function(x){x[which(x!=0)]})
  
  if(exact==F){
    p_value_wilcox = mapply(function(x,y){wilcox.test(x,y)$p.value},x=x_z,y=x_c)
  }else{
    p_value_wilcox = mapply(function(x,y){
      if((length(x)<10)|(length(y)<10)){
        wilcox.test(x,y,exact = T)$p.value
      }else{wilcox.test(x,y)$p.value}
    },x=x_z,y=x_c)
  }
  cluster_very_prob = zones[[which.min(p_value_wilcox)]]
  value_cluster     = min(p_value_wilcox)
  positions_cluster_very_prob = which(cluster_very_prob==1)
  
  value_cluster_random = c()
  if(exact==F){
    for(i in 1:simulations){
      random_values = sample(database$variable)
      
      x_z  = lapply(zones,function(x){x*random_values})          
      x_z  = lapply(x_z,function(x){x[which(x!=0)]})
      x_c  = lapply(zones,function(x){((x-1)*-1)*random_values})
      x_c  = lapply(x_c,function(x){x[which(x!=0)]})
      
      p_value_wilcox = mapply(function(x,y){wilcox.test(x,y)$p.value},x=x_z,y=x_c)
      
      value_cluster_random  = c(value_cluster_random,min(p_value_wilcox))
      message("simulation number ",i)
    }
  }else{
    for(i in 1:simulations){
      random_values = sample(database$variable)
      
      x_z  = lapply(zones,function(x){x*random_values})          
      x_z  = lapply(x_z,function(x){x[which(x!=0)]})
      x_c  = lapply(zones,function(x){((x-1)*-1)*random_values})
      x_c  = lapply(x_c,function(x){x[which(x!=0)]})
      
      p_value_wilcox = mapply(function(x,y){
        if((length(x)<10)|(length(y)<10)){
          wilcox.test(x,y,exact = T)$p.value
        }else{wilcox.test(x,y)$p.value}
      },x=x_z,y=x_c)
      
      value_cluster_random  = c(value_cluster_random,min(p_value_wilcox))
      message("simulation number ",i)
    }
  }
  w = (sum(value_cluster_random>value_cluster))
  p_value = 1-(w/(simulations+1))
  output = list(positions = positions_cluster_very_prob,
                p_value   = p_value)
  output
}
