scan_normal_kulldorff = function(Lon,Lat,variable,simulations = 499,
                                 scan_prop = 0.3,type = "both"){
  
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
  N       = dim(database)[1] 
  lambda_z = mapply(function(x_z,n_z){(X-x_z)/(N-n_z)},x_z=x_z,n_z=n_z) 
  
  indicator = rep(1,length(zones))
  if(type!="both"){
    if(type=="low"){
      indicator[which(mu_z>lambda_z)]=0
      zones = zones[c(which(indicator==1))]
      x_z  = lapply(zones,function(x){sum(x*database$variable)})
      n_z  = lapply(zones,function(x){length(which(x==1))})
      mu_z = mapply(function(x,y){x/y},x=x_z,y=n_z)
      lambda_z = mapply(function(x_z,n_z){(X-x_z)/(N-n_z)},x_z=x_z,n_z=n_z)
    }else{
      indicator[which(mu_z<lambda_z)]=0
      zones = zones[c(which(indicator==1))]
      x_z  = lapply(zones,function(x){sum(x*database$variable)})
      n_z  = lapply(zones,function(x){length(which(x==1))})
      mu_z = mapply(function(x,y){x/y},x=x_z,y=n_z)
      lambda_z = mapply(function(x_z,n_z){(X-x_z)/(N-n_z)},x_z=x_z,n_z=n_z)
    }
  }
  
  # calculating the LLR (Z)
  mu      = mean(database$variable)
  sigma2 = sum((database$variable-mu)^2)/N
  
  sigma2_z  = (1/N)*(unlist(lapply(zones,function(x){sum((x*database$variable)^2)})) -
                       (2*unlist(x_z)*mu_z) + unlist(n_z)*(mu_z^2) + 
                       unlist(lapply(zones,function(x){sum((((x-1)*-1)*database$variable)^2)})) -
                       2*(X-unlist(x_z))*lambda_z + (N-unlist(n_z))*(lambda_z^2))
  
  ln_lz = (N*log(sqrt(sigma2)) + sum(sapply(database$variable,function(x){((x-mu)^2)/(2*sigma2)})) -
             (N/2) - N*log(sqrt(sigma2_z)))
  
  cluster_very_prob = zones[[which.max(ln_lz)]]
  value_cluster     = max(ln_lz)
  positions_cluster_very_prob = which(cluster_very_prob==1)
  
  value_cluster_random = c()
  for(i in 1:simulations){
    random_values = sample(database$variable)
    
    x_z  = lapply(zones,function(x){sum(x*random_values)})          
    n_z  = lapply(zones,function(x){length(which(x==1))})         
    mu_z = mapply(function(x,y){x/y},x=x_z,y=n_z)                   
    
    lambda_z = mapply(function(x_z,n_z){(X-x_z)/(N-n_z)},x_z=x_z,n_z=n_z) 
    
    sigma2_z  = (1/N)*(unlist(lapply(zones,function(x){sum((x*random_values)^2)})) -
                         (2*unlist(x_z)*mu_z) + unlist(n_z)*(mu_z^2) + 
                         unlist(lapply(zones,function(x){sum((((x-1)*-1)*random_values)^2)})) -
                         2*(X-unlist(x_z))*lambda_z + (N-unlist(n_z))*(lambda_z^2))
    
    ln_lz = (N*log(sqrt(sigma2)) + sum(sapply(random_values,function(x){((x-mu)^2)/(2*sigma2)})) -
               (N/2) - N*log(sqrt(sigma2_z)))
    
    value_cluster_random  = c(value_cluster_random,max(ln_lz))
    message("simulation number ",i)
  }
  
  w = (sum(value_cluster_random<value_cluster))
  p_value = 1-(w/(simulations+1))
  output = list(positions = positions_cluster_very_prob,
                LLR       = value_cluster,
                p_value   = p_value)
  output
}
