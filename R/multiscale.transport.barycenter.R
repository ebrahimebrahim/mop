


multiscale.transport.barycenter <- function( X.mean, weights.mean = 1, 
                                             gmra.marginals, weights.marginals,
                                             trp.lp, p=2, nType=0, dType=1, scaleMass = TRUE, 
                                             mean.radius=10^-10, step=0.5, n.iterations=100, 
                                             rel.err=0.00001, min.step=0, ret.sol=TRUE, 
                                             store.sol = FALSE, store.prefix="./trp-", verbose=TRUE){
 
  gmra.mean = gmra.create.ikm(X=X.mean, eps=mean.radius, nKids=4, stop=3)

  if( ret.sol ){
    sol <- list() 
  }
  else{
    sol <-  NULL
  }
  prev.cost = Inf 
  cost = .Machine$double.xmax
  costs <- c()
  avg.steps  <- c()
  iteration = 0
  X.min = X.mean
  min.cost = Inf
  avg.step = Inf
  while( cost < prev.cost & 
         abs(cost - prev.cost) / cost > rel.err & 
         iteration < n.iterations & 
         avg.step > min.step){ 
    iteration = iteration+1
    prev.cost <- cost

    X.update <- as.data.table( matrix(0, nrow = nrow(X.mean), ncol = ncol(X.mean) ) )
    mass.update <- rep(0, nrow(X.mean) )
    
    cost <- 0

    for( i in 1:length(gmra.marginals) ){
      trp <- multiscale.transport.solve( trp.lp, gmra.mean, gmra.marginals[[i]], 
                                            p = p, nType=nType, 
                                            dType=dType, scaleMass=scaleMass, 
                                            w1=weights.mean, w2 = weights.marginals[[i]] )
      if( ret.sol ){
         sol[[i]] = trp
      }
      if( store.sol ){
        save( trp, file=sprintf("%s%d.Rdata", store.prefix, i) )
      }

      n = length(trp$cost)
      cost <- cost + trp$cost[[n]]
      if(trp$cost[[n]] > 0 ){
    
        map = trp$map[[n]]
        from.index = trp$fromIndex[[n]]
        from.size = trp$fromSize[[n]] 
       
        delta = as.data.table( trp$to[[n]][map[,2], ] - trp$from[[n]][map[,1], ] )
        delta = delta * map[,3]
        delta$from.id = map[,1]
        delta = delta[, lapply(.SD, sum, na.rm=TRUE), by = from.id]
        delta = delta[order(from.id)]

        delta.mass = data.table( from.id = map[,1], mass=map[,3] )
        delta.mass = delta.mass[, lapply(.SD, sum, na.rm=TRUE), by = from.id]
        delta.mass = delta.mass[order(from.id)]   

        delta = delta[,!"from.id"]

        delta = delta[ rep(1:length(from.size), from.size), ]
        delta.mass = delta.mass[ rep(1:length(from.size), from.size), ]

        
        X.update[ from.index, ]  = X.update[ from.index, ]  + delta
        mass.update[ from.index ] = mass.update[ from.index ] + delta.mass$mass

      }
    }
    
    costs = c(costs, cost)
    if( prev.cost < cost ){
      step = 0.9 * step
    }
    #if( min.cost < cost ){
      #X.mean = prev.X.mean
     # step = 0.9 * step
    #  X.mean = X.min
    #  prev.cost=cost
    #  cost = min.cost
    #}
    #else{
    #  X.min = X.mean
      X.update = X.update / mass.update
      X.update = step * X.update #/ length(gmra.marginals)
      avg.step = mean( sqrt(rowSums(X.update^2) ) )
      X.mean = X.mean + X.update
      avg.steps = c(avg.steps, avg.step) 
    #}
    gmra.delete(gmra.mean)
    gmra.mean = gmra.create.ikm(X.mean, eps=mean.radius, nKids = 4, stop=3)
    #plot( X.mean, col="#FF000099")
    
    if(verbose){
       print( sprintf( "Iteration : %i , cost=%f, step=%f",  iteration, cost, avg.step) )
    }
   # min.cost = min(min.cost, cost)
  }
 
  list( trp.plans = sol, gmra.mean = gmra.mean, X.mean = X.mean, costs=costs, avg.steps=avg.steps) 
}






barycenter.parallel.update <- function( a, b ){
  a[[1]][ b[[4]], ]  = a[[1]][ b[[4]], ]  + b[[1]]
  a[[2]][ b[[4]] ] = a[[2]][ b[[4]] ] + b[[2]]
  a[[3]] = a[[3]] + b[[3]]
  a
}

multiscale.transport.barycenter.parallel <- function( X.mean, gmra.mean.file, weights.mean = 1, 
                                             gmra.files, weights.marginals, trp.files,
                                             p=2, nType=0, dType=1, scaleMass = TRUE, 
                                             mean.radius=10^-10, step=0.5, n.iterations=100, 
                                             rel.err=0.00001, min.step=0, verbose=TRUE, 
                                             transport.type=0, massCost=0, n.parallel=-1){
  library(foreach, quietly=TRUE)
  library(doParallel, quietly=TRUE)

 
  
  gmra.tmp = gmra.create.ikm(X=X.mean, eps=mean.radius, nKids=4, stop=3)
  gmra.save.tree( gmra.tmp, gmra.mean.file)
  gmra.delete( gmra.tmp )

  if(verbose){
    print("GMRA computed")
  }

  prev.cost = Inf 
  cost = .Machine$double.xmax
  costs <- c()
  avg.steps  <- c()
  iteration = 0
  min.cost = Inf
  avg.step = Inf

  while( cost < prev.cost &
         abs(cost - prev.cost) / cost > rel.err & 
         iteration < n.iterations & 
         avg.step > min.step){ 
    iteration = iteration+1
    prev.cost <- cost
 
    if(n.parallel < 1){
      n.parallel <- detectCores() 
    } else{
      n.parallel <- min( n.parallel, detectCores() )
    }

    cl <- makeCluster( n.parallel )
    registerDoParallel( cl )
    #registerDoSEQ()

    update <- foreach( i = 1:length(gmra.files), .combine = barycenter.parallel.update  ) %dopar% {
      gmra.mean = gmra.load.tree( gmra.mean.file)
      gmra.marginal = gmra.load.tree( gmra.files[[i]] )


      trp.lp <- multiscale.transport.create.lp(oType=31, transport.type=transport.type, massCost=massCost )
      icprop <- multiscale.transport.create.iterated.capacity.propagation.strategy(3, 0)
      #multiscale.transport.add.expand.neighborhood.strategy(trp.lp, 2)
      multiscale.transport.set.propagation.strategy.1(trp.lp, icprop);

      trp <- multiscale.transport.solve( trp.lp, gmra.mean, gmra.marginal, 
                                         p = p, nType=nType, dType = dType, scaleMass = scaleMass, 
                                         w1 = weights.mean, w2 = weights.marginals[[i]] )
      
      save( trp, file = trp.files[[i]] )

      n = length(trp$cost)
    
      map = trp$map[[n]]
      from.index = trp$fromIndex[[n]]
      from.size = trp$fromSize[[n]]

       
      delta = as.data.table( trp$to[[n]][map[,2], ] - trp$from[[n]][map[,1], ] )
      delta = delta * map[,3]
      delta$from.id = map[,1]
      delta = delta[, lapply(.SD, sum, na.rm=TRUE), by = from.id]
      delta = delta[order(from.id)]
      if( nrow(delta) < length(from.size) ){
        delta.tmp <- data.table( matrix(0, nrow=length(from.size), ncol=ncol(delta)-1 ) )
        delta.tmp[delta$from.id, ] <- delta[, !"from.id"]
        delta <- delta.tmp
      }
      else{
         delta <- delta[,!"from.id"]
      }

      delta.mass = data.table( from.id = map[,1], mass=map[,3] )
      delta.mass = delta.mass[, lapply(.SD, sum, na.rm=TRUE), by = from.id]
      delta.mass = delta.mass[order(from.id)] 
      if( nrow(delta.mass) < length(from.size) ){
        delta.mass.tmp <- rep(0, length(from.size) )
        delta.mass.tmp[ delta.mass$from.id ] <- delta.mass$mass
        delta.mass <- delta.mass.tmp
      }
      else{
        delta.mass = delta.mass$mass
      }

      delta = delta[ rep(1:length(from.size), from.size), ]
      delta.mass = delta.mass[ rep(1:length(from.size), from.size) ]

      res <- list()
      res[[1]] = delta
      res[[2]] = delta.mass
      res[[3]] = trp$cost[[n]] 
      res[[4]] = from.index
      #X.update[ from.index, ]  = X.update[ from.index, ]  + delta
      #mass.update[ from.index ] = mass.update[ from.index ] + delta.mass$mass

      multiscale.transport.delete.lp(trp.lp);
      gmra.delete( gmra.mean )
      gmra.delete( gmra.marginal )
      res 
    }
    


    X.update <- update[[1]]
    mass.update <- update[[2]]
    cost <- update[[3]]
    
    costs = c(costs, cost)
    if( prev.cost < cost ){
      step = 0.9 * step
    }
    
    X.update = X.update / mass.update
    X.update[mass.update==0, ] = 0
    X.update = step * X.update # / length(gmra.files)
    avg.step = mean( sqrt(rowSums(X.update^2) ) )
    X.mean = X.mean + X.update
    avg.steps = c(avg.steps, avg.step) 
    
    gmra.tmp = gmra.create.ikm(X=X.mean, eps=mean.radius, nKids=4, stop=3)
    gmra.save.tree( gmra.tmp, gmra.mean.file)
    gmra.delete( gmra.tmp )

    update <- NULL
    if(verbose){
       print( sprintf( "Iteration : %i , cost=%f, step=%f, step.factor=%f",  iteration, cost, avg.step, step) )
    }

    stopCluster(cl)
    #plot( X.mean, col="#FF000099")
    
  }

  
  
  list( X.mean = X.mean, costs=costs, avg.steps=avg.steps, iterations=iteration) 
}

