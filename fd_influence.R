# -----------------------------------------------------------
# Influence of each species on overall Petchey Gaston FD: 
# The change in total dendrogram branch length after a species is removed
# -----------------------------------------------------------

FD.influence <- function(traitfile, scale = FALSE, traittype = NULL, dis.metric = "euclid", clus.method = "average"){
  require(reshape)
  traits <- traitfile
  
  # Create 'communities' where all species are present but one
  sp = rownames(traits)
  
  com <- melt(sapply(sp, function(x) sp[sp!=x]))[,c(3,2)]
  names(com) <- c("species", "community")
  
  distances = daisy(traits, metric = dis.metric, stand = FALSE)
  
  if(length(distances[is.na(distances)]) > 1) { NA; warning("incomparable species (check trait file)")}
  
  else { tree = hclust(distances, method = clus.method)
  xtree = Xtree(tree)
  
  FD <- tapply(com$species, com$community,
               function(x) Getlength(list(xtree[[1]], 
                                          xtree[[2]][!is.na(match(toupper(dimnames(xtree[[2]])[[1]]),toupper(x))),])))
  
  # Full dendrogram
  full <- Getlength(list(xtree[[1]], xtree[[2]]))
  
  
  absdiff <- full - FD
  reldiff <- 100*(full - FD ) / full
  }
  
  reldiff[rownames(traitfile)]
}
