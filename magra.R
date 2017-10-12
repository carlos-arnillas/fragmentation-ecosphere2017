# functions to convert a map into a graph
# and to plot the patch-size profile and the connectivity profile
require(ggplot2)
require(igraph)
require(data.table)
require(raster)

# Receives a raster map and creates the tables needed to build a graph:
#   Identify the fragments using the Queen's rule.
#   Identify the edges of each fragment
#   Measure the distance among each fragment and any fragment with larger ID
#   It returns a list containing: a data.frame of edges connecting the fragments with the distance and identity of the nodes
#                                 a data.frame containing the fragments identity and size (areas)

class.distances <- function(map, redo.clumps=FALSE, areas=NULL, mem.thre=1e7, echo=TRUE) {
  if (!class(map)=="RasterLayer") stop("Map object has to be a raster.")
  # making the fragments
  if (redo.clumps) {
    map[!is.na(map[])] <- 1
    mapC <- clump(map, directions=8)
  } else {
    mapC <- map
  }
  # getting the areas
  plonglat <- isLonLat(mapC)
  if (is.null(areas)) if (plonglat) areas <- area(mapC)
  if (is.null(areas)) dA <- data.table(freq(mapC, useNA="no")) else dA <- data.table(zonal(areas, mapC, 'sum'))
  setnames(dA,c("fragment","area"))
  # now the edges... using boundaries only
  # get the boundaries of the fragments
  fragments <- mapC * boundaries(mapC, type='inner', asNA=T)
  # get the fragments' perimeter points
  dF <- data.table(as.data.frame(fragments, xy=TRUE, na.rm=TRUE))
  rm(fragments)
  setnames(dF, c("x","y","P"))
  # organize the fragments from smaller to larger to prevent memory issues
  lfragments <- dF[,.(n=.N),by=P][order(n)]
  setkey(dF,P)
  setkey(lfragments,n)
  # for each fragment, measure the distances to the other fragments (only those with larger P code)
  dD <- data.frame(from = integer(0), to = integer(0), dist = numeric(0))
  mlp <- nrow(lfragments); pm <- 0; gc()
  for (j in 1:(mlp-1)) {
    # counter
    if (floor(j/mlp*20) != pm) {
      if (echo) {cat("+");pm <- pm+1}  
      gc()
    }
    if (echo==2) cat(j,"\n") else if (echo) cat(c("-","\\","|","/")[j%%4+1])
    # get the points of the polygon origin (p1) and every potential destiny fragment
    p1 <- dF[.(lfragments[j,P])]
    p2 <- dF[.(lfragments[(j+1):mlp,P])]
    # get the distances from p1 to every fragment in p2, if there is something in p2
    if (nrow(p2)>=1) {
      # cut the problem in small sections to be sure that memory won't collapse.
      maxP2 <- floor(mem.thre/nrow(p1))
      p2[, sec := 1:nrow(p2) %/% maxP2]
      # in each small section
      for (ii in unique(p2$sec)) {
        # counter
        if (echo==3 & ii == 1) cat(":",p1$P[1],":",sep="")
        if (echo) cat(".")
        # get the distances
        d1 <- pointDistance(p1[,.(x,y)], p2[sec==ii,.(x,y)], lonlat=plonglat, allpairs = T)
        if (echo) cat("\b")
        # get the distances from the fragment to each point
        if (nrow(p1) == 1) {
          p2[, dist := d1]
        } else {
          if (nrow(p2) == 1) {
            p2[, dist := min(d1)]
          } else {
            if (is.null(dim(d1))) browser()
            p2[sec==ii, dist := apply(d1, 2, min)]  
          }
        }
      }
      # and look for the closest pixel
      db0 <- p2[, .(dist=min(dist)), by=.(to=P)][order(to)]
      db1 <- db0[, .(from=lfragments[j,P], to=to, dist=dist)]
      dD <- rbind(dD, db1)
      # delete unneeded variables
      rm(p1,p2,d1, db0, db1)
    }
    if (echo) cat("\b")
  }
  # return the results
  return(list(edges=unique(dD), vertices=unique(dA)))
}

# convert a class i of a map stored in the variable region into a graph
# The algorithm calls class.distance
# 
# variables region, map0, areas are transfered as text just to avoid creating huge variables everytime
#
#   i: value of the class in the map to work with
#   fprefix: a prefix added to the output file
#   echo: print extra information of the advances of the function
#   region: name of the variable containing the map in the environment that calls the function
#   map0:   an empty variable used to speed up map operations
#   areas:  the area of each pixel
#   dirExp: folder to store the results
#   mem.thre: a value to prevent too much memory being used in the calculations.
map.distances <-function (i, fprefix="", echo=T, region="region", map0="map0", 
                          areas="areas", dirExp="", mem.thre=1e7) {
  # linking the variables
  if (is.character(region)) {
    region <- get0("region",-2)
    if (is.null(region)) stop("Region map doesn't exist.\n")
  }
  
  if (is.character(map0)) {
    map0 <- get0("map0",-2)
    if (is.null(map0)) {
      map0 <- map0
    }
  }
  
  if (is.character(areas)) {
    areas <- get0("areas",-2)
    if (is.null(areas)) {
      areas <- area(region)
    }
  }
  # setting some variables
  plonglat <- grepl("proj=longlat", as.character(crs(region)))
  map0[] <- NA_integer_
  
  if (echo) cat("Working in biome", i, ": ")
  # build the fragments
  map0[region == i] <- i
  cD <- class.distances(map0, redo.clumps=FALSE, areas=NULL, mem.thre=1e7)
  dD <- cbind(biome=i, cD$dD); dA <- cbind(biome=i, cD$dA)
  
  ###########################
  # if possible: export the graph
  try(save(dD, dA, file=paste0(dirExp, "toGraph.", if (fprefix=="") i else paste0(fprefix,".",i),".RData")), T)
  if (echo) cat("\n")
  # return the edges and vertices as data.frames
  return(list(edges=unique(dD), vertices=unique(dA)))
}

# read edges and vertices from files created with the previous function map.distances
read.graph <- function(fname) {
  err1 <- try(load(file=fname), silent=T)
  if (class(err1)!="try-error") return(list(edges=unique(dD), vertices=unique(dA))) else return (NULL)
}

# apply a filter to the graph, removing vertices and edges consistenlty
criteria <- function(edges0, vertices0, v.criteria=NULL, e.criteria=NULL) {
  vertices <- data.table(vertices0)
  edges <- data.table(edges0)
  if (!is.null(v.criteria)) vertices <- vertices[eval(v.criteria)]
  if (!is.null(e.criteria)) edges <- edges[eval(e.criteria)]
  edges <- edges[(edges$from %in% vertices$fragment) & (edges$to %in% vertices$fragment)]
  return(list(edges=edges, vertices=vertices))
}

# creating the profile distances
# for a given graph stored in db.dist, estimates the R-values of LEM-communities characterized by:
#   p: the parameter controlling the turnover rate (0: no turnover; 1: total turnover)
#   cuts: the number of cuts to calculate the distance
#   maxPatch: the area of the largest patch to work as a denominator in the R-value. If null, the total area is used
#   sar.z:    the exponent of the species-area relationship
#   ptype:    is p decreasing, increasing (not supported) or constant?
#   echo:     report the advances of the function?
# Returns a table with the R-values
profile.distances <- function(db.dist, p=0:10/10, cuts=50, maxPatch=NULL, sar.z = 0.25, ptype="dec", echo=T, fn=NULL) {
  # obtaining the tables
  if (class(db.dist) != 'list') {
    cat("db.dist has to be a list\n")
    return(NULL)
  }
  if (length(db.dist) != 2) {
    cat("db.dist requires two parameters\n")
    return(NULL)
  }
  
  # any file to read?
  flv <- FALSE; lpatches <- list()
  if (!is.null(fn)) {
    if (file.exists(fn)) {
      lv <- load(fn)
      if (!identical(lv, "lpatches")) stop("Wrong file:",fn)
      flv <- TRUE
      rm(lv)
    }
  }
  
  # build the graph
  db.edges <- as.data.frame(db.dist$edges)
  i <- names(db.edges) == "dist"
  if (any(i)) names(db.edges)[i] <- "weight"
  db.vertices <- as.data.frame(db.dist$vertices)
  gr <- graph_from_data_frame(db.edges[,c("from","to","weight")], F, 
                              vertices=db.vertices[,c("fragment","area","biome")])  
  # and remove any redundant node
  gr1 <- minimum.spanning.tree(gr)
  if ((vcount(gr1) -1) != ecount(gr1)) warning("Not an minimum spanning tree.\n")
  gr1.edges <- as_data_frame(gr1, "edges")
  
  # using the predifined cuts for distance or
  # using a logarithmic scale to define the cuts
  if (length(cuts) > 1) {
    lcuts <- cuts
    lcuts <- c(lcuts[lcuts < max(gr1.edges$weight)], max(gr1.edges$weight))
    lcuts <- lcuts[order(-lcuts)]
    cuts <- length(lcuts)-1L
  } else {
    rng <- log(range(gr1.edges$weight))
    lcuts <- exp(seq(rng[2],rng[1], length.out = cuts+1))
  }
  # and preparing the output table
  pCols <- paste0("p",p)
  res <- data.frame(cut = numeric(cuts+2), numPatches = integer(cuts+2), largePatch = integer(cuts+2), 
                    matrix(0,nrow=cuts+2,ncol=length(p), dimnames=list(NULL,pCols)))
  # the algorithm will begin estimating the maximum dispersal, and will remove the links that are larger
  # from each threshold
  # full dispersal
  maxP <- (sum(db.vertices$area))
  maxS <- maxP^sar.z
  res[1,2] <- 1L;res[1,c(1,3)] <- c(lcuts[1], maxP); res[1,pCols] <- maxS
  
  # changing the maximum dispersal
  db.edges1 <- db.edges
  pm <- 0
  for (i in 2:(cuts+1)) {
    # identify the cut point
    cut <- lcuts[i]
    # is the cut available in the file?
    if (flv & (paste0("c", cut) %in% names(lpatches))) {
      lP <- lpatches[[paste0("c", cut)]]
    } else {
      # no,  retain only the edges shorter to the cut
      db.edges1 <- db.edges1[db.edges1$weight <= cut,]
      if (floor((i-1)/cuts*20) != pm) {
        if (echo) {cat("+");pm <- pm+1}  
      }
      if (echo) cat(c("-","\\","|","/")[i%%4+1])
      # rebuild the tree
      gri <- graph_from_data_frame(db.edges1[,c("from","to","weight")], F, vertices=db.vertices[,c("fragment","area")])
      if (echo) cat("\b")
      # get the subsets of the graph. Each subset will be a patch
      cyc <- components(gri)
      lP <- (aggregate(db.vertices$area, cyc["membership"], FUN=sum)[,2])
      lP <- lP[order(-lP)]
      if (!is.null(fn)) lpatches[[paste0("c", cut)]] <- lP
    }
    # and calculate the results for the patches identified
    lS <- lP^sar.z
    res[i,2] <- length(lS)
    res[i,c(1,3)] <- c(cut, max(lP)); 
    #browser()
    res[i,pCols] <- switch(ptype, "dec"=apply(t(outer(p,(1:length(lS)-1),"^"))*lS,2,sum), 
                                  "const"=apply(t(cbind(1,matrix(p,nrow=length(p),ncol=length(lS)-1)))*lS,2,sum), 
                                  "inc"=stop("Not supported yet."))
  }
  # get the resultant values for null dispersal
  lP <- (db.vertices$area)
  lP <- lP[order(-lP)]
  lS <- lP^sar.z
  res[cuts+2,2] <- length(lS)
  res[cuts+2,c(1,3)] <- c(0, max(lP))
  res[cuts+2,pCols] <- switch(ptype, "dec"=apply(t(outer(p,(1:length(lS)-1),"^"))*lS,2,sum), 
                              "const"=apply(t(cbind(1,matrix(p,nrow=length(p),ncol=length(lS)-1)))*lS,2,sum), 
                              "inc"=stop("Not supported yet."))
  # look for pre-defined reference values
  if (!is.null(maxPatch)) {
    maxP <- maxPatch
    maxS <- maxP^sar.z
  }
  # adjusting to relative values
  res[,pCols] <- res[,pCols]/maxS
  res[,3] <- res[,3]/maxP
  res <- cbind(biome=db.vertices[1,c("biome")], res)
  # saving the patches
  if (!is.null(fn)) save(lpatches, file=fn)
  if (echo & !flv) cat("\n")
  return(res)
}

# Function to create a connectivity profile plot
# profile is a data.frame
# biomes is the list of biomes to plot
# id.vars is the the name of columns in profile with the Type (potential vs. remnant), biome and cut (the maximum dispersal capability)
#         cut is assumed to be in meters
# variable.name is the name that will be used in the legend
# measure.vars are name of the columns in profile that have the actual R-values
# Returns a ggplot2 object
conn.profile <- function(profile, biomes=NULL, id.vars=c("Type", "biome","cut"), 
                          variable.name="variable",
                          measure.vars=grep("^p",names(profile),value=T), plot=T, fn=NA) {
  # transform into long format and assign the names to the columns
  if (length(measure.vars)>1) prof <- melt(profile, id.vars=id.vars, measure.vars) else prof <- copy(profile)
  setDT(prof)
  setnames(prof, id.vars, c("Type","Biomes","cut"))
  if (length(measure.vars)==1) setnames(prof, measure.vars, "value")
  if (is.null(biomes)) biomes <- unique(prof$Biomes)
  # Prepare the variables: make que variable column categorical
  prof[,variable:=as.character(variable)]
  #        and transform cut to km.
  prof[, cut:=cut/1000]
  # create the plot
  q <- ggplot(data=prof[Biomes %in% biomes & grepl("^p",variable)], aes(x=cut))+
              geom_line(aes(y=value, color=sub("p","",variable))) + facet_grid(Type~Biomes) + 
              scale_y_log10(breaks=c(0.1,0.5,1,5,10,50,100,500), minor_breaks=NULL) + scale_x_log10() + 
              labs(x = "Maximum individual dispersal distance (km)", y = "R-value", color=variable.name)
  # and plot and save if needed, and return
  if (plot) print(q)
  if (!is.na(fn)) ggsave(fn, plot=q, width = 9, height=3)
  invisible(q)
}

# Function to create the patch-size profile
# One line will be in one color and any other line in a different color. each line represents a group
# areas       is a data.frame
# biome.var   is the name of the column with information about the biome
# group.var   is the name of the column used to group the values when creating the profiles
# area.var    is the name of the column with the area of each patch
# main.group  is the name of the group in the column group.var that will be high.lighted
# others.label is the name used in the label for any gropu different to the main.group
# column.var  used when facetting is needed. It will split groups using column.var as a criteria
# col.back    the color for the groups different to main.group
# col.front   the color for the main.group
# dots        if present, the list of dots to be added to the profiles that represent the split between different
#             patches of different sizes.
# Returns a ggplot2 object
ps.profile <- function(areas,biome.var="biome", group.var="type", area.var="area", main.group="present", others.label="future",
                       column.var=NA, col.back="gray", col.front="black", dots=NULL) {
  # a function to add the breaks in the x-axis
  bref <- function(x) {
    xx <- 1:ceiling(log10(max(x)))
    xx <- c(1,10^xx,10^xx*0.3)
    xx <- xx[order(xx)]
    xx <- xx[xx < max(x)]
    return(xx)
  }
  # working with the table
  areas0 <- copy(areas)
  setDT(areas0)
  # reassigning names
  setnames(areas0, c(biome.var,group.var,area.var), c("vbiome","vgroup","varea"))
  if (!is.na(column.var)) {
    setnames(areas0, column.var, "vcol")
  } else {
    areas0[,vcol:=vgroup]
  }
  # setting the key and the order
  setkey(areas0, vbiome,vgroup,vcol,varea)
  setorder(areas0, vbiome,vgroup,vcol,-varea)
  # assigning the new values: cumulative area, and grouping color
  areas0[,":="(vcs=cumsum(varea),vpos=1:(nrow(.SD)-0),mgroup=vgroup==main.group),by=.(vbiome,vgroup,vcol)]
  areas0[,mgroupf := factor(mgroup, levels=c(T,F), labels=c(main.group, others.label))]
  # now building the basic plot
  q1 <- ggplot(data=areas0, aes(x=vpos, y=vcs, group=vgroup, col=mgroupf)) + geom_line()
  # and the facets (if needed)
  if (length(unique(areas0$vbiome)) > 1) {
    if (is.na(column.var)) {
      q1 <- q1 + facet_wrap(~vbiome, scales="free") 
    } else {
      q1 <- q1 + facet_grid(vbiome ~ vcol, scales="free_y")
    }
  } else {
    if (!is.na(column.var)) {
      q1 <- q1 + facet_grid(~ vcol)
    }
  }
  # adjusting scales and format
  q1 <- q1+scale_x_log10(breaks=bref)+ylim(c(0,NA))+theme_classic() + 
        annotate("segment", y=-Inf,yend=Inf,x=0,xend=0,size=1) + 
        scale_color_manual(values=c(col.front, col.back),labels=c(main.group, others.label), name=group.var) +
        theme(strip.background = element_rect(colour = "white", fill = "white"), strip.text=element_text(size = rel(1), face = "bold")) +
        labs(x="Patch order (log)", y=~Cumulative*phantom(0)*area*phantom(0)*(Km^2))
  suppressWarnings(q1 <- q1 + annotate("segment", y=-Inf,yend=-Inf,x=0,xend=Inf,size=1))
  
  # adding the dots
  if (!is.null(dots)) {
    # order the dots to use findInterval properly
    dots <- unique(dots[order(dots)])
    areas0[,vdots:=findInterval(varea,dots),by=.(vbiome,vgroup,vcol)]
    # now, get the first record in each group, and discard the group with the largest patch
    areas0d <- areas0[,.SD[order(-varea)][1][,.(vcs=vcs,vpos,varea)],by=.(vbiome,vgroup,vcol,vdots)][vdots != length(dots)]
    # and create a factor for the legend
    areas0d$vdotsf <- factor(areas0d$vdots,levels=0:(length(dots)-1),labels=dots)
    # and adding the dots
    q1 <- q1 + geom_point(data=areas0d[vgroup != main.group],col=col.back,mapping=aes(shape=vdotsf)) + 
                    geom_point(data=areas0d[vgroup == main.group],col=col.front,mapping=aes(shape=vdotsf)) + 
                    scale_shape(name=~Marks*phantom(1)*(Km^2), breaks=dots) 
  }
  # return the plot
  return(q1)
}

