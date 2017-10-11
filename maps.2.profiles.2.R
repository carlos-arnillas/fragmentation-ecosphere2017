# From graphs to values for connectivity profiles
library(igraph)
library(data.table, verbose = F)
source("magra.R")

if (!exists("sar.z")) sar.z <- 0.25
if (!exists("ptype")) ptype <- "dec"
if (!exists("dirExp")) dirExp <- "data/"
plot.profiles <- FALSE

# and the cuts to aggregate the fragments into patches
lcuts <- unique(c(seq(2000,50000,by=500), seq(50000,100000,by=1000), seq(100000,2000000,by=10000)))

# read the information about the biomes
tbordinal <- fread(paste0(dirExp,"ordinal-eng.txt"))
# fix some biome's names
tbordinal[Descripcion=="Evergreen montane forest", Descripcion:="Evergreen forest"]
tbordinal[Descripcion=="Seasonally dry tropical montane forest", Descripcion:="Seasonal forest"]

# list of biomes
biomes <- unique(tbordinal$Codigo)

# Find the graphs
lfiles <- dir(dirExp, pattern="^toGraph\\..*(\\.4069|present\\.)")
lfiles <- unique(gsub("toGraph\\.|\\.(pot|rem)\\.[[:digit:]]+.RData","",lfiles))

# Now, with the list of cases to work with, process them one by one.
# here the loop begin
for (fbiome in lfiles) {
  cat("\nWorking with:",fbiome,"\n")
  # name to use when saving the files
  fprefix <- fbiome
  # for potential and remnant maps
  for (type in c("pot","rem")) {
    # clean the temporal maps without deleting the variables
    # if all files have been done before, read the graphs
    lfiles2 <- dir(dirExp, pattern=paste0("toGraph.",fprefix,".",type))
    if (all(biomes %in% as.integer(sub("^.*\\.([[:digit:]]{1,2})\\.RData$","\\1",lfiles2)))) {
      cat(" Reading graphs.\n")
      db.dist <- lapply(paste0(dirExp, lfiles2), read.graph)
    } else {
      warning(paste0(" Graph data missing:",fprefix, " ", type))
    }
    # if distances database exists, estimate the profiles
    if (exists("db.dist")) {
      # filter the too small fragments
      db.dist2 <- lapply(db.dist, function(x) 
        criteria(x$edges, x$vertices, v.c=expression(area >= 1e+6)))
      # apply the profile.distances with a given z value
      db.profile0 <- lapply(db.dist2, 
                            function(x) profile.distances(x, cuts=lcuts, maxPatch=1, 
                                                          sar.z = sar.z, ptype=ptype,
                                                          fn = sprintf("%s/patches.%s.%s.%s.RData",
                                                                       dirExp, fprefix, x$edges$biome[1], type)))
      # merge the tables and re-organize the table properly
      db.profiles0 <- rbindlist(db.profile0)
      db.profiles0 <- cbind(file=fprefix,type=type,z=sar.z,db.profiles0)
      # plot it if needed
      if (plot.profiles) conn.profile(db.profiles0, biomes, id.vars=c("type", "biome","cut"), plot=F, fn=paste0(dirExp, "profiles.", fprefix, ".", ptype, ".", type, ".", sar.z,".pdf"))
      # and save the data
      try(save(db.profiles0, file=paste0(dirExp, "profiles.", fprefix, ".", ptype , ".", type, ".", sar.z, ".RData")), T)
      db.profiles <- if (exists("db.profiles")) rbind(db.profiles, db.profiles0) else db.profiles0
      # clean the memory
      rm(db.dist,db.dist2, db.profile0, db.profiles0)
    }
  }
}

# from wide to long format
db.profiles0 <- copy(db.profiles)

db.profiles <- melt(db.profiles0, id.vars=c("file", "type", "biome","cut", 
                                            "numPatches", "largePatch", "z"), 
                    variable.factor=F, value.name = "S")
# variable column have now the turnover rate 
# is everything ok? only one value per biome, file, type, variable, cut and z value
if (nrow(db.profiles[,.(n=nrow(.SD)),by=.(biome,file,type,variable,cut,z)][n>1])) 
  warning("something went wrong!")

# get the richness of the unitary tolerant community at the present.
db.profiles[,SuT:=mean(NA^(!(file=="present" & type=="pot" & variable=="p0" & numPatches==1))*S,na.rm=TRUE),by=biome]
# and get the R values
db.profiles[,R:=S/SuT]
# organize columns
db.profiles[,c("Model","Scenario","Year"):=tstrsplit(file,".",fixed=T)]

# quick print of the R values for the extreme communities
# db.profiles[file=="present" & ((numPatches==1 & variable=="p0") | (cut == 0 & variable %in% c("p0","p1"))),R,by=.(biome,type,cut,z)]

# and fix some names
db.profiles <- tbordinal[,.(Codigo,Descripcion,Abreviacion)][db.profiles, on=c(Codigo="biome")]
setnames(db.profiles,"Abreviacion","biome")

# save the file
save(db.profiles, file=paste0(dirExp,"db.profiles.4069.",sar.z, ".", ptype,".RData"))

