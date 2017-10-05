# Analysis of the biomes' graphs
source("codigos/magra.R")
library(data.table, verbose = F)
if (!exists("sar.z")) sar.z <- 0.25

# read the relevant file
if (!exists("db.profiles")) {
  dirExp <- "data2/temp/"
  if (!exists("ptype")) ptype <- "dec"
    
  # read the information about the biomes
  tbordinal <- fread(paste0(dirExp, "ordinal-eng.txt"))
  # fix some biome's names
  tbordinal[Descripcion=="Evergreen montane forest", Descripcion:="Evergreen forest"]
  tbordinal[Descripcion=="Seasonally dry tropical montane forest", Descripcion:="Seasonal forest"]

  lprofiles <- dir(dirExp, pattern="^db.profiles\\.[[:digit:].]*(dec|const)\\.RData$")
  # choose file to focus on
  if (ptype=="dec") try(load (paste0(dirExp,"db.profiles.4069.",sar.z,".RData")),silent=TRUE)
  if (!exists("db.profiles")) load (paste0(dirExp,"db.profiles.4069.",sar.z,".",ptype,".RData")) 
}

sar.z <- unique(db.profiles$z)
# confirm only 1 sar.z value exists
if (length(na.omit(sar.z)) != 1) stop("No z-value, or more than 1.")

# get the connectivity profiles for three biomes
q3 <- conn.profile(db.profiles[file=="present" & biome %in% c("eMF","Sh","xP"),
                                .(type2=factor(type, levels=c("pot","rem"),labels=c("Tolerant","Intolerant")),
                                  Descripcion=factor(Descripcion,c("Evergreen forest", "Montane shrublands", "Xeric puna")),
                                  cut,variable,R)], 
                    id.vars=c("type2","Descripcion","cut"), measure.vars = "R", variable.name = "Species\nturnover")
ggsave(paste0("conn.profile.emf.sh.xp.",sar.z,".",ptype,".pdf"), plot=q3, width = 7, height=3.5)

# create the table of extreme LEM-communities
tabExtremes <- dcast(db.profiles[file=="present" & ((numPatches==1 & variable=="p0") | (cut == 0 & variable %in% c("p0","p1"))),
                                 .(Descripcion,biome,type,unitary=(cut>1),nested=(variable=="p0"),R)][
                     ,.(Descripcion,biome,type=ifelse(type=="pot", "Tolerant", "Intolerant"),
                         comm=ifelse(unitary,"Unitary",ifelse(nested,"Nested","Disjoint")),R)],
                     Descripcion+biome~comm+type)
setcolorder(tabExtremes,c("Descripcion","biome","Unitary_Tolerant", "Unitary_Intolerant", 
                          "Nested_Tolerant", "Nested_Intolerant", "Disjoint_Tolerant", "Disjoint_Intolerant"))
# export the extreme values
write.csv(tabExtremes,file=paste0(dirExp, "extremes.",sar.z,".",ptype,".csv"))


# Do not include present time conditions, unitary communities or fully disjoint ones
# and use only ten cuts to fit linear models per biome
uncert2 <- db.profiles[file!="present" & numPatches != 1 & variable!="p1",
                       { # get the cuts
                         lc <- unique(.SD[,cut])
                         # sort them
                          lc <- lc[order(lc)]
                          # and sample them homogeneously
                          lc <- lc[seq(1,length(lc),length.out=10)]
                          # report it
                          cat("Cuts used in biome", unique(biome), ":", lc)
                          # fit the model
                          x <- anova(lm(R ~ Scenario+ Model + type + as.factor(cut) * variable, data=.SD[cut %in% lc]))
                          # recover the key values
                          colnames(x) <- c("Df","ssq","mssq","F","p")
                          cat("\n")
                          # return the anova 
                          cbind(component=factor(row.names(x),levels=row.names(x)),x[,c("Df","ssq","mssq")])
                       }, by=biome]

# get the proportion and export it    
uncert2[,pssq := ssq/sum(ssq)*100,by=biome]
write.csv(rbindlist(list(
  cbind(value.var="Df",   dcast(tbordinal[uncert2,on=c(Abreviacion="biome")],component~Descripcion, value.var = "Df")),
  cbind(value.var="pssq", dcast(tbordinal[uncert2,on=c(Abreviacion="biome")],component~Descripcion, value.var = "pssq")),
  cbind(value.var="mssq", dcast(tbordinal[uncert2,on=c(Abreviacion="biome")],component~Descripcion, value.var = "mssq")))),
  file=paste0(dirExp, "uncertainty.",sar.z,".",ptype,".csv"))
