# Patch size profiles for the andean biomes
source("magra.R")
dirExp <- "data2/temp/"

# read the information about the biomes
tbordinal <- fread(paste0(dirExp, "ordinal-eng.txt"))
# fix some biome's names
tbordinal[Descripcion=="Evergreen montane forest", Descripcion:="Evergreen forest"]
tbordinal[Descripcion=="Seasonally dry tropical montane forest", Descripcion:="Seasonal forest"]


# read the graph files
lfiles <- dir(dirExp, pattern="^toGraph\\..*\\.RData$")
db.areas0 <- lapply(paste0(dirExp, lfiles), read.graph)
names(db.areas0) <- gsub("^toGraph\\.(.*)\\.RData$","\\1", lfiles)

# get the areas for each fragment in each biome
db.areas1 <- lapply(names(db.areas0), function(x) cbind(map=x, db.areas0[[x]]$vertices))
# and combine them in a single table
db.areas <- rbindlist(db.areas1)
# clean the memory
rm(db.areas0, db.areas1)

# recover the information stored in the file name
db.areas[,c("Model","Scenario","Period","Type","Biome"):=tstrsplit(map, split="\\.")]
db.areas[Model=="present",":="(Scenario="Present",Period=NA,Model="Present",Type=Scenario,Biome=Period)]
db.areas[,Biome:=as.integer(Biome)]

# Building present time profiles
# get the information from the table
tbo1 <- tbordinal[order(factor(Abreviacion, c("eMF","Sh","xP")))][Abreviacion %in% c("eMF","Sh","xP"),.(biome=Codigo,Abreviacion,Desc=paste(c("a)","b)","c)"), Descripcion))]
# create the figure
q1 <- ps.profile(db.areas[Model=="Present"][tbo1, on="biome"][,.(Desc,Type=factor(Type,c("pot","rem"), c("Potential","Remnant")),area=area/1e6)],
                 biome.var="Desc", group.var="Type", area.var="area", main.group="Potential", others.label="Remnant",
                 column.var=NA, col.back="gray", col.front="black", dots=c("1"=1.5,"10"=10.5))
# update the x/y labels
q1 <- q1 + labs(y=~Cumulative*phantom(0)*area*phantom(0)*(1000*phantom(0)*Km^2)) + 
  scale_y_continuous(labels=function(x) x/1000,limits=c(0,NA))
# save the file!
ggsave("ps.emf.sh.xp.pdf", q1, width=8, height=3.5) ## Fig.5

# again for the biomes in the appendix
tbo2 <- tbordinal[order(factor(Abreviacion, c("Par","hP","GCr", "dMF","PrP")))][(Abreviacion %in% c("hP","Par","GCr", "dMF","PrP")),.(biome=Codigo,Abreviacion,Desc=paste(c("a)","b)","c)","d)","e)"), Descripcion))]
# one name is too long
tbo2[,Desc := sub(" and "," and\n",Desc)]
# create the figure
q2 <- ps.profile(db.areas[Model=="Present"][tbo2, on="biome"][,.(Desc,Type=factor(Type,c("pot","rem"), c("Potential","Remnant")),area=area/1e6)],
                 biome.var="Desc", group.var="Type", area.var="area", main.group="Potential", others.label="Remnant",
                 column.var=NA, col.back="gray", col.front="black", dots=c("1"=1.5,"10"=10.5))
# update x/y labels
q2 <- q2 + labs(y=~Cumulative*phantom(0)*area*phantom(0)*(1000*phantom(0)*Km^2)) + scale_y_continuous(labels=function(x) x/1000,limits=c(0,NA))
# save it!
ggsave("ps.hp.par.gcr.dmf.prp.pdf", q2, width=8, height=4)  ## Fig.S5
# restore the names
tbo2[,Desc:=gsub("\n"," ",Desc)]

# now the future profiles
# split the table in two
tbo1 <- tbordinal[order(factor(Abreviacion, c("eMF","dMF","PrP")))][Abreviacion %in% c("eMF","dMF","PrP"),.(biome=Codigo,Abreviacion,Desc=paste(c("a)","b)","c)"), Descripcion))]
tbo2 <- tbordinal[order(factor(Abreviacion, c("Par","GCr","Sh","xP","hP")))][(Abreviacion %in% c("Par","GCr","Sh","xP","hP")),.(biome=Codigo,Abreviacion,Desc=paste(c("a)","b)","c)","d)","e)"), Descripcion))]
# and combine them
ltbo <- list(tbo1,tbo2)
# this is to save the plots
lplots <- list()
# iterate in each list
for (jtb in ltbo) {
  for (i in 1:nrow(jtb)) {
    # create the figure
    q3 <- ps.profile(db.areas[Scenario %in% c("Present","A1B")][jtb[i], on="biome"]
                     [,.(Desc,Period=Model,Type=factor(Type,c("pot","rem"), c("Potential","Remnant")),area=area/1e6)],
                     biome.var="Desc", group.var="Period", area.var="area", main.group="Present", others.label="2040-2069",
                     column.var="Type", col.back="gray", col.front="black", dots=c("1"=1.5,"10"=10.5))
    # update the labels
    q3 <- q3 + labs(y=~Cumulative*phantom(0)*area*phantom(0)*(1000*phantom(0)*Km^2)) + 
      scale_y_continuous(labels=function(x) x/1000,limits=c(0,NA)) + labs(title=jtb[i,Desc],colour=~Period) + 
      theme(title=element_text(size = rel(1), face = "bold"))
    # store the plot in the list
    lplots[[1+length(lplots)]] <- q3
  }
}

library(gridExtra)
# organize the plots and their legends
qA <- grid.arrange(lplots[[1]]+theme(legend.position = "none"),lplots[[2]]+theme(legend.position = "none"),
                   lplots[[3]]+theme(legend.box = "horizontal",legend.position = c(0.3,0.29),
                                     legend.background = element_rect(fill = "gray90", colour = NA),
                                     legend.box.background = element_rect(fill = "gray90", colour = NA)),ncol=1)
ggsave("ps.future.A.pdf",qA,width=8, height=9.6)  ## Fig. 7

qB <- grid.arrange(lplots[[4]]+theme(legend.position = "none"),lplots[[5]]+theme(legend.position = "none"),
                   lplots[[6]]+theme(legend.box = "horizontal",legend.position = c(0.35,0.2),
                                     legend.background = element_rect(fill = "gray90", colour = NA),
                                     legend.box.background = element_rect(fill = "gray90", colour = NA)),ncol=1)
ggsave("ps.future.B.pdf",qB,width=10, height=12)  ## Fig. S6a

qC <- grid.arrange(lplots[[7]]+theme(legend.position = "none"),
                   lplots[[8]]+theme(legend.box = "horizontal",legend.position = c(0.35,0.2),
                                     legend.background = element_rect(fill = "gray90", colour = NA),
                                     legend.box.background = element_rect(fill = "gray90", colour = NA)),ncol=1)
ggsave("ps.future.C.pdf",qC,width=10, height=8)   ## Fig. S6b
