# Building the patch-size profile
#rm(list=ls())
source("clumps.ps-prof.R")

# Building the profiles
rm(list=ls())
lsar.z <- c(0.1,0.25,1, 0.5, 0.4,0.6,0.9)
lptype <- c("const","dec")
for (ptype in lptype) {
  for (sar.z in lsar.z) {
    source("maps.2.profiles.2.R")
    source("clumps.conn-prof.uncertainty.R")
    rm(list=setdiff(ls(),c("ptype", "sar.z", "lptype","lsar.z")))
  }
}
rm(ptype,sar.z,lptype,lsar.z)


# pulling together uncertainties
library(data.table)
dirExp <- "data/"
lptype <- c("const","dec")
for (ptype in lptype) {
  lfilesU <- dir(dirExp, pattern=sprintf("uncertainty\\..*\\.%s\\.csv",ptype))
  uncertainties0 <- rbindlist(lapply(lfilesU, function(x) cbind(z0=gsub(sprintf("uncertainty\\.|\\.%s\\.csv",ptype), "", x),
                                                               fread(paste0(dirExp, x)))))
  uncertainties0[,V1:=NULL]
  uncertainties0[,z := as.numeric(z0)]
  uncertainties0[,Component:=factor(component, levels=c("Scenario", "Model", "type", 
                                                        "as.factor(cut)", "variable", 
                                                        "as.factor(cut):variable", "Residuals"),
                                               labels=c("Climate scenario", "Climate model",
                                                        "Tolerance", "Dispersal", "Turnover",
                                                        "Dispersal:Turnover", "Residuals"))]
  
  uncertainties <- melt(uncertainties0, id.vars=c("z","z0","value.var","component","Component"), variable.factor = FALSE,
                        variable.name = "Biome")
  
  library(ggplot2)
  ggplot(uncertainties[value.var=="pssq"], aes(x=z, y=value, colour=Component)) + geom_line() + facet_wrap(~Biome) +
         geom_vline(xintercept = 0.25, colour="gray30", linetype="dashed") +
         geom_vline(xintercept = 0.4, colour="gray30", linetype="dotted") + labs(y="Proportion explained", x="Exponent (z) of the SAR model")
  ggsave(sprintf("proportion explained_z-range_%s.pdf",ptype), width=9, height=5) # Fig. S7a and b
}
