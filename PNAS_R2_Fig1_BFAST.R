### Tested with R 4.2.0, BFAST V2.2

#library(R.matlab)
library(rhdf5)
library(bfast)
library(zoo)
library(stats)

mat.name=paste("./Americas_continent_tropical_belt_Fig1_forBFAST_Vancutsemratio95",".mat",sep="") #toString(i),
#mat.name=paste("./Africa_continent_tropical_belt_Fig1_forBFAST_Vancutsemratio95",".mat",sep="") #toString(i),
#mat.name=paste("./Asia_continent_tropical_belt_Fig1_forBFAST_Vancutsemratio95",".mat",sep="") #toString(i),


data2m <- h5read(mat.name, "/Final_time_series_FullRadar_allpixel_average_anomaly")
data2m=as.vector(data2m)

ts2 <- ts(data2m, frequency = 12, start=c(1992,10))
fit2=bfast(ts2,h=0.15,season="harmonic", breaks = 5, max.iter=10, decomp="stlplus") ### dummy  harmonic


tiff(file=paste("./Am_PNAS_R2_BFAST.tif",sep=""),
     width=4, height=6, units="in", res=300)

plot(fit2,type = c("components", "all", "data", "seasonal",
                   "trend", "noise"))

dev.off()
