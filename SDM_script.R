### In this script you will use a subset of the occurrence records to run 
### the ENM Evaluation methods presented

## new ENMEval models 
library("rgdal")
library(ape)
library(ecospat)
library(ENMeval)
library(raster)
library(RColorBrewer)
library(maptools)
library(tidyverse)
library(rasterVis)
data("wrld_simpl")

## replace this pathway to the folder holding the downloaded data
input_dir = "~/Desktop/SDMsubset_data/"
occs_dir = "~/Desktop/SDMsubset_data/ebnswi.csv"
enviro_dir = "~/Desktop/SDMsubset_data/enviro_data/"
output_dir = "~/Desktop/SDMsubset_data/results/"

## load occurrence data --------------------------------------------------
setwd(input_dir)
occs_data <- read.csv(file.path(input_dir,"ebnswi_presence.csv"))
head(occs_data,20)
occs_data <- occs_data[2:5] #remove unecessary first row
head(occs_data,20)

#subset for presence
presence <- subset(occs_data, species_observed==TRUE)
dim(presence)
coords<-presence[,2:3]
head(coords) #presence coordinates!

# load environmental layers -----------------------------------------------
setwd(enviro_dir)
files<-list.files(file.path(enviro_dir), pattern=".asc")
files #30 variables
layers<-stack(files)
crs(layers)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
layers #check resolution and crs

## Random Sampling background points--------------------------------------------
## make sure seed is set for repeatability of background point generation
set.seed(13)
bg.buff <- dismo::randomPoints(layers[[11]], p = coords, n = 30000) %>% as.data.frame()
colnames(bg.buff) <- colnames(coords)
dim(bg.buff)
head(bg.buff)
#30,000 random sample points
#make sure there is good coverage across study region!

# check that every cell is covered by background points
plot(layers[[11]], main = names(layers[[11]]))
points(bg.buff, pch = 20, cex = 0.2)
points(coords, pch=16, cex=0.2, col="red")


# create subsets of layers for different hypotheses ---------------------------
## these are layers identified by select07 method as uncorrelated
## list of layers is saved in that script so input can be read in here

## SELECT 07
pred_sel<-read.table(paste0("~/Desktop/swiftlet_data/enviro_data/ebnger_uncorrelated_var.txt"), header=FALSE)
pred_sel
layers

pred_sel<-layers[[c("CHELSA_bio4_1981","CHELSA_bio2_1981","CHELSA_bio5_1981", "glc_shv10_03",
                    "glc_shv10_05","glc_shv10_04","CHELSA_bio8_1981","CHELSA_bio9_1981","glc_shv10_01","glc_shv10_02",
                    "CHELSA_bio18_1981","glc_shv10_07","glc_shv10_09","glc_shv10_06","CHELSA_bio16_1981","glc_shv10_08")]]

## subset the rasterstack of layers to combinations for each hypothesis
#all layers
all_predictors<-pred_sel
#human layers only
human_predictors<- layers[[c("glc_shv10_04","glc_shv10_02","glc_shv10_03",
                             "glc_shv10_09","glc_shv10_06","glc_shv10_07",
                             "glc_shv10_05","glc_shv10_01","glc_shv10_08")]]
#climate layers only
climate_predictors<-layers[[c("CHELSA_bio16_1981","CHELSA_bio7_1981",
                              "CHELSA_bio11_1981","CHELSA_bio15_1981",
                              "CHELSA_bio10_1981")]]

#check
names(all_predictors)
names(climate_predictors)
names(human_predictors)
## ur a queen continue 

#plot
plot(all_predictors[[9]], main=names(all_predictors[[9]])) #checking layers
points(bg.buff, pch = 20, cex = 0.2) #checking background points
points(coords, pch=16, cex=0.2, col="red") #checking coordinates on top
#always be checking

## NA values (more checking)
plot(all_predictors[[11]], main=names(all_predictors[[11]]), colNA="red")
plot(wrld_simpl, add=TRUE)

### running models -------------------------------------------------------
# fc are feature classes and
# rm is the 
tune.args=list(fc = c("L", "LQ","H", "LQP", "LQHP", "LQHPT"), rm = c(seq(0.5,4,0.5),5:6))

##these take a while to run, make sure to adjust the number of cores you are 
##running this function in parallel. 
numCores() # do 1 or less than this number. (will slow down your computer)

#SEL7
Sel7_all <- ENMevaluate(coords, all_predictors, bg.buff, partitions="block", tune.args=tune.args, algorithm='maxent.jar', parallel=TRUE, numCores = 7)
Sel7_climate <- ENMevaluate(coords, climate_predictors, bg.buff, partitions="block", tune.args=tune.args, algorithm='maxent.jar', parallel=TRUE, numCores = 7)
Sel7_hyde <- ENMevaluate(coords, human_predictors, bg.buff, partitions="block", tune.args=tune.args, algorithm='maxent.jar', parallel=TRUE, numCores = 7)


## save models with closer cropped extent
output_dir<-"~/Desktop/SDMsubset_data/results"
save(Sel7_all, file=paste0(output_dir, "Sel7_all_ebnger.RData"))
save(Sel7_climate, file=paste0(output_dir, "Sel7_climate_ebnger.RData"))
save(Sel7_hyde, file=paste0(output_dir, "Sel7_hyde_ebnger.RData"))

### if reading in from here make sure load your .RData --------------------
load(output_dir, "Sel7_all_ebnger.RData")
#Sel7
model<-Sel7_all #all layers Sel7
#model<-Sel7_climate #climate only Sel7
#model<-Sel7_hyde #hyde only Sel7

# plotting
cols <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
countries <- maps::map("world", plot=FALSE) 
countries <- map2SpatialLines(countries, proj4string = CRS("+proj=longlat"))

##### model selection 
res <- eval.results(model)
model
res
############# filtering for max boyce
opt.boyce<- res%>% filter(cbi.val.avg==max(cbi.val.avg))
ENMeval_all_opt.boyce<-opt.boyce
### variable importance plots
mod.seq <- eval.models(model)[[opt.boyce$tune.args]]
mod.seq@lambdas
plot(mod.seq, type = "cloglog") #save as png 
par(mfrow=c(1,1)) #adjust viewing pane
pred.boyce <- eval.predictions(model)[[opt.boyce$tune.args]]
plot(pred.boyce)
plot(wrld_simpl, add=T) #check it lines up

##dimensions to save png: 880, 535
levelplot(pred.boyce, margin=FALSE, main="Edible-nest Swiftlet Select07 All Predictors", col.regions=cols, at=seq(0, 1, len=20), 
          xlab = "Longitude", ylab = "Latitude") 

##all done!! congrats b

