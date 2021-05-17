#############################################################
## PPM For Foraging Model using Sample GPS data #############
#############################################################
#This code uses a subset of real GPS data used for analysis as an example of how to run code, and so results differ

library(INLA, quietly=TRUE)
library(sp, quietly=TRUE)
library(rgdal, quietly=TRUE)
library(rgeos, quietly=TRUE)
library(raster, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(moveHMM, quietly=TRUE)
library(spatialEco, quietly=TRUE)
library(spdep, quietly=TRUE)

source("Functions/MakeSpatialRegion.R")
source("Functions/MakeIntegrationStack.R")
source("Functions/GetNearestCovariate.R")
source("Functions/MakeBinomStack.R")
source("Functions/FitModel.R")
source("Functions/MakeIntegrationStackPts.R")
source("Functions/MakeProjectionGrid.R")
source("Functions/MakeProjectionGridPts.R")
source("Functions/MarkedMakePointsStack.R")
source("Functions/MarkedMakePointsStackPtsOnly.R")
source("Functions/FitModelPtsOnly.R")
##########################################################################

load("R:/rsrch/cb751/phd/nmp513/Chapter1/Github/HMM_sub.all.Rdata")                    ## reads in HMM model (from HMM script)
load("R:/rsrch/cb751/phd/nmp513/Chapter1/Github/GPSdata_withcovariates.Rdata")     ## loads the GPS point locations cleaned and with covariates

###########################
######    Prep data    ####
###########################

### Compute the most likely activities at each time point:
GPSdata_withcovariates$states <- viterbi(HMM_thermals_time_altitude)

# Remove NAs
vult.locs <- GPSdata_withcovariates[!is.na(GPSdata_withcovariates$Long) & !is.na(GPSdata_withcovariates$Lat),]


# Add Season 
getSeason <- function(DATES) {
  Wet <- as.Date("2012-12-01", format = "%Y-%m-%d") 
  Dry <- as.Date("2012-5-31",  format = "%Y-%m-%d") 
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2012-%m-%d"))
  ifelse (d >= Wet | d < Dry, "Wet", "Dry")
}

vult.locs$season <- (getSeason(vult.locs$Date) == "Wet")* 1  # Wet = 1 Dry =0

# Create data for season/id plots
GPSdata_withcovariates$season <- getSeason(GPSdata_withcovariates$Date)
data_summary <- aggregate(rep(1, NROW(GPSdata_withcovariates)), by = list(state = GPSdata_withcovariates$states, ID = GPSdata_withcovariates$ID,
                                                                          season = GPSdata_withcovariates$season), FUN = sum)

##turn data into an sp object
NDVI  <- raster("R:/rsrch/cb751/phd/nmp513/Chapter1/Data/Covariate_TIFFS/TZFULL_NDVI.tif")
NDVI[NDVI<0] <- NA
coordinates(vult.locs) <- ~ Long + Lat
proj4string(vult.locs) <- proj4string(NDVI)

## Add a spatial interpolation / density plot of the vulture locations - all fixes for observer effect. This takes a long time so better to save and skip 
pt_dens <- sp.kde(vult.locs,  bw = 0.1, newdata = aggregate(NDVI, 10), standardize = FALSE, mask = TRUE)

vult.locs <- spTransform(vult.locs, "+init=EPSG:21037")  ## to reproject to UTM 

# Get spatial projection for boundary
proj <- crs(vult.locs)

## Import covariate data

tz.protectedarea <- raster("R:/rsrch/cb751/phd/nmp513/Chapter1/Data/Covariate_TIFFS/ProtectedAreas/TZandMOZ.tif")
tz.livestock <- raster("R:/rsrch/cb751/phd/nmp513/Chapter1/Data/Covariate_TIFFS/Livestock/AF_TLU_ruminants.tif")
tz.tree.cover   <- raster("R:/rsrch/cb751/phd/nmp513/Chapter1/Data/Covariate_TIFFS/TZFULL_TreeCover.tif")
tz.albedoB <- raster("R:/rsrch/cb751/phd/nmp513/Chapter1/Data/Covariate_TIFFS/FULLTZ_BlackSky_Albedo.tif")
tz.emissivity <- raster("R:/rsrch/cb751/phd/nmp513/Chapter1/Data/Covariate_TIFFS/TZFULL_Emissivity.tif")
tz.selous <- raster("R:/rsrch/cb751/phd/nmp513/Chapter1/Data/Covariate_TIFFS/Selousarea.tif")
rivers.dist<- raster("R:/rsrch/cb751/phd/nmp513/Chapter1/Data/Covariate_TIFFS/Rivers/Dist_to_Rivers.tif")

##Build Mesh
coord<- as.matrix(coordinates(vult.locs))
obs.window <- inla.nonconvex.hull(coord, -0.03, -0.05, resolution = c(110, 110))
obs.window.poly <- SpatialPolygons(list(Polygons(list(Polygon(obs.window$loc[1:which(obs.window$idx[,2] == 1), ])), 1)), proj4string = proj)
Meshpars <- list(max.edge = c(16000, 55000), offset = c(11000,110000), cutoff = 22000)

Mesh <- MakeSpatialRegion(
  data = NULL,
  bdry = obs.window.poly,        
  meshpars = Meshpars, 
  proj = proj
)

mesh_coords <- data.frame(Mesh$mesh$loc[,1:2])  # this section added
names(mesh_coords) <- c("x", "y")
coordinates(mesh_coords) <- ~ x + y

## create a grid of points for projection:
grid.extent <- extent(vult.locs)
grid.resolution <- 10000  ## 10km
projection.grid <- SpatialPoints(expand.grid(x = seq(grid.extent@xmin, grid.extent@xmax, by = grid.resolution),
                                             y = seq(grid.extent@ymin, grid.extent@ymax, by = grid.resolution)),
                                 proj4string = proj)

## then create a covariates dataset that is the extraction of the covariates (and recorder effort) at all these grid locations
covariates <- as(projection.grid, "SpatialPointsDataFrame")
covariates@data <- data.frame(PA = rep(0, NROW(projection.grid)))

# Extract the covariates to grid
covariates$rivers <- extract(rivers.dist, projection.grid)
covariates$PA <- extract(tz.protectedarea, projection.grid)
covariates$PA[is.na(covariates$PA)] <- 4                          # NA points outside of PA areas set to 4
covariates$PA <- as.factor(covariates$PA)       
covariates$PA2 <- (covariates$PA == 2)*1
covariates$PA3 <- (covariates$PA == 3)*1
covariates$PA4 <- (covariates$PA == 4)*1
covariates$tlu <- extract(tz.livestock, projection.grid)
covariates$tlu[covariates$tlu > 100] = 100 #sets values over 100 to NA
covariates$log.tlu <- log(covariates$tlu +1)
covariates$treec   <- extract(tz.tree.cover, projection.grid)
covariates$treec <- abs(covariates$treec)                         #imports as negative. changes to positive
covariates$emissivity <- extract(tz.emissivity, projection.grid)
covariates$black_alb <- extract(tz.albedoB, projection.grid)
covariates$obs_eff <- extract(pt_dens, projection.grid)
covariates$selous <- extract(tz.selous, projection.grid)
covariates$selous[is.na(covariates$selous)] <- 0
covariates <- covariates[!is.na(rowSums(covariates@data[,-1])),] #delete all NAs in dataset

# quadratic variables
covariates$emissivity2 <- covariates$emissivity ^ 2
covariates$log.tlu2 <- covariates$log.tlu ^ 2
covariates$black_alb2 <- covariates$black_alb ^ 2

## then scale covariates (not the effort or PA)

covariates$log.tlu.s   <- c(scale(covariates$log.tlu))
covariates$treec.s   <- c(scale(covariates$treec))   
covariates$emissivity.s <- c(scale(covariates$emissivity))
covariates$black_alb.s <- c(scale(covariates$black_alb))
covariates$emissivity2.s <- c(scale(covariates$emissivity2))
covariates$log.tlu2.s <- c(scale(covariates$log.tlu2))
covariates$black_alb2.s <- c(scale(covariates$black_alb2))
covariates$rivers.s <- c(scale(covariates$rivers))

#################################################
###         Obs matrix for foraging          ###        
#################################################
state <- 3

#skip below after running and saving
load("class.nonclass.sub.Rdata")   

#############################################################################
class.locs <- vult.locs[vult.locs$states == state,]

## and exclude any duplicated points:
class.locs <- class.locs[!duplicated(coordinates(class.locs)),]

knn1 <- knearneigh(class.locs, k=1)$nn  
class.locs$dist <- 0

nonclass.locs <- vult.locs[vult.locs$states != state,]

## remove any points where the birds have been recorded in this activity type:
nonclass.locs <- nonclass.locs[!apply(coordinates(nonclass.locs), 1, paste, collapse = "_") %in% apply(coordinates(class.locs), 1, paste, collapse = "_"),]
## and remove any duplicated points:
nonclass.locs <- nonclass.locs[!duplicated(coordinates(nonclass.locs)),]

nonclass.locs$dist <- 0
nonclass.locs <- nonclass.locs[sample(1:NROW(nonclass.locs), NROW(class.locs), replace = FALSE),]   ### takes a sample of nonclass.locs

for(i in 1:NROW(class.locs)) class.locs$dist[i] <- gDistance(class.locs[i,], class.locs[knn1[i,1],])
for(i in 1:NROW(nonclass.locs)) {nonclass.locs$dist[i] <- gDistance(nonclass.locs[i,], class.locs); if(i %% 1000 == 0) print(paste(i, "completed"))}
#

save(class.locs, nonclass.locs, file="class.nonclass.sub.Rdata")

#####################################################################
###############         Build Model               ###################


# Make stack for background mesh
stk.ip <- MakeIntegrationStackPts(
  mesh = Mesh$mesh,
  data = covariates, 
  area = Mesh$w,
  tag = "ip",
  InclCoords = TRUE,
  coordnames = c("x", "y")
)


if(!exists("Nxy.scale")) Nxy.scale <- 10000   ## make sure scale is appropraite for UTM 

Boundary <- Mesh$mesh$loc[Mesh$mesh$segm$int$idx[, 2], ]
Nxy.size <- c(diff(range(Boundary[, 1])), diff(range(Boundary[, 2])))
Nxy <- round(Nxy.size / Nxy.scale)

# Make stack for projections
stk.pred <- MakeProjectionGridPts(
  nxy = Nxy,
  mesh = Mesh$mesh,
  data = covariates,
  tag = "pred.cur",  
  boundary = Boundary, #shapefile
  InclCoords = TRUE,
  coordnames = c("x", "y")
)


##  make a stack for the actual activity points
stk.classlocs <- MakeMarkedPointsStackPts(
  presences = class.locs,
  data = covariates,
  mesh = Mesh$mesh,
  tag = "vp",
  InclCoords = TRUE,
  coordnames = c("x", "y")
)


#save stacks and mesh to use in modelling - can be loaded without having to process it all again
save(stk.ip, stk.pred, stk.classlocs, Mesh, covariates, file = "R:/rsrch/cb751/phd/nmp513/Chapter1/Github/VultureForaging.RData")

########################################################################################

##### Fit model

load("R:/rsrch/cb751/phd/nmp513/Chapter1/Github/VultureForaging.RData")

scl <- 1000  # This will scale everything to units of kilometres (not m, the UTM resolution) adjust to make best!

# create priors 
C.F. <- list(
  mean = list(int.vp = 0),
  mean.intercept = 0,
  prec = list(int.vp = 1),
  prec.intercept = 0.001
)

Mesh$mesh$loc <- Mesh$mesh$loc / scl   #scale Mesh
Mesh$w <- Mesh$w / (scl^2)

spde <- inla.spde2.pcmatern(
  mesh = Mesh$mesh,
  alpha = 2,
  prior.range = c(10000/scl, 0.5),  ## in UTM
  prior.sigma = c(5, 0.1)  
)


## define some linear combinations to use to generate effect plots: ####
lc.PArd = inla.make.lincombs(Intercept = rep(1, 4),
                             PA2 = c(0,1,0,0),
                             PA3 = c(0,0,1,0),
                             PA4 = c(0,0,0,1)
)
names(lc.PArd) <- paste0("PArd", 1:4)

lc.PAs = inla.make.lincombs(Intercept = rep(1, 4),
                            selous = rep(1,4),
                            PA2 = c(0,1,0,0),
                            PA3 = c(0,0,1,0),
                            PA4 = c(0,0,0,1),
                            'selous:PA2' = c(0,1,0,0),
                            'selous:PA3' = c(0,0,1,0),
                            'selous:PA4' = c(0,0,0,1)
)
names(lc.PAs) <- paste0("PAs", 1:4)
lc.PArw = inla.make.lincombs(Intercept = rep(1, 4),
                             season = rep(1,4),
                             PA2 = c(0,1,0,0),
                             PA3 = c(0,0,1,0),
                             PA4 = c(0,0,0,1),
                             'PA2:season' = c(0,1,0,0),
                             'PA3:season' = c(0,0,1,0),
                             'PA4:season' = c(0,0,0,1)
)
names(lc.PArw) <- paste0("PArw", 1:4)
######################################
## R=ruaha d= dry season

lc.riversRd = inla.make.lincombs(Intercept = rep(1, 100),
                                 rivers.s = seq(min(covariates$rivers.s, na.rm = T), max(covariates$rivers.s, na.rm = TRUE), length = 100)) 
names(lc.riversRd) <- paste0("riversRd", 1:100)

lc.riversS = inla.make.lincombs(Intercept = rep(1, 100),
                                selous = rep(1,100),
                                rivers.s = seq(min(covariates$rivers.s, na.rm = T), max(covariates$rivers.s, na.rm = TRUE), length = 100),
                                'selous:rivers.s' = seq(min(covariates$rivers.s, na.rm = T), max(covariates$rivers.s, na.rm = TRUE), length = 100))
names(lc.riversS) <- paste0("riversS", 1:100)
lc.riversRw = inla.make.lincombs(Intercept = rep(1, 100),
                                 season = rep(1,100),
                                 rivers.s = seq(min(covariates$rivers.s, na.rm = T), max(covariates$rivers.s, na.rm = TRUE), length = 100),
                                 'rivers.s:season' = seq(min(covariates$rivers.s, na.rm = T), max(covariates$rivers.s, na.rm = TRUE), length = 100))
names(lc.riversRw) <- paste0("riversRw", 1:100)

######################################

lc.treeRd = inla.make.lincombs(Intercept = rep(1, 100),
                               treec.s = seq(min(covariates$treec.s, na.rm = T), max(covariates$treec.s, na.rm = TRUE), length = 100)) 
names(lc.treeRd) <- paste0("treeRd", 1:100)

lc.treeS = inla.make.lincombs(Intercept = rep(1, 100),
                              selous = rep(1,100),
                              treec.s = seq(min(covariates$treec.s, na.rm = T), max(covariates$treec.s, na.rm = TRUE), length = 100),
                              'selous:treec.s' = seq(min(covariates$treec.s, na.rm = T), max(covariates$treec.s, na.rm = TRUE), length = 100))
names(lc.treeS) <- paste0("treeS", 1:100)
lc.treeRw = inla.make.lincombs(Intercept = rep(1, 100),
                               season = rep(1,100),
                               treec.s = seq(min(covariates$treec.s, na.rm = T), max(covariates$treec.s, na.rm = TRUE), length = 100),
                               'treec.s:season' = seq(min(covariates$treec.s, na.rm = T), max(covariates$treec.s, na.rm = TRUE), length = 100))
names(lc.treeRw) <- paste0("treeRw", 1:100)
#######################################

lc.tluRd = inla.make.lincombs(Intercept = rep(1, 100),
                              log.tlu.s = seq(min(covariates$log.tlu.s, na.rm = T), max(covariates$log.tlu.s, na.rm = TRUE), length = 100), 
                              log.tlu2.s = seq(min(covariates$log.tlu2.s, na.rm = T), max(covariates$log.tlu2.s, na.rm = TRUE), length = 100)) 
names(lc.tluRd) <- paste0("tluRd", 1:100)

lc.tluS = inla.make.lincombs(Intercept = rep(1, 100),
                             selous = rep(1,100),
                             log.tlu.s = seq(min(covariates$log.tlu.s, na.rm = T), max(covariates$log.tlu.s, na.rm = TRUE), length = 100),
                             log.tlu2.s = seq(min(covariates$log.tlu2.s, na.rm = T), max(covariates$log.tlu2.s, na.rm = TRUE), length = 100),
                             'selous:log.tlu.s' = seq(min(covariates$log.tlu.s, na.rm = T), max(covariates$log.tlu.s, na.rm = TRUE), length = 100),
                             'selous:log.tlu2.s' = seq(min(covariates$log.tlu2.s, na.rm = T), max(covariates$log.tlu2.s, na.rm = TRUE), length = 100))
names(lc.tluS) <- paste0("tluS", 1:100)

lc.tluRw = inla.make.lincombs(Intercept = rep(1, 100),
                              season = rep(1,100),
                              log.tlu.s = seq(min(covariates$log.tlu.s, na.rm = T), max(covariates$log.tlu.s, na.rm = TRUE), length = 100),
                              log.tlu2.s = seq(min(covariates$log.tlu2.s, na.rm = T), max(covariates$log.tlu2.s, na.rm = TRUE), length = 100),
                              'log.tlu.s:season' = seq(min(covariates$log.tlu.s, na.rm = T), max(covariates$log.tlu.s, na.rm = TRUE), length = 100),
                              'log.tlu2.s:season' = seq(min(covariates$log.tlu2.s, na.rm = T), max(covariates$log.tlu2.s, na.rm = TRUE), length = 100))
names(lc.tluRw) <- paste0("tluRw", 1:100)
#######################################

lc.emisRd = inla.make.lincombs(Intercept = rep(1, 100),
                               emissivity.s = seq(min(covariates$emissivity.s, na.rm = T), max(covariates$emissivity.s, na.rm = TRUE), length = 100),
                               emissivity2.s = seq(min(covariates$emissivity2.s, na.rm = T), max(covariates$emissivity2.s, na.rm = TRUE), length = 100))
names(lc.emisRd) <- paste0("emisRd", 1:100)

lc.emisS = inla.make.lincombs(Intercept = rep(1, 100),
                              selous = rep(1,100),
                              emissivity.s = seq(min(covariates$emissivity.s, na.rm = T), max(covariates$emissivity.s, na.rm = TRUE), length = 100),
                              emissivity2.s = seq(min(covariates$emissivity2.s, na.rm = T), max(covariates$emissivity2.s, na.rm = TRUE), length = 100),
                              'selous:emissivity.s' = seq(min(covariates$emissivity.s, na.rm = T), max(covariates$emissivity.s, na.rm = TRUE), length = 100),
                              'selous:emissivity2.s' = seq(min(covariates$emissivity2.s, na.rm = T), max(covariates$emissivity2.s, na.rm = TRUE), length = 100))
names(lc.emisS) <- paste0("emisS", 1:100)

lc.emisRw = inla.make.lincombs(Intercept = rep(1, 100),
                               season = rep(1,100),
                               emissivity.s = seq(min(covariates$emissivity.s, na.rm = T), max(covariates$emissivity.s, na.rm = TRUE), length = 100),
                               emissivity2.s = seq(min(covariates$emissivity2.s, na.rm = T), max(covariates$emissivity2.s, na.rm = TRUE), length = 100),
                               'emissivity.s:season' = seq(min(covariates$emissivity.s, na.rm = T), max(covariates$emissivity.s, na.rm = TRUE), length = 100),
                               'emissivity2.s:season' = seq(min(covariates$emissivity2.s, na.rm = T), max(covariates$emissivity2.s, na.rm = TRUE), length = 100))
names(lc.emisRw) <- paste0("emisRw", 1:100)
######################################

lc.albRd = inla.make.lincombs(Intercept = rep(1, 100),
                              black_alb.s = seq(min(covariates$black_alb.s, na.rm = T), max(covariates$black_alb.s, na.rm = TRUE), length = 100),
                              black_alb2.s = seq(min(covariates$black_alb2.s, na.rm = T), max(covariates$black_alb2.s, na.rm = TRUE), length = 100)) 
names(lc.albRd) <- paste0("albRd", 1:100)

lc.albS = inla.make.lincombs(Intercept = rep(1, 100),
                             selous = rep(1,100),
                             black_alb.s = seq(min(covariates$black_alb.s, na.rm = T), max(covariates$black_alb.s, na.rm = TRUE), length = 100),
                             black_alb2.s = seq(min(covariates$black_alb2.s, na.rm = T), max(covariates$black_alb2.s, na.rm = TRUE), length = 100),
                             'selous:black_alb.s' = seq(min(covariates$black_alb.s, na.rm = T), max(covariates$black_alb.s, na.rm = TRUE), length = 100),
                             'selous:black_alb2.s' = seq(min(covariates$black_alb2.s, na.rm = T), max(covariates$black_alb2.s, na.rm = TRUE), length = 100))
names(lc.albS) <- paste0("albS", 1:100)

lc.albRw = inla.make.lincombs(Intercept = rep(1, 100),
                              season = rep(1,100),
                              black_alb.s = seq(min(covariates$black_alb.s, na.rm = T), max(covariates$black_alb.s, na.rm = TRUE), length = 100),
                              black_alb2.s = seq(min(covariates$black_alb2.s, na.rm = T), max(covariates$black_alb2.s, na.rm = TRUE), length = 100),
                              'black_alb.s:season' = seq(min(covariates$black_alb.s, na.rm = T), max(covariates$black_alb.s, na.rm = TRUE), length = 100),
                              'black_alb.s:season' = seq(min(covariates$black_alb2.s, na.rm = T), max(covariates$black_alb2.s, na.rm = TRUE), length = 100))
names(lc.albRw) <- paste0("albRw", 1:100)

######################################

lc.all <- c(lc.PArd, lc.riversRd, lc.treeRd, lc.tluRd,lc.emisRd,lc.albRd,
            lc.emisS, lc.riversS, lc.tluS,lc.treeS, lc.PAs, lc.albS,
            lc.emisRw, lc.riversRw, lc.tluRw,lc.treeRw, lc.PArw, lc.albRw
)
#############################################################

form <- formula(resp ~ 0 + selous * (PA2 + PA3 + PA4 + rivers.s + treec.s + log.tlu.s +  log.tlu2.s + emissivity.s + emissivity2.s + black_alb.s + black_alb2.s) +
                  season * (PA2 + PA3 + PA4 + rivers.s + treec.s + log.tlu.s +  log.tlu2.s + emissivity.s + emissivity2.s + black_alb.s + black_alb2.s) +
                  Intercept 
                + f(ID, model = "iid") 
                + x + y 
                + obs_eff  
                + int.vp 
                + f(i, model = spde)
)

dat.stk <- inla.stack(stk.ip, stk.classlocs, stk.pred$stk)  # combined stacks
dat.stk$data$data$e <- dat.stk$data$data$e / (scl^2)  # scale areas


bird_model <- FitModelPts(
  dat.stk,
  formula = form,
  CovNames = NULL,       #leave
  mesh = Mesh$mesh,
  lincombs = lc.all,
  verbose = TRUE,
  predictions.cur = TRUE,   
  predictions.fut = FALSE,  #for future
  control.fixed = C.F.,     #list of priors
  waic = TRUE,              #model fit
  nthreads = 1              #cpu threads
)


save(
  C.F., spde, form, bird_model,
  file = "R:/rsrch/cb751/phd/nmp513/Chapter1/Github/VultureForagingModel.RData")
