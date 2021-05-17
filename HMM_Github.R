##########################################################
## Hidden Markov Model Fora  Sample GPS data #############
##########################################################
#This code uses a subset of real GPS data used for analysis as an example of how to run code, and so results differ

library(moveHMM)

load("GPSdata_withcovariates.Rdata")     #subset of GPS data cleaned and covariates added                                

########################################################################################

vultsCOV <- prepData(GPSdata_withcovariates, type="LL", coordNames=c("Long","Lat"))  #Warning from NA's in height  


###########################################################
##                   Prep Covariates                     ##
###########################################################
#####           declare factors and scale                ##
############################################################

## Group habitats

vultsCOV$habitat <- as.factor(vultsCOV$habitat)
#table(vultsCOV$habitat)  ### shows which habitats and which are most used
vultsCOV$habitat.type <- as.character(vultsCOV$habitat)

vultsCOV$habitat.type[vultsCOV$habitat.type %in% c("4")]  <- "Forest"
vultsCOV$habitat.type[vultsCOV$habitat.type %in% c("9", "10")]  <- "Savanna"
vultsCOV$habitat.type[vultsCOV$habitat.type %in% c("12")]  <- "Other"

vultsCOV$habitat.type <- as.factor(vultsCOV$habitat.type)

hab.indicators <- model.matrix(vultsCOV$ID ~ -1 + vultsCOV$habitat.type)
vultsCOV$Forest <- hab.indicators[,"vultsCOV$habitat.typeForest"]
vultsCOV$Savanna <- hab.indicators[,"vultsCOV$habitat.typeSavanna"]
vultsCOV$Other <- hab.indicators[,"vultsCOV$habitat.typeOther"]

## Scale Covariates 

vultsCOV$easting <- sin(2*pi*vultsCOV$aspect/360)  ## angles in radians switched to easting
vultsCOV$easting.s <- c(scale(vultsCOV$easting))

vultsCOV$slope.s <- c(scale(vultsCOV$slope))
vultsCOV$height.s <- c(scale(vultsCOV$height.above.ground))

vultsCOV$emissivity.s <- c(scale(vultsCOV$emissivity))
vultsCOV$albedoW.s <- c(scale(vultsCOV$albedo_white))
vultsCOV$albedoB.s <- c(scale(vultsCOV$albedo_black))

vultsCOV$aspect.s <- c(scale(vultsCOV$aspect))
vultsCOV$elevation.s <- c(scale(vultsCOV$elevation))
vultsCOV$ndvi.s <- c(scale(vultsCOV$ndvi))
vultsCOV$treecover.s <- c(scale(vultsCOV$treecover))

############################################################
###       Change time to time since 6am (dawn):           ##
############################################################

vultsCOV$Time.sincedawn <- as.numeric(format(vultsCOV$timestamp, "%H")) -6
vultsCOV$Time.sincedawn[vultsCOV$Time.sincedawn == -6] <- 18
vultsCOV$Time.sincedawn[vultsCOV$Time.sincedawn == -5] <- 19
vultsCOV$Time.sincedawn[vultsCOV$Time.sincedawn == -4] <- 20
vultsCOV$Time.sincedawn[vultsCOV$Time.sincedawn == -3] <- 21
vultsCOV$Time.sincedawn[vultsCOV$Time.sincedawn == -2] <- 22
vultsCOV$Time.sincedawn[vultsCOV$Time.sincedawn == -1] <- 23

## save to skip later

#save(vultsCOV, file = "vultsCOV_sub.Rdata")
load("vultsCOV_sub.Rdata")

##########################################################
###                start analysis                     ####
##########################################################


#########################################################
####        4 state chosen parameters                 ###
#########################################################

mu0           <- c(0.2, 10, 25, 45)  
sigma0        <- c(0.5, 5, 5, 10) 
zeroMass      <- c(0.1, 0.1, 0.1, 0.1)
stepPar0      <- c(mu0, sigma0, zeroMass) 
kappa0        <- c(0.1, 1, 1, 0.5)   
angleMean     <- c(0, 0, 0, 0)       
anglePar0     <- c(angleMean, kappa0)

## Create different HMM model combinations 
HMMno_cov <- fitHMM(vultsCOV,
                    nbStates = 4,
                    stepPar0,
                    anglePar0,
                    formula = ~1,
                    stepDist = "gamma",
                    angleDist = "vm",
                    verbose = 1)
HMM_time <- fitHMM(vultsCOV,
                   nbStates = 4,
                   stepPar0,
                   anglePar0,
                   formula = ~Time.sincedawn,
                   stepDist = "gamma",
                   angleDist = "vm",
                   verbose = 1)
HMM_habitat <- fitHMM(vultsCOV,
                      nbStates = 4,
                      stepPar0,
                      anglePar0,
                      formula = ~Forest + Savanna + Other,           
                      stepDist = "gamma",
                      angleDist = "vm",
                      verbose = 1)
HMM_ndvi <- fitHMM(vultsCOV,
                   nbStates = 4,
                   stepPar0,
                   anglePar0,
                   formula = ~ndvi.s,  
                   stepDist = "gamma",
                   angleDist = "vm",
                   verbose = 1)
HMM_treec <- fitHMM(vultsCOV,
                    nbStates = 4,
                    stepPar0,
                    anglePar0,
                    formula = ~treecover.s,  
                    stepDist = "gamma",
                    angleDist = "vm",
                    verbose = 1)
HMM_ndvi_habitat <- fitHMM(vultsCOV,
                           nbStates = 4,
                           stepPar0,
                           anglePar0,
                           formula = ~ndvi.s + Forest + Savanna + Other,            
                           stepDist = "gamma",
                           angleDist = "vm",
                           verbose = 1) 
HMM_treec_habitat <- fitHMM(vultsCOV,
                            nbStates = 4,
                            stepPar0,
                            anglePar0,
                            formula = ~treecover.s + Forest + Savanna + Other,                                             
                            stepDist = "gamma",
                            angleDist = "vm",
                            verbose = 1)
HMM_elevation <- fitHMM(vultsCOV,
                        nbStates = 4,
                        stepPar0,
                        anglePar0,
                        formula = ~slope.s + easting.s,
                        stepDist = "gamma",
                        angleDist = "vm",
                        verbose = 1)
HMM_altitude <- fitHMM(vultsCOV,
                       nbStates = 4,
                       stepPar0,
                       anglePar0,
                       formula = ~height.s,    
                       angleDist = "vm",
                       verbose = 1)
HMM_height <- fitHMM(vultsCOV,
                     nbStates = 4,
                     stepPar0,
                     anglePar0,
                     formula = ~slope.s + easting.s + height.s,    
                     stepDist = "gamma",
                     angleDist = "vm",
                     verbose = 1)
HMM_slope_time <- fitHMM(vultsCOV,
                         nbStates = 4,
                         stepPar0,
                         anglePar0,
                         formula = ~slope.s + Time.sincedawn,
                         stepDist = "gamma",
                         angleDist = "vm",
                         verbose = 1)
HMM_habitat_slope_time <- fitHMM(vultsCOV,
                                 nbStates = 4,
                                 stepPar0,
                                 anglePar0,
                                 formula = ~Forest + Savanna + Other + slope.s + Time.sincedawn,                                
                                 stepDist = "gamma",
                                 angleDist = "vm",
                                 verbose = 1)
HMM_thermals <- fitHMM(vultsCOV,         
                       nbStates = 4,
                       stepPar0,
                       anglePar0,
                       formula = ~emissivity.s + albedoW.s + albedoB.s,  
                       stepDist = "gamma",
                       angleDist = "vm",
                       verbose = 1)
HMM_thermals_time <- fitHMM(vultsCOV,   
                            nbStates = 4,
                            stepPar0,
                            anglePar0,
                            formula = ~emissivity.s + albedoW.s + albedoB.s + Time.sincedawn,  
                            stepDist = "gamma",
                            angleDist = "vm",
                            verbose = 1)
HMM_thermals_time_altitude <- fitHMM(vultsCOV,    
                                     nbStates = 4,
                                     stepPar0,
                                     anglePar0,
                                     formula = ~emissivity.s + albedoW.s + albedoB.s + Time.sincedawn + height.s,  
                                     stepDist = "gamma",
                                     angleDist = "vm",
                                     verbose = 1)

## Test to see which is best fitted model (lowest AIC)
AIC(HMMno_cov, HMM_time, HMM_habitat, HMM_ndvi, HMM_treec, HMM_ndvi_habitat, HMM_treec_habitat, HMM_elevation, HMM_altitude, HMM_height, HMM_slope_time, HMM_habitat_slope_time, HMM_thermals, HMM_thermals_time, HMM_thermals_time_altitude)

##Investigate best model
print(HMM_thermals_time_altitude)
plot(HMM_thermals_time_altitude)

#Save for use in PPM scripts
save(HMMno_cov, HMM_time, HMM_habitat, HMM_ndvi, HMM_treec, HMM_ndvi_habitat, HMM_treec_habitat, HMM_elevation, HMM_altitude, HMM_height, HMM_slope_time, HMM_habitat_slope_time, HMM_thermals, HMM_thermals_time, HMM_thermals_time_altitude, file = "HMM_sub.all.Rdata")

