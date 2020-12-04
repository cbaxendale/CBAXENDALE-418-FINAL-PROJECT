install.packages("spgwr")
install.packages("spatstat")
install.packages("tmap")
install.packages("gstat")
install.packages("sf")
install.packages("raster")
install.packages("rgdal")
install.packages("e1071")
install.packages("spdep")

library(spgwr)
library(spatstat)
library(tmap)
library(gstat)
library(sf)
library(raster)
library(rgdal)
library(e1071)
library(spdep)
library(sp)

dir <- "C:/Users/cbaxe/Desktop/Final"
setwd(dir)
getwd()

###DATA PREP###
#Read in PM2.5 data
pm2.5 <- readOGR(".", "Pm25Sample") 
pm2.5 <- spTransform(pm2.5, CRS("+init=epsg:26910"))

#Read in census income data
income <- read.csv("./Income.csv")

#Select only ID and Income columns:
colnames(income) <- c("DAUID", "Income")

#Read in dissemination tract shapefile:
census.tracts <- readOGR(".", "BC_DA") 

#Merge income and dissemination data:
income.tracts <- merge(census.tracts,income, by = "DAUID") 

#Determine the number of columns in the dataframe:
nrow(income.tracts)

#Remove NA values:
income.tracts <- income.tracts[!is.na(income.tracts$Income),]

#Reproject the data:
income.tracts <- spTransform(income.tracts, CRS("+init=epsg:26910"))

tmaptools::palette_explorer()

map_Income <- tm_shape(income.tracts) +
  tm_polygons(col = "Income",
              title = "Median Income",
              style = "jenks",
              palette = "viridis", n = 6) +
  tm_legend(legend.position = c("LEFT", "BOTTOM"))

map_Income

###Descriptive Stats###

meanincome <- mean(income$Income, na.rm = TRUE)

###Income Segregation###
#Global Moran's I
vri.nb <- poly2nb(income.tracts)
vri.net <- nb2lines(vri.nb, coords=coordinates(income.tracts))
crs(vri.net) <- crs(income.tracts)

tm_shape(income.tracts) + tm_borders(col='lightgrey') + 
  tm_shape(vri.net) + tm_lines(col='red')


vri.lw <- nb2listw(vri.nb, zero.policy = TRUE, style = "W")
print.listw(vri.lw, zero.policy = TRUE)


mi <- moran.test(income.tracts$Income, vri.lw, zero.policy = TRUE)
mi

moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(vri.lw)

mI <- mi$estimate[[1]]
eI <- mi$estimate[[2]]
var <- mi$estimate[[3]]

z <- (mI-eI)/(sqrt(var))
z

#Local Moran's I
lisa.test <- localmoran(income.tracts$Income, vri.lw, zero.policy = TRUE)

income.tracts$Ii <- lisa.test[,1]
income.tracts$E.Ii<- lisa.test[,2]
income.tracts$Var.Ii<- lisa.test[,3]
income.tracts$Z.Ii<- lisa.test[,4]
income.tracts$P<- lisa.test[,5]

moran.plot(income.tracts$Income, vri.lw, zero.policy=TRUE, spChk=NULL, labels=NULL, xlab="Median Income", 
           ylab="Spatially Lagged Income", quiet=NULL)


###Spatial Interpolation of PM2.5
#Create an empty grid where n is the total number of cells
#Create a grid called grd to use in your interpolation
grd <- as.data.frame(spsample(pm2.5, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  
fullgrid(grd)    <- TRUE  
proj4string(grd) <- proj4string(income.tracts)


P.idw <- gstat::idw(PM25 ~ 1, pm2.5, newdata=grd, idp=4)
r       <- raster(P.idw)
r.m     <- mask(r, income.tracts)

tm_shape(r.m) + 
  tm_raster(n=10,palette = "YlOrRd",
            title="Predicted PM2.5 Concentration (ppm)") + 
  tm_shape(income.tracts) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

###Point Pattern Analysis to Determine if PM2.5 is random
kma <- pm2.5
kma$x <- coordinates(kma)[,1]
kma$y <- coordinates(kma)[,2]


zd <- zerodist(kma)
zd

kma <- remove.duplicates(kma)

kma.ext <- as.matrix(extent(kma)) 

window <- as.owin(list(xrange = kma.ext[1,], yrange = kma.ext[2,]))

kma.ppp <- ppp(x = kma$x, y = kma$y, window = window)

###NND
nearestNeighbour <- nndist(kma.ppp)

nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))

colnames(nearestNeighbour) = "Distance"
n <- nrow(nearestNeighbour)

nnd = (sum(nearestNeighbour$Distance))/n

studyArea <- sum(area(income.tracts))
pointDensity <- n/studyArea
studyArea

r.nnd = 1/(2*(sqrt(pointDensity)))

d.nnd = (1.07453/(sqrt(pointDensity)))

R = nnd/r.nnd

SE.NND <- 0.26136/(sqrt(n*pointDensity))

z = (nnd-r.nnd)/SE.NND
z

##QUADRAT ANALYSIS

quads <- 10

qcount <- quadratcount(kma.ppp, nx = quads, ny = quads)

plot(kma.ppp, pch = "+", cex = 0.5)
plot(qcount, add = T, col = "red")

qcount.df <- as.data.frame(qcount)

qcount.df <- plyr::count(qcount.df,'Freq')

colnames(qcount.df) <- c("x","f")

sum.f.x2 <- sum((qcount.df$f)*((qcount.df$x)^2))

M <- sum(qcount.df$f)

N <- sum(qcount.df$x)

sum.fx.2 <- sum(((qcount.df$f)*(qcount.df$x))^2)

VAR <- ((sum.f.x2)-((sum.fx.2)/M))/M-1

MEAN <- N/M

VMR <- VAR/MEAN

ChiSq <- (VMR*(M-1))
ChiSq

VAR

###COMBINING DATA
#Convert your interpolation into a raster and map it:
r <- raster(r.m)
sufaceMap <- tm_shape(r.m) + 
  tm_raster(n=5,palette = "viridis",
            title="PM 2.5 \n(in ppm)") +
  tm_shape(pm2.5) + tm_dots(size=0.2)

sufaceMap

income.tracts$Pm2.5 <- round(extract(r.m, income.tracts, fun = mean)[,1], 5)

###LINEAR REGRESSION

income.tracts.no0 <-  income.tracts[which(income.tracts$Pm2.5 > 0), ]

plot(income.tracts.no0$Income~income.tracts.no0$Pm2.5)

lm.model <- lm(income.tracts.no0$Income~income.tracts.no0$Pm2.5)

plot(income.tracts.no0$Income~income.tracts.no0$Pm2.5)
abline(lm.model, col = "red")

summary(lm.model)

income.tracts.no0$predictlm <- lm.model$fitted.values

income.tracts.no0$residuals <- residuals.lm(lm.model)

head(income.tracts.no0)

map_resid <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "residuals",
              title = "Median Income",
              style = "fisher",
              palette = "viridis", n = 6,
              midpoint = 0, border.alpha = 0.1) +
  tm_legend(legend.outside = TRUE)

map_resid

###GEOGRAPHICALL WEIGHTED REGRESSION
income.tracts.no0.coords <- sp::coordinates(income.tracts.no0)

head(income.tracts.no0.coords)
income.tracts.no0$X <- income.tracts.no0.coords[,1]
income.tracts.no0$Y <- income.tracts.no0.coords[,2]

GWRbandwidth <- gwr.sel(income.tracts.no0$Income~income.tracts.no0$Pm2.5, 
                        data=income.tracts.no0, coords=cbind(income.tracts.no0$X,income.tracts.no0$Y),adapt=T) 

gwr.model = gwr(income.tracts.no0$Income~income.tracts.no0$Pm2.5, 
                data=income.tracts.no0, coords=cbind(income.tracts.no0$X,income.tracts.no0$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

gwr.model

results<-as.data.frame(gwr.model$SDF)
head(results)

income.tracts.no0$localr <- results$localR2

#Create choropleth map of r-square values
map_r2 <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "localr",
              title = "R2 values",
              style = "fisher",
              palette = "viridis", n = 6)
map_r2


income.tracts.no0$coeff <- results$income.tracts.no0.Pm2.5
map_coef <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "coeff",
              title = "Coefficients",
              style = "fisher",
              palette = "viridis", n = 6,
              midpoint = 0, border.alpha = 0.1)
map_coef
