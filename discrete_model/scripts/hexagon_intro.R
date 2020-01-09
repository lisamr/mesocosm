require(sp)
data(meuse.riv)
meuse.sr = SpatialPolygons(list(Polygons(list(Polygon(meuse.riv)), "x")))
plot(meuse.sr)

library(rgeos)
meuse.large = gBuffer(meuse.sr, width = 2000)
HexPts <-spsample(meuse.large, type="hexagonal", cellsize=1000)
HexPols <- HexPoints2SpatialPolygons(HexPts)
plot(HexPols[meuse.sr,], add=TRUE, col='blue')
plot(HexPols, add=TRUE)
plot(meuse.sr, add=TRUE)


#check out the hex polys
coordinates(HexPols)
plot(coordinates(HexPols)[,1], coordinates(HexPols)[,2])
