library(leaflet)
library(terra)
library(PoPS)
library(terra)
library(folderfun)
library(doParallel)
library(raster)
# library(plyr)

setff("In", "C:/Users/cmjone25/Desktop/SOD_OR/")
infected_file <- ffIn("End of Year Infections/end_inf_2021_eu1.tif")
eu1_2021 <- rast(infected_file)

crs(eu1_2021)
# raster::crs(eu1_2021)
crs(eu1_2021) <- "EPSG:3857"
s <- raster(eu1_2021)
crs(s) <- CRS('+init=EPSG:26710')
plet(eu1_2021)
s[s == 0] <- NA

pal <- colorNumeric(palette = "YlOrRd", values(s),
                    na.color = "transparent")

leaflet() %>% addProviderTiles(providers$Esri.WorldGrayCanvas) %>%
  addRasterImage(s, colors = pal) %>%
  addLegend(pal = pal, values = values(s), title = "Infections") 
# Layers control
addLayersControl(
  overlayGroups = c("Quakes", "Outline"),
  options = layersControlOptions(collapsed = FALSE)
)

t <- raster(eu1_sim)
crs(t) <- CRS('+init=EPSG:26710')
t[t == 0] <- NA


# pal <- colorNumeric(palette = "YlOrRd", values(t), na.color = "transparent")
leaflet() %>% addProviderTiles(providers$Esri.WorldGrayCanvas) %>%
  addRasterImage(t, colors = pal) %>%
  addLegend(pal = pal, values = values(t), title = "Infections") 
