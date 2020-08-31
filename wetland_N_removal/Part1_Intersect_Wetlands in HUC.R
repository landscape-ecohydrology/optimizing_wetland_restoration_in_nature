# =====================================================================
# Wetlandscape Analysis - Preparing Wetland Tables for Analysis
# DATE: 26 February 2019
# 
# This code intersects the NWI wetlands with its HUC2 region
# This code was repeated for each HUC2 region, and appending the 
# appropriate states from the NWI (sometimes in several parts due to
# file size constraints)
# =====================================================================
cat("\014") ; rm(list=ls());            # clear console, clear plots, clear memory

library(sf)
library(dplyr)

# Get Watershed boundaries ----------------------------------------------------------
setwd("D:/GIS DATA FILES/Watershed Boundary Dataset/Reprojected/HUC8_15")
HUC_Boundary <-  st_read('HUC8_15.shp')

# Get wetlands by state -------------------------------------------------------------
setwd("D:/GIS DATA FILES/National Wetlands Inventory/NV_shapefile_wetlands")
wetland.shp <-  st_read('NV_Wetlands_South.shp')
a <- st_is_valid(wetland.shp)
wetland.shp$valid <- a
wetland.shp <- filter(wetland.shp, valid == T)
wetland.shp <- select(wetland.shp, -valid)
intersect_test <- st_intersection(wetland.shp, HUC_Boundary)


# Get wetlands
setwd("D:/GIS DATA FILES/National Wetlands Inventory/NM_shapefile_wetlands")
wetland.shp <-  st_read('NM_Wetlands.shp')
a <- st_is_valid(wetland.shp)
wetland.shp$valid <- a
wetland.shp <- filter(wetland.shp, valid == T)
wetland.shp <- select(wetland.shp, -valid)
intersect_test_b <- st_intersection(wetland.shp, HUC_Boundary)

# Get wetlands
setwd("D:/GIS DATA FILES/National Wetlands Inventory/UT_shapefile_wetlands")
wetland.shp <-  st_read('UT_Wetlands.shp')
a <- st_is_valid(wetland.shp)
wetland.shp$valid <- a
wetland.shp <- filter(wetland.shp, valid == T)
wetland.shp <- select(wetland.shp, -valid)
intersect_test_c <- st_intersection(wetland.shp, HUC_Boundary)


intersect_test_b <- select(intersect_test_b, -GLOBALID)


# Combine each state's wetlands to one giant table
HUC8_15a_Wetlands.shp <- rbind(intersect_test, intersect_test_b, intersect_test_c)


# Save Outputs ---------------------------------------------------------------------
setwd("D:/GIS DATA FILES/6 Wetlandscape Analysis/R_WetlandsByHUC_geospatial/OLD")
 saveRDS(HUC8_15a_Wetlands.shp, file = "HUC8_15a_Wetlands.rds")       # save for future in case i need geospatially referenced wetlands

# Save a pared down version of the table; just list of wetlands and area with HUC8 code + areas
col_index <- c('WETLAND_TY','ACRES','HUC8')
HUC8_15a_Wetlands.table <- st_drop_geometry(HUC8_15a_Wetlands.shp[,col_index ])
setwd("D:/GIS DATA FILES/6 Wetlandscape Analysis/R_WetlandsByHUC_table/OLD")
saveRDS(HUC8_15a_Wetlands.table, file = "HUC8_15a_Wetlands.rds")


