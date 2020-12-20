library(here)
library(data.table)
library(lubridate)
library(stringr)
library(dplyr)
library(sp)
library(raster)
library(ggplot2)
library(suncalc)
library(tidyr)
library(purrr)
library(sf)
library(matlib)
library(tmap)

ag_poly <- st_read(here::here("data", "crop_map", "i15_Crop_Mapping_2016.shp")) %>%
  st_transform(4326)

#clip and flatten ag polygon to stady extent for faster computation
study_extent <- st_read(here::here("data", "Fresno_polygon", "Extent_polygon.shp"))
ag_poly <- st_intersection(ag_poly, study_extent)
built <- ag_poly %>% filter(Crop2016 == "Urban")
ag_poly <- ag_poly %>% filter(Crop2016 != "Urban")

ndvi2018 <- raster(here::here("data", "Landsat_fresno_july", "NDVI", paste0("NDVI", 2018, ".tif")))
ndvi2019 <- raster(here::here("data", "Landsat_fresno_july", "NDVI", paste0("NDVI", 2019, ".tif")))
ndvi2020 <- raster(here::here("data", "Landsat_fresno_july", "NDVI", paste0("NDVI", 2020, ".tif")))

getrow <- function(file, kind, correction = "corrected"){
  print(file)
  
  date <- ymd(substr(file, 14, 23))
  hhmmss <- str_extract(file, regex('[0-9]{2}:{1}[0-9]{2}:{1}[0-9]{2}'))
  dt <- ymd_hms(paste(date, hhmmss), tz = "America/Los_Angeles")
  year <- year(date)
  
  ras <- raster(paste0(here::here("data", "June-Sept_2018-2020", kind, correction), "/",
                       file))
  
  mean <- cellStats(ras, "mean")
  
  meanfrommask <- function(poly, inv = FALSE){
    for (i in 1:length(inv)){
      ras <- raster::mask(ras,
                          poly[[i]],
                          inverse = inv[i])
    }
    y <- cellStats(ras, "mean", na.rm = TRUE)
  }
  
  #AG
  print("ag")
  ag <- meanfrommask(list(ag_poly), FALSE)
  nag <- meanfrommask(list(ag_poly), TRUE)
  
  #Veg
  print("veg")
  vegmask <- resample(get(paste0("ndvi", year)), ras)
  vegmask[vegmask[] < .5] <- NA
  vegmask[vegmask[] >= .5] <- 1
  
  veg <- meanfrommask(list(vegmask), FALSE)
  nveg <- meanfrommask(list(vegmask), TRUE)
  
  #Built
  print("built")
  urban <- meanfrommask(list(built), FALSE)
  
  print("mix")
  urb_veg <- meanfrommask(list(vegmask, built), c(FALSE, FALSE))
  urb_nveg <- meanfrommask(list(vegmask, built), c(TRUE, FALSE))
  ag_veg <- meanfrommask(list(vegmask, ag_poly), c(FALSE, FALSE))
  ag_nveg <- meanfrommask(list(vegmask, ag_poly), c(TRUE, FALSE))
  
  print("row")
  row <- data.frame(date, dt, kind, correction, mean, ag, nag, veg, nveg, urban, urb_veg, urb_nveg, ag_veg, ag_nveg)
}


get_time_series_df <- function(kind, correction = "corrected"){
  files <- list.files(here::here("data", "June-Sept_2018-2020", kind, correction))
  rows <- lapply(files, getrow, kind, correction)
  rbindlist(rows)
}

# biting_df <- get_time_series_df("tarsalis_biting_rate", "corrected")
# write.csv(biting_df, here::here("data", "June-Sept_2018-2020", "tars_biting_df.csv"), row.names = FALSE)
# transmission_df <- get_time_series_df("tarsalis_transmission", "corrected")
# write.csv(transmission_df, here::here("data", "June-Sept_2018-2020", "tars_transmission_df.csv"), row.names = FALSE)
# 
# temperature_df <- get_time_series_df("temperature", "corrected")
# write.csv(temperature_df, here::here("data", "June-Sept_2018-2020", "temperature_df.csv"), row.names = FALSE)
lst_df <- get_time_series_df("temperature", "not_corrected")
write.csv(lst_df, here::here("data", "June-Sept_2018-2020", "lst_df.csv"), row.names = FALSE)
