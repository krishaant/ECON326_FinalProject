library(tidyverse) #importing libraries, not all used yet
library(gtools)
library(modelsummary)
library(ggeffects)
library(data.table)
library(sf)
library(jsonlite)
library(data.table)
library(tmap)
library(ggplot2)
library(terra)
library(dplyr)
library(tidyr)

main_dir <- "C:/Users/zaneb/Documents/Vault/UBC/24W/Term 2/ECON 326/FinalProject" #equivalent to getwd(), saved for ease

regions <- c("atl", "ont", "pnr", "pyr") #the four areas of data collection in AQHI files

d_AQHI <- lapply(setNames(regions, regions), function(region) { #lapply applies a function to a list
  dir_path <- file.path(main_dir, "Data/AQHI", region) #folder of data
  files <- list.files(path = dir_path, pattern = "\\.csv$", full.names = TRUE) #list of files in the folder
  rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE) #reading all of the folders, then binding them together under separate lists in d_AQHI
})

naps <- st_read(file.path(main_dir, "Data/aqhi_station.geojson")) #read the geojson guide to the different research locations

valid_cgndb <- unlist(lapply(regions, function(region){
  setdiff(names(d_AQHI[[region]]), c("Date", "Hour (UTC)")) #generate a list of cgndb values in d_AQHI
}))

CMA <- st_read(file.path(main_dir, "Data/Spatial/CMA_shape.shp")) #import CMA

PR <- st_read(file.path(main_dir, "Data/Spatial/PR_shape_simplified.shp")) #import PR

CMA <- CMA |> filter(PRUID != "24") #filter out QC

naps <- st_transform(naps, 3347) #CRS correction

naps_CMA <- naps |> filter(lengths(st_intersects(naps, CMA)) > 0) #create a naps dataframe that includes only sites within a CMA
naps_CMA <- st_join(naps_CMA, CMA["CMANAME"], join = st_intersects) #add CMA value

naps_df <- as.data.table(naps_CMA)[, .(cgndb, CMANAME, name, longitude = st_coordinates(geometry)[,1], latitude = st_coordinates(geometry)[,2])] #prep for merge, grab name & longitude and latitude values

naps_df <- naps_df[naps_df$cgndb %in% valid_cgndb] #remove irrelevant stations

naps_df[, name := lapply(name, function(x) fromJSON(x)$en)] #adjust names to english

cgndb_key <- setNames(naps_df$name, naps_df$cgndb) #create a key of place names and cgndb values

placenames <- unlist(lapply(regions, function(region){
  setdiff(names(d_AQHI[[region]]), c("Date", "Hour (UTC)")) #making a list of cgndb now in d_AQHI
}))

d_Census <- fread(file.path(main_dir, "Data/Census/CMA_Data.csv")) #import census data

CMA_naps <- CMA |> filter(lengths(st_intersects(CMA, naps)) > 0) #create a CMA dataframe that includes only sites with a naps

CMA_nonaps <- CMA |> filter(lengths(st_intersects(CMA, naps)) == 0) #create a CMA with sites that don't have naps

naps_CMA <- st_join(naps_CMA, CMA_naps["CMAUID"], join = st_intersects)#add CMAID values to naps

d_CMA <- as.data.table(naps_CMA)[, .(cgndb, CMAUID, CMANAME)] #short df of cdngb values with CMA locations & ID

aqhi_cgndb <- unlist(lapply(regions, function(region){
  setdiff(names(d_AQHI[[region]]), c("Date", "Hour (UTC)")) #making a list of placenames now in d_AQHI
}))

d_AQHI <- rbindlist(d_AQHI, use.names = TRUE, fill = TRUE, idcol = "region") #break apart regional grouping

d_AQHI_l <- melt(d_AQHI, id.vars = c("region", "Date", "Hour (UTC)"), variable.name = "Station", value.name = "AQHI") #shift from wide to long

d <- d_AQHI_l |>
  inner_join(d_CMA, by = c("Station" = "cgndb")) #add CMA data: stations will be labelled by CMA

d <- d |>
  group_by(CMAUID, CMANAME, Date) |>
  summarise(mean_AQHI = mean(AQHI, na.rm = TRUE), .groups = "drop") #regroup by average daily AQHI for a CMA at a given day

d <- d |>
  group_by(CMAUID, CMANAME) |> #group by CMA
  summarise(
    total_days = n(), #total days with data
    high_AQHI_days = sum(mean_AQHI > 3, na.rm = TRUE), #number of days with high AQHI
    percent_high_AQHI = (high_AQHI_days/total_days)*100, #% of days with high AQHI
    average_AQHI = sum(mean_AQHI, na.rm = TRUE)/total_days #average of daily AQHI averages
  )

d_Census <- d_Census |>
  rename(CMAUID=`Census Metro Area ID`) #change for merge

d_Census <- d_Census |>
  mutate(CMAUID = str_pad(as.character(CMAUID), width = 3, side = "left", pad = "0")) #fix formatting

d <- d |>
  left_join(d_Census, by = "CMAUID") #combine the census and AQHI data!

d <- d |>
  left_join(CMA_naps, by = "CMAUID") |> #add the geographic data!
  select(-CMANAME.y) |> 
  rename(CMANAME = CMANAME.x) #removing doubled CMANAME column

d <- st_as_sf(d) #make d recognized as spatial again