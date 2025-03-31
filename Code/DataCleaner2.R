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
library(splines2)

main_dir <- "D:/LeStoreDownload/ECON326_FinalProject-main/ECON326_FinalProject-main" #equivalent to getwd(), saved for ease

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

d_Census <- fread(file.path(main_dir, "Data/Census/CMA_data_updated.csv")) #import census data

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

colnames(d)[8] <- "I0" #naming by I matrix as described in our math analysis
colnames(d)[9] <- "I1"
colnames(d)[10] <- "I2"
colnames(d)[11] <- "I3"

placeholder_matrix <- matrix(NA, nrow = 2, ncol = 5, dimnames = list(c("x", "y"), paste0(1:5))) #placeholder matrix to be added to each CMA
d$IP <- replicate(nrow(d),placeholder_matrix,simplify = FALSE) #add the matrix to all rows of d

fIP <- function(i){ #iterable
  mat <- d$IP[[i]] #grab empty IP matrix
  I0 <- d$I0[i]
  I1 <- d$I1[i]
  I2 <- d$I2[i]
  I3 <- d$I3[i] #grab I vector values
  xvals <- sort(c(0,0.5,1-(I1/100),1-(I1*I2/10000),1-(I1*I2*I3/1000000))) #calculate IP matrix x values
  yvals <- sort(c(0,I0/300000,100000/300000,200000/300000,1)) #calculate IP matrix y values
  mat <- matrix(c(xvals, yvals), nrow = 2, byrow = TRUE) #assign values to matrix
  return(mat) #assign matrix to d$IP
}

d$IP <- lapply(seq_len(nrow(d)), fIP) #apply the function above to all rows in d

fix_IP <- function(mat){
  mat[2, ] <- pmin(mat[2, ],mat[1, ])  # force IPy to be less than or equal to IPx for income cases with erratic behavior
  return(mat)
}
d$IP <- lapply(d$IP, fix_IP)

d$P <- lapply(d$IP, function(mat){ #list apply to create P function for all rows in d
  splinefun(mat[1,],mat[2,],method = "hyman") #cubic spline function, using "hyman" method for monotonically increasing natural cubic spline
})

d$E <- lapply(d$IP, function(mat){ #list apply to create E function for all rows in d
  slope <- mat[2,5]/mat[1,5]  # IPy[5]/IPx[5]
  return(function(x) slope*x) #return function
})

d$EP <- lapply(seq_len(nrow(d)), function(i){ #list apply to create EP function for all rows in d
  E <- d$E[[i]]
  P <- d$P[[i]]
  return(function(x) E(x) - P(x)) #return function
})

d$Gini <- sapply(seq_len(nrow(d)), function(i){ #numerical values
  areaA <- integrate(d$EP[[i]],d$IP[[i]][1,1],d$IP[[i]][1,5])
  areaB <- integrate(d$E[[i]],d$IP[[i]][1,1],d$IP[[i]][1,5])
  return(areaA$value/areaB$value) #will assign this value to all rows of d$Gini
})

Lorenz <- function(i){ #function to plot Lorenz curves: some Gini coefficients were returning negative, need to check plots
  IPx <- d$IP[[i]][1,]
  IPy <- d$IP[[i]][2,]
  P <- d$P[[i]]
  E <- function(x) (IPy[5] / IPx[5]) * x
  x <- seq(IPx[1], IPx[5], length.out = 1000) #linspace of x
  y  <- P(x) #approximation spline
  Ey <- E(x) #equality line
  plot(IPx, IPy, pch = 19, main = paste("Lorenz Curve of", d$CMANAME[i]),
       xlab = "Population Share", ylab = "Income Share", ylim = c(0, 1))
  lines(x, y, col = "blue", lwd = 2)
  lines(x, Ey, col = "green", lwd = 2) #plotting
}

######################################

######################################
df_final <- d %>%
  st_drop_geometry() %>%  # 对于OLS，不需要geometry
  dplyr::select(
    CMAUID,
    CMANAME,
    Gini,
    average_AQHI,
    percent_high_AQHI,
    PRUID,
    `% Male Population by Age | Males, 2024`,
    `% Household Population by Visible Minority Status | Household Population For Visible Minority | Visible Minority Total, 2024`,
    `% Household Population 15 Years or Over by Marital Status | Household Population 15 Years Or Over | Married Or Living With A Common-Law Partner, 2024`,
    `% Household Population 15 Years or Over by Labour Force Activity | Employment Rate, 2024`,
    `% Household Population by Income (Current Year) | Household Population 15 Years Or Over For Income | With Income, 2024`
    # 
  ) %>%
  # 
  rename(
    pct_male        = `% Male Population by Age | Males, 2024`,
    pct_vismin      = `% Household Population by Visible Minority Status | Household Population For Visible Minority | Visible Minority Total, 2024`,
    pct_married     = `% Household Population 15 Years or Over by Marital Status | Household Population 15 Years Or Over | Married Or Living With A Common-Law Partner, 2024`,
    employment_rate = `% Household Population 15 Years or Over by Labour Force Activity | Employment Rate, 2024`,
    pct_with_income = `% Household Population by Income (Current Year) | Household Population 15 Years Or Over For Income | With Income, 2024`
  ) %>%
  # 
  filter(
    !is.na(Gini),
    !is.na(average_AQHI),
    !is.na(pct_male),
    !is.na(pct_vismin),
    !is.na(pct_married),
    !is.na(employment_rate),
    !is.na(pct_with_income)
  )

library(lmtest)
library(sandwich)

######################################

######################################
mod_main <- lm(
  Gini ~ average_AQHI +
    pct_male +
    pct_vismin +
    pct_married +
    employment_rate +
    pct_with_income,
  data = df_final
)

# 
summary(mod_main)  # OLS
coeftest(mod_main, vcov=vcovHC(mod_main, type="HC1"))  # 


mod_alt1 <- lm(
  Gini ~ percent_high_AQHI +
    pct_male +
    pct_vismin +
    pct_married +
    employment_rate +
    pct_with_income,
  data=df_final
)

coeftest(mod_alt1, vcov=vcovHC(mod_alt1, type="HC1"))


df_no_outliers <- df_final %>%
  filter(!CMANAME %in% c("Toronto", "Vancouver"))

mod_alt2 <- lm(
  Gini ~ average_AQHI +
    pct_male +
    pct_vismin +
    pct_married +
    employment_rate +
    pct_with_income,
  data=df_no_outliers
)

coeftest(mod_alt2, vcov=vcovHC(mod_alt2, type="HC1"))


df_final$region_factor <- as.factor(df_final$PRUID)

mod_alt3 <- lm(
  Gini ~ average_AQHI +
    pct_male +
    pct_vismin +
    pct_married +
    employment_rate +
    pct_with_income +
    region_factor,
  data=df_final
)

coeftest(mod_alt3, vcov=vcovHC(mod_alt3, type="HC1"))


library(stargazer)
######################################

######################################
stargazer(
  mod_main,
  mod_alt1,
  mod_alt2,
  mod_alt3,
  type = "text",   # 
  title = "Main vs. Robustness Checks",
  digits = 4,
  column.labels = c("Main", "Alt1", "Alt2", "Alt3")
)


# 
# install.packages("car")
library(car)

# 
vif_values_main <- vif(mod_main)
vif_values_main



# 
df_final$average_AQHI_sq <- df_final$average_AQHI^2

# 
mod_poly <- lm(
  Gini ~ average_AQHI + average_AQHI_sq +
    pct_male +
    pct_vismin +
    pct_married +
    employment_rate +
    pct_with_income,
  data = df_final
)

# 
summary(mod_poly)
library(lmtest); library(sandwich)
coeftest(mod_poly, vcov = vcovHC(mod_poly, type="HC1"))


df_final$log_AQHI <- log(df_final$average_AQHI)

mod_log <- lm(
  Gini ~ log_AQHI +
    pct_male +
    pct_vismin +
    pct_married +
    employment_rate +
    pct_with_income,
  data = df_final
)

coeftest(mod_log, vcov=vcovHC(mod_log, type="HC1"))

