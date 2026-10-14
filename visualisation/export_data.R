# Regenerates visualisation/data/cma.geojson and stations.geojson
# Run from the repo root: Rscript visualisation/export_data.R
# Reproduces the same join/Gini pipeline as Code/DataCleaner2.R, limited to
# what the 3D map needs (geometry + AQHI + Gini + income), and writes it out
# as GeoJSON for the browser instead of keeping it in-memory.

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

main_dir <- getwd() # run this script from the repo root
out_dir <- file.path(main_dir, "visualisation/data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

regions <- c("atl", "ont", "pnr", "pyr")

d_AQHI <- lapply(setNames(regions, regions), function(region) {
  dir_path <- file.path(main_dir, "Data/AQHI", region)
  files <- list.files(path = dir_path, pattern = "\\.csv$", full.names = TRUE)
  rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
})

naps <- st_read(file.path(main_dir, "Data/aqhi_station.geojson"), quiet = TRUE)

valid_cgndb <- unlist(lapply(regions, function(region){
  setdiff(names(d_AQHI[[region]]), c("Date", "Hour (UTC)"))
}))

CMA <- st_read(file.path(main_dir, "Data/Spatial/CMA_shape.shp"), quiet = TRUE)
PR <- st_read(file.path(main_dir, "Data/Spatial/PR_shape_simplified.shp"), quiet = TRUE)
CMA <- CMA |> filter(PRUID != "24")
naps <- st_transform(naps, 3347)

naps_CMA <- naps |> filter(lengths(st_intersects(naps, CMA)) > 0)
naps_CMA <- st_join(naps_CMA, CMA["CMANAME"], join = st_intersects)

naps_df <- as.data.table(naps_CMA)[, .(cgndb, CMANAME, name,
  longitude = st_coordinates(geometry)[,1], latitude = st_coordinates(geometry)[,2])]
naps_df <- naps_df[naps_df$cgndb %in% valid_cgndb]
naps_df[, name := lapply(name, function(x) fromJSON(x)$en)]

d_Census <- fread(file.path(main_dir, "Data/Census/CMA_data_updated.csv"))

CMA_naps <- CMA |> filter(lengths(st_intersects(CMA, naps)) > 0)

naps_CMA <- st_join(naps_CMA, CMA_naps["CMAUID"], join = st_intersects)
d_CMA <- as.data.table(naps_CMA)[, .(cgndb, CMAUID, CMANAME)]

d_AQHI <- rbindlist(d_AQHI, use.names = TRUE, fill = TRUE, idcol = "region")
d_AQHI_l <- melt(d_AQHI, id.vars = c("region", "Date", "Hour (UTC)"), variable.name = "Station", value.name = "AQHI")

d <- d_AQHI_l |>
  inner_join(d_CMA, by = c("Station" = "cgndb"))

d <- d |>
  group_by(CMAUID, CMANAME, Date) |>
  summarise(mean_AQHI = mean(AQHI, na.rm = TRUE), .groups = "drop")

d <- d |>
  group_by(CMAUID, CMANAME) |>
  summarise(
    total_days = n(),
    high_AQHI_days = sum(mean_AQHI > 3, na.rm = TRUE),
    percent_high_AQHI = (high_AQHI_days/total_days)*100,
    average_AQHI = sum(mean_AQHI, na.rm = TRUE)/total_days,
    .groups = "drop"
  )

d_Census <- d_Census |>
  rename(CMAUID = `Census Metro Area ID`) |>
  mutate(CMAUID = str_pad(as.character(CMAUID), width = 3, side = "left", pad = "0"))

d <- d |>
  left_join(d_Census, by = "CMAUID")

d <- d |>
  left_join(CMA_naps, by = "CMAUID") |>
  select(-CMANAME.y) |>
  rename(CMANAME = CMANAME.x)

d <- st_as_sf(d)

# Positional rename matches Appendix A's I-vector: I0 median income, I1-I3 income thresholds
colnames(d)[8] <- "I0"
colnames(d)[9] <- "I1"
colnames(d)[10] <- "I2"
colnames(d)[11] <- "I3"

names(d)[names(d) == "Household Population by Income (Current Year) | Average Income ($), 2024"] <- "avg_income_dollars"
d$median_income_dollars <- d$I0

# Gini coefficient via monotone spline Lorenz curve approximation (see paper Appendix A)
placeholder_matrix <- matrix(NA, nrow = 2, ncol = 5, dimnames = list(c("x", "y"), paste0(1:5)))
d$IP <- replicate(nrow(d), placeholder_matrix, simplify = FALSE)

fIP <- function(i){
  mat <- d$IP[[i]]
  I0 <- d$I0[i]; I1 <- d$I1[i]; I2 <- d$I2[i]; I3 <- d$I3[i]
  xvals <- sort(c(0,0.5,1-(I1/100),1-(I1*I2/10000),1-(I1*I2*I3/1000000)))
  yvals <- sort(c(0,I0/300000,100000/300000,200000/300000,1))
  matrix(c(xvals, yvals), nrow = 2, byrow = TRUE)
}
d$IP <- lapply(seq_len(nrow(d)), fIP)

fix_IP <- function(mat){
  mat[2, ] <- pmin(mat[2, ], mat[1, ])
  return(mat)
}
d$IP <- lapply(d$IP, fix_IP)

d$P <- lapply(d$IP, function(mat){
  splinefun(mat[1,], mat[2,], method = "hyman")
})
d$E <- lapply(d$IP, function(mat){
  slope <- mat[2,5]/mat[1,5]
  function(x) slope*x
})
d$EP <- lapply(seq_len(nrow(d)), function(i){
  E <- d$E[[i]]; P <- d$P[[i]]
  function(x) E(x) - P(x)
})
d$Gini <- sapply(seq_len(nrow(d)), function(i){
  areaA <- integrate(d$EP[[i]], d$IP[[i]][1,1], d$IP[[i]][1,5])
  areaB <- integrate(d$E[[i]], d$IP[[i]][1,1], d$IP[[i]][1,5])
  areaA$value / areaB$value
})

d <- d |> filter(Gini > 0.02)       # excludes Wood Buffalo, Yellowknife (see paper, Data Description)
d <- d |> filter(average_AQHI > 1)  # excludes Brandon (insufficient AQHI coverage)

d <- d |> mutate(income_cat = factor(ntile(avg_income_dollars, 3), labels = c("Low", "Mid", "High")))

cat("Final CMA count:", nrow(d), "\n")

sf_use_s2(FALSE)

d_out <- d |>
  select(CMAUID, CMANAME, average_AQHI, percent_high_AQHI, Gini,
         avg_income_dollars, median_income_dollars, income_cat, LANDAREA) |>
  st_simplify(dTolerance = 1000, preserveTopology = TRUE) |>
  st_make_valid() |>
  st_transform(4326) |>
  st_make_valid() |>
  mutate(Gini = round(Gini, 4), average_AQHI = round(average_AQHI, 3),
         percent_high_AQHI = round(percent_high_AQHI, 1),
         avg_income_dollars = round(avg_income_dollars, 0))

st_write(d_out, file.path(out_dir, "cma.geojson"), delete_dsn = TRUE, quiet = TRUE)

naps_out <- naps_CMA |>
  st_transform(4326) |>
  select(cgndb, CMANAME)
st_write(naps_out, file.path(out_dir, "stations.geojson"), delete_dsn = TRUE, quiet = TRUE)

cat("Wrote cma.geojson and stations.geojson to", out_dir, "\n")
