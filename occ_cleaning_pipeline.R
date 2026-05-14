# --- Required Libraries --- #
library(sf)
library(terra)
library(dplyr)
library(purrr)
library(CoordinateCleaner)
library(tidyr)
library(ggplot2)
library(tidyterra)
library(rnaturalearth)

# --- PHASE I: Data Loading & Pre-processing --- #

# Use relative paths for GitHub compatibility
csv_paths  <- sort(list.files("data/raw_csv/", pattern = "\\.csv$", full.names = TRUE))
iucn_paths <- sort(list.files("data/iucn_shp/", pattern = "\\.shp$", full.names = TRUE))
exp_paths  <- sort(list.files("data/expert_shp/", pattern = "\\.shp$", full.names = TRUE))

process_raw_data <- function(df_raw) {
  # Helper to validate decimal precision
  is_valid_decimal <- function(x) {
    txt <- formatC(x, format = "f", digits = 10, drop0trailing = FALSE)
    txt <- sub("0+$", "", txt)
    grepl("^[-+]?[0-9]+\\.[0-9]{1,}$", txt)
  }
  
  df_raw %>%
    mutate(across(contains(c("ID", "Number", "Code")), as.character),
           decimalLatitude = as.numeric(decimalLatitude),
           decimalLongitude = as.numeric(decimalLongitude),
           year = as.numeric(year)) %>%
    mutate(
      .coord_ok = !is.na(decimalLatitude) & !is.na(decimalLongitude) & 
        is_valid_decimal(decimalLatitude) & is_valid_decimal(decimalLongitude),
      .year_ok  = !is.na(year),
      .sour_ok  = !is.na(basisOfRecord) & basisOfRecord != "FOSSIL_SPECIMEN"
    )
}

apply_coord_cleaner <- function(df) {
  if (!any(df$.coord_ok)) return(mutate(df, .cc_ok = FALSE))
  
  valid_coords <- df %>% filter(.coord_ok)
  bad_coords   <- df %>% filter(!.coord_ok)
  
  flags <- clean_coordinates(x = valid_coords, lon = "decimalLongitude", 
                             lat = "decimalLatitude", species = "species_co",
                             tests = c("zeros", "sea", "equal"))
  
  valid_coords$.cc_ok <- as.logical(flags$.summary)
  bad_coords$.cc_ok   <- FALSE
  bind_rows(valid_coords, bad_coords)
}

# --- PHASE II: Spatial Functions (Terra) --- #

# Calculates a dynamic buffer to increase area by a specific proportion (e.g., 10%)
calculate_dynamic_buffer <- function(iucn_p, exp_p, prop = 0.10, metric_crs = "EPSG:5070") {
  v_union <- aggregate(rbind(vect(iucn_p), vect(exp_p)))
  v_proj  <- makeValid(project(v_union, metric_crs))
  
  area_orig   <- expanse(v_proj)
  target_area <- area_orig * (1 + prop)
  
  # Optimization to find buffer distance 'd'
  opt_fun <- function(d) {
    b_temp <- aggregate(buffer(v_proj, width = d))
    return(expanse(makeValid(b_temp)) - target_area)
  }
  
  d_root <- uniroot(opt_fun, interval = c(0, 500000), extendInt = "yes", tol = 0.1)$root
  
  b_final <- project(aggregate(buffer(v_proj, width = d_root)), "EPSG:4326")
  return(list(buffer = b_final, dist_m = d_root, area_orig = area_orig))
}

spatial_validation <- function(df, iucn_p, exp_p, prop = 0.10) {
  res <- calculate_dynamic_buffer(iucn_p, exp_p, prop)
  pts <- vect(df, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  
  # Point-in-polygon test
  df$.union <- as.logical(is.related(pts, res$buffer, "intersects"))
  return(df)
}

# --- PHASE III: Execution Pipeline --- #

run_species_pipeline <- function(df, iucn_p, exp_p) {
  df_clean <- apply_coord_cleaner(process_raw_data(df))
  df_eval  <- df_clean %>% filter(.coord_ok & .cc_ok)
  
  if (nrow(df_eval) == 0) return(mutate(df_clean, .union = FALSE, .is_optimal = FALSE))
  
  df_val <- spatial_validation(df_eval, iucn_p, exp_p)
  
  df_clean %>%
    left_join(select(df_val, gbifID, .union), by = "gbifID") %>%
    mutate(.union = replace_na(.union, FALSE),
           .is_optimal = .coord_ok & .year_ok & .sour_ok & .cc_ok & .union)
}

# Execute for all species
raw_data_list <- lapply(csv_paths, read.csv)
results <- map(seq_along(raw_data_list), 
               ~run_species_pipeline(raw_data_list[[.x]], iucn_paths[[.x]], exp_paths[[.x]]))

full_results <- bind_rows(results)

# --- PHASE IV: Reporting & Export --- #

summary_table <- full_results %>%
  group_by(species_co) %>%
  summarise(n_total = n(),
            pct_optimal = round(100 * sum(.is_optimal) / n(), 1))

write.csv(summary_table, "output/filtering_summary.csv", row.names = FALSE)

# Export clean dataset
final_clean_data <- full_results %>% filter(.is_optimal == TRUE)
write.csv(final_clean_data, "output/cleaned_occurrence_data.csv", row.names = FALSE)