  # ==============================================================================
  # SPECIES DISTRIBUTION MODELLING (SDM) - SPECIES OCCURRENCE FILTERING PIPELINE
  # ==============================================================================
# Description: This script loads raw occurrence data (CSV), applies spatial and
# temporal quality filters, validates coordinates against dynamic biogeographical
# buffer polygons (IUCN & Expert maps), and exports a curated database.
#
# INSTRUCTIONS FOR USERS:
# 1. Directory Setup: Ensure the following subfolders exist in your working directory:
#    - "base/"     : Contains raw CSV occurrence files (one file per species/group).
#    - "pol_iucn/" : Contains IUCN range polygon shapefiles (.shp).
#    - "pol_exp/"  : Contains Expert range polygon shapefiles (.shp).
# 2. File Alignment: Ensure files in "base/", "pol_iucn/", and "pol_exp/" follow
#    an identical alphabetical sorting order to match species correctly.
# 3. Output: Generates "curated_species_database.csv" in UTF-8 format.

# =================================================================================

# --- DEPENDENCIES ---
library(sf)
library(terra)
library(dplyr)
library(purrr)
library(CoordinateCleaner)
library(tidyr)
library(scales)
library(readr)

# --- WORKING DIRECTORY SETUP ---
setwd("")

# ==============================================================================
# PHASE I: ENVIRONMENT SETUP AND DATA LOADING
# ==============================================================================

# List and sort input file paths
csv_paths    <- sort(list.files("base/", pattern = "\\.csv$", full.names = TRUE))
iucn_paths   <- sort(list.files("pol_iucn/", pattern = "\\.shp$", full.names = TRUE))
expert_paths <- sort(list.files("pol_exp/", pattern = "\\.shp$", full.names = TRUE))

# Load CSV files enforcing character type to avoid type-mismatch errors across files
raw_csv_list <- lapply(csv_paths, function(file) {
  read.csv(file, colClasses = "character", stringsAsFactors = FALSE)
})

# Function to clean and flag basic attribute quality
process_raw_data <- function(df_raw, sp_index = 1) {
  
  # Helper function: check if coordinates have at least 1 explicit decimal
  has_valid_decimals <- function(x) {
    txt <- formatC(x, format = "f", digits = 10, drop0trailing = FALSE)
    txt <- sub("0+$", "", txt)
    grepl("^[-+]?[0-9]+\\.[0-9]{1,}$", txt)
  }
  
  # Standardize species name column
  # 'species_co' contains the verified species name based on manual taxonomic curation.
  # Fallback to 'species' or NA if 'species_co' is missing in raw input.
  if (!"species_co" %in% names(df_raw)) {
    if ("species" %in% names(df_raw)) {
      df_raw$species_co <- df_raw$species
    } else {
      df_raw$species_co <- NA_character_
    }
  }
  
  df_clean <- df_raw %>%
    mutate(
      species_co       = trimws(as.character(species_co)),
      decimalLatitude  = as.numeric(decimalLatitude),
      decimalLongitude = as.numeric(decimalLongitude),
      year             = if("year" %in% names(.)) as.numeric(year) else NA_real_,
      row_uid          = paste0("sp", sp_index, "_row", row_number())
    ) %>%
    filter(!is.na(species_co) & species_co != "" & species_co != "NA") %>%
    mutate(  
      # Flag 1: Valid coordinate ranges and decimal precision
      .coord_ok = !is.na(decimalLatitude) & !is.na(decimalLongitude) & 
        decimalLatitude >= -90 & decimalLatitude <= 90 & 
        decimalLongitude >= -180 & decimalLongitude <= 180 &
        has_valid_decimals(decimalLatitude) & 
        has_valid_decimals(decimalLongitude),
      # Flag 2: Presence of temporal information
      .year_ok  = !is.na(year),
      # Flag 3: Exclusion of non-suitable record sources (e.g., fossils)
      .sour_ok  = !is.na(basisOfRecord) & basisOfRecord != "FOSSIL_SPECIMEN"
    )
  
  return(df_clean)
}

# Function to apply automated geographic cleaning algorithms via CoordinateCleaner
apply_coordinate_cleaner <- function(df_raw) {
  if (!any(df_raw$.coord_ok)) {
    df_raw$.cc_ok <- FALSE 
    return(df_raw) 
  } 
  
  coords_ok <- df_raw %>% filter(.coord_ok) 
  
  flags <- clean_coordinates( 
    x = coords_ok, 
    lon = "decimalLongitude", 
    lat = "decimalLatitude", 
    species = "species_co", 
    tests = c("zeros", "seas", "equal"),
    value = "spatialvalid",
    verbose = FALSE
  ) 
  
  coords_ok$.cc_ok <- as.logical(flags$.summary) 
  
  df_raw %>%
    left_join(
      coords_ok %>% select(row_uid, .cc_ok),
      by = "row_uid"
    ) %>%
    mutate(.cc_ok = ifelse(is.na(.cc_ok), FALSE, as.logical(.cc_ok)))
}

# Execute Phase I transformations
processed_list <- map(seq_along(raw_csv_list), ~ process_raw_data(raw_csv_list[[.x]], sp_index = .x))
post_cc_list   <- lapply(processed_list, apply_coordinate_cleaner)

# ==============================================================================
# PHASE II: SPATIAL ANALYSIS AND DYNAMIC BUFFERS (Terra)
# ==============================================================================

# Calculates an area-proportional dynamic buffer around range maps
calculate_dynamic_buffer <- function(iucn_path, expert_path, prop = 0.10, metric_crs = "EPSG:5070") {
  v_iucn   <- vect(iucn_path)
  v_expert <- vect(expert_path)
  
  # Merge baseline range geometries
  merged_original <- aggregate(rbind(v_iucn, v_expert))
  projected_union <- project(merged_original, metric_crs)
  projected_union <- makeValid(projected_union) 
  
  base_area <- expanse(projected_union)
  if(base_area == 0) stop("The baseline extent area is 0. Check shapefile topology.")
  
  target_area <- base_area * (1 + prop)
  
  # Objective function to solve for exact buffer distance yielding target expansion
  optimization_function <- function(d) {
    temp_buffer <- aggregate(buffer(projected_union, width = d))
    temp_buffer <- makeValid(temp_buffer)
    return(expanse(temp_buffer) - target_area)
  }
  
  buffer_distance <- uniroot(optimization_function, 
                             interval = c(0, 500000), 
                             extendInt = "yes", 
                             tol = 0.1)$root
  
  final_buffer <- aggregate(buffer(projected_union, width = buffer_distance))
  final_buffer <- makeValid(final_buffer)
  
  return(list(
    buffer = project(final_buffer, "EPSG:4326")
  ))
}

# Validates occurrence point intersections against calculated spatial buffer
validate_spatial_bounds <- function(df_species, iucn_path, expert_path, prop = 0.10, metric_crs = "EPSG:5070") {
  buffer_results <- calculate_dynamic_buffer(iucn_path, expert_path, prop, metric_crs)
  final_buffer   <- buffer_results$buffer 
  
  species_points <- vect(df_species, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  points_inside  <- is.related(species_points, final_buffer, "intersects")
  
  df_species$.union <- as.logical(points_inside)
  return(df_species)
}

# ==============================================================================
# PHASE III: PIPELINE EXECUTION AND EXPORT
# ==============================================================================

process_species_pipeline <- function(df_post_cc, iucn_path, expert_path, prop = 0.10) {
  
  # BIOGEOGRAPHIC EVALUATION: Applied to all records with valid coordinates (.coord_ok)
  # regardless of CoordinateCleaner or temporal/source filter status.
  coords_subset <- df_post_cc %>% filter(.coord_ok)
  
  if (nrow(coords_subset) == 0) {
    df_post_cc$.union      <- FALSE
    df_post_cc$.is_optimal <- FALSE
    return(df_post_cc)
  }
  
  validated_df <- validate_spatial_bounds(coords_subset, iucn_path, expert_path, prop)
  
  df_post_cc %>%
    left_join(validated_df %>% select(row_uid, .union), by = "row_uid") %>%
    mutate(
      .union      = ifelse(is.na(.union), FALSE, as.logical(.union)),
      .is_optimal = as.logical(.coord_ok & .year_ok & .sour_ok & .cc_ok & .union)
    )
}

# Execute full pipeline across species list
curated_results <- map(
  seq_along(post_cc_list),
  ~ process_species_pipeline(post_cc_list[[.x]], iucn_paths[[.x]], expert_paths[[.x]])
)

# Combine datasets and derive final taxonomy features
all_final <- bind_rows(curated_results) %>%
  mutate(
    species_name = as.character(species_co),
    genero       = sub(" .*", "", species_name)
  )

# Export complete curated database using readr to guarantee UTF-8 encoding
readr::write_csv(all_final, "sciuridae_database.csv")