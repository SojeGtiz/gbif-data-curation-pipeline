# --- Required Libraries ---
library(sf)
library(terra)
library(dplyr)
library(purrr)
library(CoordinateCleaner)
library(tidyr)

# --- PHASE I: Environment Setup and Data Loading ---
csv_paths    <- sort(list.files("data/raw_csv/", pattern = "\\.csv$", full.names = TRUE))
iucn_paths   <- sort(list.files("data/iucn_shp/", pattern = "\\.shp$", full.names = TRUE))
expert_paths <- sort(list.files("data/expert_shp/", pattern = "\\.shp$", full.names = TRUE))

# Load raw occurrences datasets into a list
raw_csv_list <- lapply(csv_paths, read.csv)

# Define standard core target columns to retain
target_columns <- c("gbifID", "species_co", "decimalLatitude", "decimalLongitude", "year", "month", "basisOfRecord", "identifiedBy", "institutionCode", "locality")

# Parse, filter, and structure raw occurrence inputs
process_raw_data <- function(df_raw) {
  
  # Structural helper to evaluate geographic coordinate decimal precision
  has_valid_decimals <- function(x) {
    txt <- formatC(x, format = "f", digits = 10, drop0trailing = FALSE)
    txt <- sub("0+$", "", txt)
    grepl("^[-+]?[0-9]+\\.[0-9]{1,}$", txt)
  }
  
  df_raw %>%
    select(intersect(names(.), target_columns)) %>%
    mutate(
      gbifID           = as.character(gbifID),
      decimalLatitude  = as.numeric(decimalLatitude),
      decimalLongitude = as.numeric(decimalLongitude),
      year             = as.numeric(year)
    ) %>%
    mutate(  
      .coord_ok = !is.na(decimalLatitude) & !is.na(decimalLongitude) & 
        has_valid_decimals(decimalLatitude) & 
        has_valid_decimals(decimalLongitude),
      .year_ok  = !is.na(year),
      .sour_ok  = !is.na(basisOfRecord) & basisOfRecord != "FOSSIL_SPECIMEN"
    )
}

# Flag point dataset coordinate quality records using CoordinateCleaner
apply_coordinate_cleaner <- function(df_raw) {
  if (!any(df_raw$.coord_ok)) {
    df_raw$.cc_ok <- FALSE 
    return(df_raw) 
  } 
  
  coords_ok  <- df_raw %>% filter(.coord_ok) 
  coords_bad <- df_raw %>% filter(!.coord_ok) 
  
  flags <- clean_coordinates( 
    x = coords_ok, 
    lon = "decimalLongitude", 
    lat = "decimalLatitude", 
    species = "species_co", 
    tests = c("zeros", "sea", "equal") 
  ) 
  
  coords_ok$.cc_ok  <- as.logical(flags$.summary) 
  coords_bad$.cc_ok <- logical(nrow(coords_bad)) 
  
  bind_rows(coords_ok, coords_bad) 
}

# Execute Phase I cleaning workflows
processed_list <- lapply(raw_csv_list, process_raw_data)
post_cc_list    <- lapply(processed_list, apply_coordinate_cleaner)

# --- PHASE II: Spatial Functions (Terra) ---

# Calculates an automated range buffer optimized for a targeted area expansion
calculate_dynamic_buffer <- function(iucn_path, expert_path, prop = 0.10, metric_crs = "EPSG:5070") {
  
  # 1. Load range layout shapefiles as SpatVectors
  v_iucn   <- vect(iucn_path)
  v_expert <- vect(expert_path)
  
  # 2. Merge shapefiles and dissolve internal borders to handle fragmentation
  merged_original <- aggregate(rbind(v_iucn, v_expert))
  
  # Project layout to metric CRS for optimized area/buffer math
  projected_union <- project(merged_original, metric_crs)
  projected_union <- makeValid(projected_union) 
  
  # 3. Compute baseline extent metrics
  base_area <- expanse(projected_union)
  if(base_area == 0) stop("The baseline extent area is 0. Check shapefile topology.")
  
  target_area <- base_area * (1 + prop)
  
  # 4. Optimization function to solve for distance parameter 'd'
  optimization_function <- function(d) {
    temp_buffer <- aggregate(buffer(projected_union, width = d))
    temp_buffer <- makeValid(temp_buffer)
    return(expanse(temp_buffer) - target_area)
  }
  
  # 5. Root-finding optimization boundary execution
  buffer_distance <- uniroot(optimization_function, 
                             interval = c(0, 500000), 
                             extendInt = "yes", 
                             tol = 0.1)$root
  
  # 6. Build the final buffer vector layout and re-project to WGS84
  final_buffer <- aggregate(buffer(projected_union, width = buffer_distance))
  final_buffer <- makeValid(final_buffer)
  
  final_buffer_area   <- expanse(final_buffer)
  percentage_increase <- ((final_buffer_area / base_area) - 1) * 100
  
  return(list(
    buffer        = project(final_buffer, "EPSG:4326"), 
    increase_pct  = percentage_increase,
    distance_m    = buffer_distance,
    base_area_m2  = base_area,
    final_area_m2 = final_buffer_area
  ))
}

# Perform spatial intersection mapping against the dynamic target bounds
validate_spatial_bounds <- function(df_species, iucn_path, expert_path, prop = 0.10, metric_crs = "EPSG:5070") {
  
  # Retrieve optimized polygon parameters
  buffer_results <- calculate_dynamic_buffer(iucn_path, expert_path, prop, metric_crs)
  final_buffer   <- buffer_results$buffer 
  
  # Map dataframe entities to point vectors
  species_points <- vect(df_species, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  
  # Point-in-polygon spatial evaluation matrix check
  points_inside <- is.related(species_points, final_buffer, "intersects")
  
  df_species$.union <- as.logical(points_inside)
  return(df_species)
}

# --- PHASE III: Pipeline Execution ---

# Evaluates verification steps across a single selected data profile
process_species_pipeline <- function(df_post_cc, iucn_path, expert_path, prop = 0.10) {
  filtered_eval <- df_post_cc %>% filter(.coord_ok & .cc_ok)
  
  if (nrow(filtered_eval) == 0) {
    df_post_cc$.union      <- FALSE
    df_post_cc$.is_optimal <- FALSE
    return(df_post_cc)
  }
  
  validated_df <- validate_spatial_bounds(filtered_eval, iucn_path, expert_path, prop)
  
  df_post_cc %>%
    left_join(validated_df %>% select(gbifID, .union), by = "gbifID") %>%
    mutate(
      .union      = ifelse(is.na(.union), FALSE, .union),
      .is_optimal = .coord_ok & .year_ok & .sour_ok & .cc_ok & .union
    )
}

# Loop core calculation models through entire species dataset stack
curated_results <- map(
  seq_along(post_cc_list),
  ~ process_species_pipeline(post_cc_list[[.x]], iucn_paths[[.x]], expert_paths[[.x]])
)

# --- PHASE IV: Curated Dataset Exporting ---

# Verify and create output destination tracking folders
output_directory <- "curated_output/"
if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

# Batch export isolated outputs
purrr::walk(curated_results, ~ {
  species_name <- unique(.x$species_co)
  
  if(length(species_name) == 0 || is.na(species_name)) {
    species_name <- "Unknown_Species"
  }
  
  # Clean up string syntax values (converting spacing blocks to clean characters)
  clean_name  <- gsub("[[:space:]]|[/]", "_", species_name)
  output_file <- paste0(output_directory, clean_name, "_curated.csv")
  
  write.csv(.x, output_file, row.names = FALSE)
  
  # Print script progression update notice
  message(paste("Curated file successfully exported:", output_file))
})
