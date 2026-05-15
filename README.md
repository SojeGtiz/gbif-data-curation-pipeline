# Pipeline for data curation #

# Description
This repository contains an automated R pipeline designed to clean, filter, and validate species occurrence data using taxonomic, quality, and biogeographic criteria.

# Project Structure
To run the pipeline, your local directory should be organized as follows:

```text
├── data/
│   ├── raw_csv/       # Put your GBIF .csv files here
│   ├── iucn_shp/      # Put IUCN range .shp files here
│   └── expert_shp/    # Put Expert-validated .shp files here
├── output/            # Script will automatically save results here
├── main_pipeline.R    # The R script with the code
└── README.md          # This documentation file

# Data Sources 
The workflow relies on three main types of input files:
1. GBIF Occurrence Dataset (.csv)
  The files located in data/raw_csv/ correspond to occurrence records downloaded from the Global Biodiversity Information Facility (GBIF).
  Format: Darwin Core (DwC) standard.
  Key Fields Required: decimalLatitude, decimalLongitude, year, basisOfRecord, species_co, gbifID.

2. IUCN Range Polygons (.shp)
  The files in data/iucn_shp/ represent the species distribution ranges published by the IUCN Red List of Threatened Species.
  Source: Extracted from the IUCN main polygon dataset.

3. Expert-Validated Polygons (.shp)
  The files in data/expert_shp/ are custom shapefiles created by experts.
  Role: These polygons complement global datasets by adding specialized regional knowledge.

# Prerequisites
install.packages(c("sf", "terra", "dplyr", "purrr", "CoordinateCleaner", "tidyr", "ggplot2", "tidyterra", "rnaturalearth"))

# Running the Pipeline
  - Clone this repository or download the files.
  - Place your species data in the respective folders inside data/ (ensure the files match alphabetically by species name so the loop pairs them correctly).
  - Open main_pipeline.R and run the script.

# Expected Outputs
After a successful run, the following files will be generated in the output/ folder:
  - filtering_summary.csv: A table showing the percentage of clean vs. rejected data per species.
  - cleaned_occurrence_data.csv: The final, clean database ready for ecological modeling (e.g., MaxEnt).
  - quality_analysis.tiff: A high-resolution plot (300+ DPI) summarizing the validation status for publication.

# Author:
María José Morán-Gutiérrez, Centro de Investigaciones Biológicas del Noroeste, S.C.
Contact: mjmoran@pg.cibnor.mx
