# Pipeline for data curation #

# Description
This repository contains an automated R pipeline designed to clean, filter, and validate species occurrence data using taxonomic, quality, and biogeographic criteria.

# Project Structure
To run the pipeline, your local directory should be organized as follows:

```text
├── data/
│   ├── base/                      # Put your GBIF .csv files here
│   ├── pol_iucn/                  # Put IUCN range .shp files here
│   └── pol_exp/                   # Put Expert-validated .shp files here
├── output/                        # Script will automatically save results here
├── Filtering_pipeline_script.R    # The R script with the code
└── README.md                      # This documentation file
```

# Data Sources 
The workflow relies on three main types of input files:
1. **GBIF Occurrence Dataset (.csv)**
  The files located in data/base/ correspond to occurrence records downloaded from the Global Biodiversity Information Facility (GBIF).
  Format: Darwin Core (DwC) standard.
  Key Fields Required: decimalLatitude, decimalLongitude, year, basisOfRecord, gbifID, and species_co.
_Note: **species_co** is created post-taxonomic harmonization and can be replaced by the standard **species** column if necessary._

2. **IUCN Range Polygons (.shp)**
  The files in data/pol_iucn/ represent the species distribution ranges published by the IUCN Red List of Threatened Species.
  Source: Extracted from the IUCN main polygon dataset.

3. **Expert-Validated Polygons (.shp)**
  The files in data/pol_exp/ are custom shapefiles created by experts.
  Role: These polygons complement global datasets by adding specialized regional knowledge.

# Prerequisites
install.packages(c("sf", "terra", "dplyr", "purrr", "CoordinateCleaner", "tidyr", "scales", "readr"))

# Running the Pipeline
  - Clone this repository or download the files.
  - Place your species data in the respective folders inside data/.
  - **Ensure that files across all three folders match alphabetically by species name so the loop pairs them correctly.**
    _Example:_ Genus_species.csv, Genus_species.shp (IUCN), and Genus_species.shp (Expert).
  - Open Filtering_pipeline_script.R and run the script.

# Expected Outputs
After a successful run, the following files will be generated in the output/ folder:
  - sciuridae_database.csv: The final, clean database ready for ecological modeling (e.g., MaxEnt).

# Author:
María José Morán-Gutiérrez, Centro de Investigaciones Biológicas del Noroeste, S.C.
Contact: mjmoran@pg.cibnor.mx
Raúl Octavio Martínez-Rincón, Centro de Investigaciones Biológicas del Noroeste, S.C.
Contact: rrincon@cibnor.mx
