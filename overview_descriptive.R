# intro -------

pacman::p_load(
  rio,        # importing data  
  here,       # relative file pathways  
  janitor,    # data cleaning and tables
  lubridate,  # working with dates
  matchmaker, # dictionary-based cleaning
  epikit,     # age_categories() function
  tidyverse,  # data management and visualization
  skimr,
  readxl,      # reads excel datasets
  units,
  kableExtra,        # pivoting
  circlize          # colorRamp2 heatmap
)


# IMPORT 
linelist_raw <- import("life_expectancy.csv")

linelist <- linelist_raw %>%
## name cleaning
janitor::clean_names() 

# Overview dataset => plausibility-check
## show
head(linelist)
## general
glimpse(linelist)
## structure
str(linelist)
## Variable names
names(linelist)

# descriptive statistics
## summary => plausibility check
summary(linelist)
## value cleaning
skimr::skim(linelist)
## additional counts

