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
linelist_ari <- import("ari.csv")
linelist_Sc <- import("Sc.csv")
linelist_Y <- import("Y.csv")
linelist_Y_death <- import("Y_death.csv")



# Overview dataset => plausibility-check
## show
head
## general
glimpse(linelist_raw)
## structure
str(linelist_raw)
## Variable names
names(linelist_raw)

# descriptive statistics
## summary => plausibility check
summary(linelist_raw)
## skimr
skimr::skim(linelist_raw)
## additional counts
linelist_raw %>%
  count(sex)

linelist_raw %>%
  count(smoking_status)

linelist_raw %>%
  count(complication)


# Visualisation

## ggplot
### age
ggplot(linelist_raw, aes(x = age)) +
  geom_histogram(bins = 30, fill = "steelblue") +
  theme_minimal() +
  labs(title = "Distribution of Age")
### bmi
ggplot(linelist_raw, aes(x = bmi)) +
  geom_histogram(bins = 30, fill = "darkgreen") +
  theme_minimal() +
  labs(title = "Distribution of BMI")
