# ANALYSIS: SOCIOECONOMIC DETERMINANTS OF LIFE EXPECTANCY --------
# Research Question: What is the influence of socioeconomic factors and
#                    education on life expectancy?
# Note on how we developed this: We have made it such that you can put the cursor
# on line 18 and click on run continuously as it flows through all the way to the
# end of the script. Keep an eye on the comments and, console and output
# Sections:
#   1. Setup & Data Preparation
#   2. Descriptive Statistics
#   3. Correlation Analysis
#   4. Group Comparisons
#   5. Temporal Analysis
#   6. Scatter plots with subgroups
#   7. Multiple Linear Regression
#   8. Interaction / Moderation Analysis
# -----------------------------------------------------------------------------

# 1. SETUP & DATA PREPARATION ---------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  janitor,        # clean_names()
  tidyverse,      # dplyr, ggplot2, tidyr, purrr
  psych,          # describe()
  ggcorrplot,     # correlation heatmaps
  car,            # vif()
  lmtest,         # bptest(), dwtest()
  dunn.test,      # Dunn post-hoc
  emmeans,        # estimated marginal means
  scales,         # percent_format()
  viridis,        # viridis scales
  broom,          # tidy()
  nortest,        # ad.test() — kept for reference, not used in tables
  moments,        # skewness / kurtosis
  effectsize,     # cohens_d(), eta_squared(), cramers_v()
  parameters,     # model_parameters()
  gtsummary,      # publication Table 1
  gt,             # gt tables
  rstatix,        # tidy test wrappers
  ggpubr,         # stat_compare_means()
  gganimate,      # animated GIFs
  gifski,         # GIF renderer for gganimate
  transformr,      # smooth transitions in gganimate
  purrr,
  GGally,
  ppcor,
  ggrepel,
  pROC,
  caret
)
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate

## ── Global ggplot theme -------------------------------------------------------
theme_pub <- theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, colour = "grey40"),
    axis.title      = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text      = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    panel.grid.major = element_line(colour = "grey92"),
    panel.grid.minor = element_blank()
  )
theme_set(theme_pub)

palette_status <- c("Developed" = "#2166AC", "Developing" = "#D6604D")

## mode calculation formula -------
get_mode <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}



## ── 1.1  Load & clean ---------------------------------------------------------
df_raw <- read.csv("life_expectancy.csv") %>%
  clean_names() %>%
  rename(
    life_exp      = life_expectancy,
    adult_mort    = adult_mortality,
    pct_expend    = percentage_expenditure,
    tot_expend    = total_expenditure,
    under5_deaths = under_five_deaths,
    thinness_1_19 = thinness_1_19_years,
    thinness_5_9  = thinness_5_9_years,
    income_comp   = income_composition_of_resources
  ) %>%
  mutate(status = factor(status, levels = c("Developed", "Developing")))

cat("Raw dimensions:", nrow(df_raw), "rows x", ncol(df_raw), "columns\n")

## ── 1.2  Missing data summary -------------------------------------------------
missing_summary <- df_raw %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(100 * n_missing / nrow(df_raw), 1)) %>%
  filter(n_missing > 0) %>%
  arrange(desc(n_missing))

cat("\nMissing data (variables with any NA):\n")
print(missing_summary)

## ── 1.3  Drop rows missing outcome or key predictors; impute remainder --------
df_dropped <- df_raw %>%
  filter(!is.na(life_exp), !is.na(schooling),
         !is.na(income_comp), !is.na(gdp)) %>%
  mutate(across(where(is.numeric),
                ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))
# VARIATION 2: Grouped median imputation (Specific to Country)
df_med <- df_raw %>%
  # Group by country so mutations happen within each nation's timeline
  group_by(country) %>%
  mutate(across(where(is.numeric), function(x) {
    # Calculate the median for this specific country
    country_median <- median(x, na.rm = TRUE)
    
    # If the country has valid data, use its median. 
    # If it's missing entirely, leave it as NA for the next step.
    if (!is.na(country_median)) {
      ifelse(is.na(x), country_median, x)
    } else {
      x 
    }
  })) %>%
  ungroup() %>% # Always ungroup after group_by modifications
  
  # CRITICAL FALLBACK STEP: 
  # For countries completely missing GDP/Population across all 15 years,
  # impute the global median baseline for that specific year.
  group_by(year) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  ungroup()

## lets see the changes-------------
### 1. Create the shadow flag columns to map where original NAs are located
df_flags <- df_raw %>%
  # Use c() to combine selection rules and -any_of() to safely exclude "year"
  mutate(across(c(where(is.numeric), -any_of("year")), 
                ~ as.numeric(is.na(.)), 
                .names = "{.col}_imputed")) 

### 2. Perform the country-specific median imputation
df_imputed_highlighted <- df_flags %>%
  group_by(country) %>%
  # Exclude "year" and our newly created "_imputed" columns from the math
  mutate(across(c(where(is.numeric), -ends_with("_imputed"), -any_of("year")), function(x) {
    country_median <- median(x, na.rm = TRUE)
    if (!is.na(country_median)) {
      ifelse(is.na(x), country_median, x)
    } else {
      x 
    }
  })) %>%
  ungroup() %>%
  
  # Fallback step: Global median baseline for that specific year
  group_by(year) %>%
  mutate(across(c(where(is.numeric), -ends_with("_imputed"), -any_of("year")), 
                ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  ungroup()

### Example A: View rows where GDP was imputed to inspect the values
df_imputed_highlighted %>%
  filter(gdp_imputed == 1) %>%
  select(country, year, gdp, gdp_imputed) %>%
  head()

### Example B: Calculate exactly how many values were imputed per variable
imputation_summary <- df_imputed_highlighted %>%
  summarise(across(ends_with("_imputed"), sum)) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Total_Imputed_Rows")

print(imputation_summary)
cat("\nCleaned dimensions:", nrow(df_med), "rows x", ncol(df_med), "columns\n")
cat("Status distribution:\n"); print(table(df_med$status))

## ── 1.4  Derived variables ----------------------------------------------------
df <- df_med %>% #change to df_med for median imputed values
  mutate(
    log_gdp       = log1p(gdp),
    log_population = log1p(population),
    schooling_cat = cut(schooling,
                        breaks = quantile(schooling, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                        labels = c("Low", "Medium", "High"),
                        include.lowest = TRUE),
    income_cat    = cut(income_comp,
                        breaks = quantile(income_comp, c(0, 0.25, 0.75, 1), na.rm = TRUE),
                        labels = c("Low", "Middle", "High"),
                        include.lowest = TRUE),
    life_exp_cat  = cut(life_exp,
                        breaks = c(0, 60, 70, 80, Inf),
                        labels = c("<60", "60-70", "70-80", ">80"),
                        right  = FALSE),
    year_bracket  = cut(year,
                        breaks = c(1999, 2004, 2009, 2015),
                        labels = c("2000-2004", "2005-2009", "2010-2015"),
                        include.lowest = TRUE)
  )



# 2. DESCRIPTIVE STATISTICS ----------

# select variables
socio_edu_vars <- c("life_exp", "schooling", "income_comp", "gdp",
                    "log_gdp", "tot_expend", "adult_mort",
                    "infant_deaths", "hiv_aids", "bmi", "alcohol")

## mean, mode, median, range-------
target_years <- c(2000, 2005, 2010, 2015)
target_vars <- c("life_exp", "gdp", "bmi")

df_milestone_summary <- df %>%
  # 1. Filter for your specific years
  filter(year %in% target_years) %>%
  select(country, year, all_of(target_vars)) %>%
  
  # 2. Pivot to long format so we can group by year AND variable simultaneously
  pivot_longer(cols = all_of(target_vars), names_to = "variable", values_to = "value") %>%
  
  # 3. Drop NAs so our math and country extractions don't fail
  filter(!is.na(value)) %>%
  
  # 4. Group and calculate all statistics
  group_by(year, variable) %>%
  summarise(
    Mean   = round(mean(value), 1),
    SD     = round(sd(value), 1),
    Median = round(median(value), 1),
   # Mode = get_mode(value),
    Range  = paste0(round(min(value), 1), " to ", round(max(value), 1)),
    
    # THE TRICK: Extract the country name that matches the minimum and maximum values
    Lowest_Country  = country[which.min(value)],
    Highest_Country = country[which.max(value)],
    
    .groups = "drop"
  ) %>%
  
  # 5. Optional: Sort cleanly by variable then year
  arrange(variable, year)

# Print the final summary table
print(df_milestone_summary)

tbl_milestones <- df_milestone_summary %>%
  
  # The Magic Step: This automatically creates distinct sections for 2000, 2005, etc.
  gt(groupname_col = "year") %>%
  
  # 1. Add Titles and Headers
  tab_header(
    title = md("**Global Health & Economic Milestones**"),
    subtitle = "Summary statistics and extreme values across 15 years"
  ) %>%
  
  # 2. Clean up the column names for the audience
  cols_label(
    variable        = "Indicator",
    Lowest_Country  = "Lowest (Country)",
    Highest_Country = "Highest (Country)"
  ) %>%
  
  # 3. Format specific column strings (Capitalize the variable names)
  text_transform(
    locations = cells_body(columns = variable),
    fn = function(x) str_to_title(str_replace_all(x, "_", " "))
  ) %>%
  
  # 4. Alignment
  cols_align(align = "left", columns = c(variable, Lowest_Country, Highest_Country)) %>%
  cols_align(align = "center", columns = c(Mean, SD, Median, Range)) %>%
  
  # 5. Styling to make it pop on a slide
  tab_style(
    style = cell_text(weight = "bold", color = "#2C3E50"),
    locations = cells_row_groups() # Bolds the Year headers
  ) %>%
  tab_options(
    table.width = pct(100),
    heading.align = "left",
    column_labels.font.weight = "bold",
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    row_group.background.color = "gray95" # Adds a subtle background to the Year dividers
  )

# Render the table
tbl_milestones

## Sample characteristics (gtsummary) -------------------------------
# One row per country-year; for a cleaner Table 1, use one row per country
# (most recent year) so values are not repeated across years
df_t1 <- df %>%
  group_by(country) %>%
  slice_max(year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # SAFEGUARD: Force R to treat all your predictors as numeric
  # This prevents character/factor bugs if "NA" strings snuck into the data
  mutate(across(-c(country, status), as.numeric)) %>%
  select(status, life_exp, schooling, income_comp, gdp,
         tot_expend, adult_mort, infant_deaths, hiv_aids, bmi, alcohol)

tbl1 <- df_t1 %>%
  tbl_summary(
    by = status,
    
    # THE FIX: Explicitly force gtsummary to treat everything as continuous
    # (It automatically ignores 'status' since it is the 'by' grouping variable)
    type = list(everything() ~ "continuous"), 
    
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})"
      # Removed the categorical statistic since you have no categorical predictors here
    ),
    digits  = all_continuous() ~ 2,
    label   = list(
      life_exp      ~ "Life Expectancy (years)",
      schooling     ~ "Schooling (years)",
      income_comp   ~ "Income Composition Index",
      gdp       ~ "GDP per capita)",
      tot_expend    ~ "Health Expenditure (% GDP)",
      adult_mort    ~ "Adult Mortality (per 1,000)",
      infant_deaths ~ "Infant Deaths (per 1,000)",
      hiv_aids      ~ "HIV/AIDS Deaths (per 1,000)",
      bmi           ~ "Mean BMI",
      alcohol       ~ "Alcohol Consumption (L/capita)"
    ),
    missing = "no"
  ) %>%
  add_p(
    test = list(
      all_continuous()  ~ "wilcox.test"
      # Removed categorical test since everything is forced to continuous
      # opted to add the U test here as the rendered table visualizes better
      # tests for the variables run later in the code, results in the *!console
    ),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Sample characteristics by country development status** (most recent year per country)")

print(tbl1)

## PREP THE DATA FOR MULTI-YEAR ANALYSIS-----------------------
df_years <- df %>%
  # Select milestone years to keep the table readable
  filter(year %in% c(2000, 2005, 2010, 2015)) %>%
  
  # ERROR FIX: Strip NAs and force 'status' into a clean character string 
  # to destroy any phantom factor levels causing the 2-level error
  filter(!is.na(status) & status != "") %>%
  mutate(status = str_trim(as.character(status))) %>%
  
  # Force predictors to numeric to prevent the Chi-Square bug
  mutate(across(-c(country, status, year), as.numeric)) %>%
  
  # Select final columns
  select(year, status, life_exp, schooling, income_comp, gdp,
         tot_expend, adult_mort, infant_deaths, hiv_aids, bmi, alcohol)


## GENERATE THE STRATIFIED YEAR-BY-YEAR TABLE --*better------
tbl_yearly <- df_years %>%
  # tbl_strata splits the table into side-by-side columns based on the 'year'
  tbl_strata(
    strata = year,
    .tbl_fun = ~ .x %>%
      tbl_summary(
        by = status,
        type = list(everything() ~ "continuous"), 
        statistic = list(all_continuous() ~ "{mean} ({sd})"),
        digits  = all_continuous() ~ 2,
        label   = list(
          life_exp      ~ "Life Expectancy (years)",
          schooling     ~ "Schooling (years)",
          income_comp   ~ "Income Composition Index",
          gdp           ~ "GDP per capita",
          tot_expend    ~ "Health Expenditure (% GDP)",
          adult_mort    ~ "Adult Mortality (per 1,000)",
          infant_deaths ~ "Infant Deaths (per 1,000)",
          hiv_aids      ~ "HIV/AIDS Deaths (per 1,000)",
          bmi           ~ "Mean BMI",
          alcohol       ~ "Alcohol Consumption (L/capita)"
        ),
        missing = "no"
      ) %>%
      add_p(
        test = list(all_continuous() ~ "wilcox.test"),
        pvalue_fun = ~ style_pvalue(.x, digits = 3)
      ) %>%
      modify_header(label ~ "**Variable**")
  ) %>%
  # Add a clean caption to the overall stratified table
  modify_caption("**Longitudinal characteristics by country development status** (2000 vs 2007 vs 2015)")

# Print the table
# Investigated the reason we have 151 in developing countries for hours!
# Turns out there are 10 countries with data for only one year, 2013
print(tbl_yearly)

## ── Normality tests -- uses all observ--------------------------------
key_vars <- c("life_exp", "schooling", "income_comp", "gdp", "adult_mort", "log_gdp") #can use log_gdp

tbl2 <- purrr::map(key_vars, function(v) {
  
  x <- as.numeric(na.omit(df[[v]])) 
  
  # 1. THE FIX: Safety check for minimum sample size required by shapiro.test
  if (length(x) >= 3) {
    sw <- shapiro.test(sample(x, min(5000, length(x))))
    sw_w   <- round(sw$statistic, 4)
    sw_p   <- round(sw$p.value, 4)
    normal <- ifelse(sw$p.value > 0.05, "Yes", "No")
  } else {
    # If N < 3, gracefully return NAs instead of crashing
    sw_w   <- NA
    sw_p   <- NA
    normal <- "N/A"
  }
  
  tibble(
    Variable  = v,
    N         = length(x), 
    Mean      = round(mean(x, na.rm = TRUE), 2),
    SD        = round(sd(x, na.rm = TRUE), 2),
    Skewness  = round(moments::skewness(x, na.rm = TRUE), 3),
    Kurtosis  = round(moments::kurtosis(x, na.rm = TRUE), 3),
    SW_W      = sw_w,
    SW_p      = sw_p,
    Normal    = normal
  )
  
}) %>% 
  purrr::list_rbind() %>%   
  mutate(
    Variable = dplyr::recode_values(Variable,
                                    "life_exp"    ~ "Life Expectancy", 
                                    "schooling"   ~ "Schooling",
                                    "income_comp" ~ "Income Comp. Index",      
                                    "gdp"     ~ "GDP",
                                    "log_gdp" ~ "log GDP",
                                    "adult_mort"  ~ "Adult Mortality",         
                                    default = Variable
    )
  )
# skewness - indicating where majority of observations lie relative to the mean
# and how big the tail is. negative means more values to the right of the mean
# and tail on the left is long vice versa is true
# kurtosis measures how the data peaks, positive values mean a high peak, negative
# mean a flat peak
# both of these have a normal value of 0 (zero)
tbl2_gt <- tbl2 %>%
  gt() %>%
  tab_header(title = "Table 2. Normality Assessment of Key Continuous Variables") %>%
  tab_spanner(label = "Descriptive", columns = c(N, Mean, SD, Skewness, Kurtosis)) %>%
  tab_spanner(label = "Shapiro-Wilk", columns = c(SW_W, SW_p, Normal)) %>%
  cols_label(SW_W = "W", SW_p = "p-value", Normal = "Normal?") %>%
  fmt_number(columns = c(Mean, SD, Skewness, Kurtosis, SW_W), decimals = 3) %>%
  fmt_number(columns = SW_p, decimals = 4) %>%
  tab_style(
    style = cell_text(color = "red"),
    locations = cells_body(columns = Normal, rows = Normal == "No")
  ) %>%
  tab_footnote("Shapiro-Wilk W; p < 0.05 indicates departure from normality. Sample drawn (n=5000) for large datasets.")

# best to have a visual inspection, and perform the shapiro-wilk test 
# test is too sensitive with large datasets
# another identifiable distribution? eg lognormal; outliers in the data; 
# if the deviation is small we can ignore that
# in this case it wasnt small to be ignored so we chose non-parametric tests
# we still compare with the parametric tests anyway, as the decision is
# not supposed to be automated.

print(tbl2_gt)

## ── Figure 1: Distribution panel ----------------------------------------------
p_hist <- df %>%
  select(all_of(key_vars)) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = dplyr::recode(Variable,
    life_exp    = "Life Expectancy",
    schooling   = "Schooling (years)",
    income_comp = "Income Comp. Index",
    gdp     = "GDP per capita",
    log_gdp = "log(GDP per capita)",
    adult_mort  = "Adult Mortality"
  )) %>%
  ggplot(aes(Value)) +
  geom_histogram(bins = 35, fill = "#2166AC", colour = "white", alpha = 0.85) +
  facet_wrap(~ Variable, scales = "free", ncol = 3) +
  labs(title = "Figure 1. Distribution of Key Variables",
       x = NULL, y = "Frequency")
print(p_hist)

## ── Figure 2: Violin-boxplot by status ----------------------------------------
p_violin_status <- ggplot(df, aes(status, life_exp, fill = status)) +
  geom_violin(alpha = 0.35, colour = NA) +
  geom_boxplot(width = 0.22, outlier.shape = 21, outlier.size = 1.5,
               outlier.alpha = 0.5) +
  scale_fill_manual(values = palette_status) +
  labs(title = "Figure 2. Life Expectancy by Country Development Status",
       x = "Development Status", y = "Life Expectancy (years)") +
  guides(fill = "none")
print(p_violin_status)

## ── Figure 3: Countries by status (bar) ---------------------------------------
p_bar_status <- df %>%
  distinct(country, status) %>%
  count(status) %>%
  ggplot(aes(status, n, fill = status)) +
  geom_col(width = 0.5, colour = "white") +
  geom_text(aes(label = n), vjust = -0.4, fontface = "bold") +
  scale_fill_manual(values = palette_status, guide = "none") +
  labs(title = "Figure 3. Number of Countries by Development Status",
       x = "Status", y = "Count")
print(p_bar_status)


## MAP THE DATA (Filter years and fix country names)------------------
df_map_data <- df %>%
  filter(year %in% c(2000, 2005, 2010, 2015), !is.na(life_exp)) %>%
  select(country, year, life_exp) %>%
  # Translate WHO formal names to match the 'maps' package dictionary
  mutate(country = case_match(country,
                              "United States of America"                             ~ "USA",
                              "United Kingdom of Great Britain and Northern Ireland" ~ "UK",
                              "Russian Federation"                                   ~ "Russia",
                              "Republic of Korea"                                    ~ "South Korea",
                              "Democratic People's Republic of Korea"                ~ "North Korea",
                              "Iran (Islamic Republic of)"                           ~ "Iran",
                              "Bolivia (Plurinational State of)"                     ~ "Bolivia",
                              "Venezuela (Bolivarian Republic of)"                   ~ "Venezuela",
                              "Viet Nam"                                             ~ "Vietnam",
                              "Syrian Arab Republic"                                 ~ "Syria",
                              "United Republic of Tanzania"                          ~ "Tanzania",
                              "Côte d'Ivoire"                                        ~ "Ivory Coast",
                              "Lao People's Democratic Republic"                     ~ "Laos",
                              "Congo"                                                ~ "Republic of Congo",
                              .default = as.character(country)
  ))

# Extract the base world map
world_base <- map_data("world") %>%
  filter(region != "Antarctica") # Drop Antarctica to save space

# To facet polygons properly without dropping countries that lack data in a specific 
# year, we cross the base map with our target years.
world_expanded <- expand_grid(
  world_base,
  year = c(2000, 2005, 2010, 2015)
)

# Join our cleaned WHO data to the expanded map coordinates
map_final <- world_expanded %>%
  left_join(df_map_data, by = c("region" = "country", "year" = "year"))

# maps
map_2000 <- map_final %>% filter(year == 2000)

ggplot(map_2000, aes(x = long, y = lat, group = group, fill = life_exp)) +
  geom_polygon(color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "gray85", 
                       name = "Life Expectancy", limits = c(40, 85)) +
  theme_void(base_size = 14) + 
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))

map_2005 <- map_final %>% filter(year == 2005)

ggplot(map_2005, aes(x = long, y = lat, group = group, fill = life_exp)) +
  geom_polygon(color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "gray85", 
                       name = "Life Expectancy", limits = c(40, 85)) +
  theme_void(base_size = 14) + 
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))

map_2010 <- map_final %>% filter(year == 2010)

ggplot(map_2010, aes(x = long, y = lat, group = group, fill = life_exp)) +
  geom_polygon(color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "gray85", 
                       name = "Life Expectancy", limits = c(40, 85)) +
  theme_void(base_size = 14) + 
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))

map_2015 <- map_final %>% filter(year == 2015)

ggplot(map_2015, aes(x = long, y = lat, group = group, fill = life_exp)) +
  geom_polygon(color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "gray85", 
                       name = "Life Expectancy", limits = c(40, 85)) +
  theme_void(base_size = 14) + 
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))

# 3. CORRELATION ANALYSIS - Alena-------

corr_vars <- c("life_exp", "schooling", "income_comp", "log_gdp",
               "tot_expend", "adult_mort", "infant_deaths",
               "hiv_aids", "bmi", "alcohol", "hepatitis_b", "polio", "diphtheria")

corr_df      <- df %>% select(any_of(corr_vars)) %>% drop_na()
pearson_mat  <- cor(corr_df, method = "pearson",  use = "complete.obs")
spearman_mat <- cor(corr_df, method = "spearman", use = "complete.obs")

## scatter plots--------
# Aggregate all variables to the 193 country baselines
target_vars <- c("life_exp", "gdp", "schooling", "bmi", "income_comp", "adult_mort")

df_corr_prep <- df %>%
  filter(!is.na(status) & status != "") %>%
  mutate(status = str_trim(as.character(status))) %>%
  # Get the historical baseline for every country across ALL key variables
  group_by(country, status) %>%
  summarise(across(all_of(target_vars), ~mean(., na.rm = TRUE)), .groups = "drop") %>%
  # Pivot to a long format so we can facet the predictors
  pivot_longer(
    cols = -c(country, status, life_exp),
    names_to = "predictor",
    values_to = "value"
  ) %>%
  # Clean names for the presentation panel headers
  mutate(
    clean_predictor = case_match(predictor,
                                 "gdp"         ~ "GDP per Capita",
                                 "schooling"   ~ "Schooling (Years)",
                                 "bmi"         ~ "Mean BMI",
                                 "income_comp" ~ "Income Composition",
                                 "adult_mort"  ~ "Adult Mortality",
                                 .default = predictor
    )
  ) %>%
  # Drop NAs so the plots render cleanly
  filter(!is.na(value), !is.na(life_exp))

# GENERATE THE SCATTER PLOT PANEL
p_scatter <- ggplot(df_corr_prep, aes(x = value, y = life_exp)) +
  
  # Add the scatter points, colored by development status
  geom_point(aes(color = status), alpha = 0.6, size = 2) +
  
  # Add a smoothing line to detect non-linear curves (critical for correlation prep)
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1) +
  
  # Facet the grid, allowing the X-axis to adapt to each variable's unique scale
  facet_wrap(~ clean_predictor, scales = "free_x", ncol = 3) +
  
  scale_color_manual(values = c("Developed" = "#2C3E50", "Developing" = "#E74C3C")) +
  
  labs(
    title = "Predictor Relationships with Baseline Life Expectancy",
    subtitle = "Visualizing linearity, clusters, and outliers prior to correlation analysis (N = 193)",
    x = "Predictor Average (2000 - 2015)",
    y = "Mean Life Expectancy (Years)",
    color = "Status"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray90", color = NA),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Render the plot
print(p_scatter)

# data heavily skewed, comparing both correlations
# schooling
cor.test(df$life_exp, df$schooling, method = "pearson")
cor.test(df$life_exp, df$schooling, method = "spearman")
# adult mortality
cor.test(df$life_exp, df$adult_mort, method = "pearson")
cor.test(df$life_exp, df$adult_mort, method = "spearman")
# income comp
cor.test(df$life_exp, df$income_comp, method = "pearson")
cor.test(df$life_exp, df$income_comp, method = "spearman")
# gdp
cor.test(df$life_exp, df$gdp, method = "pearson")
cor.test(df$life_exp, df$gdp, method = "spearman")
# bmi
cor.test(df$life_exp, df$bmi, method = "pearson")
cor.test(df$life_exp, df$bmi, method = "spearman")

## ── Table 3: Pearson r and Spearman rho with life expectancy -----------------
tbl3 <- map_dfr(setdiff(names(corr_df), "life_exp"), function(v) {
  pt <- cor.test(corr_df$life_exp, corr_df[[v]], method = "pearson")
  st <- cor.test(corr_df$life_exp, corr_df[[v]], method = "spearman", exact = FALSE)
  tibble(
    Variable     = v,
    r_pearson    = round(pt$estimate, 3),
    p_pearson    = round(pt$p.value, 4),
    rho_spearman = round(st$estimate, 3),
    p_spearman   = round(st$p.value, 4)
  )
}) %>%
  mutate(
    sig_pearson  = case_when(p_pearson  < 0.001 ~ "***", p_pearson  < 0.01 ~ "**",
                             p_pearson  < 0.05  ~ "*",   TRUE              ~ ""),
    sig_spearman = case_when(p_spearman < 0.001 ~ "***", p_spearman < 0.01 ~ "**",
                             p_spearman < 0.05  ~ "*",   TRUE              ~ "")
  ) %>%
  arrange(desc(abs(r_pearson))) %>%
  mutate(Variable = dplyr::recode(Variable,
    schooling     = "Schooling",       income_comp   = "Income Comp. Index",
    log_gdp       = "log(GDP)",        tot_expend    = "Health Expenditure",
    adult_mort    = "Adult Mortality", infant_deaths = "Infant Deaths",
    hiv_aids      = "HIV/AIDS",        bmi           = "BMI",
    alcohol       = "Alcohol",         hepatitis_b   = "Hepatitis B",
    polio         = "Polio",           diphtheria    = "Diphtheria"
  ))

tbl3_gt <- tbl3 %>%
  select(Variable, r_pearson, sig_pearson, p_pearson, rho_spearman, sig_spearman, p_spearman) %>%
  gt() %>%
  tab_header(title = "Table 3. Correlations with Life Expectancy") %>%
  tab_spanner(label = "Pearson",  columns = c(r_pearson, sig_pearson, p_pearson)) %>%
  tab_spanner(label = "Spearman", columns = c(rho_spearman, sig_spearman, p_spearman)) %>%
  cols_label(r_pearson = "r", sig_pearson = "", p_pearson = "p",
             rho_spearman = "rho", sig_spearman = "", p_spearman = "p") %>%
  fmt_number(columns = c(r_pearson, rho_spearman), decimals = 3) %>%
  fmt_number(columns = c(p_pearson, p_spearman), decimals = 4) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = r_pearson, rows = abs(r_pearson) >= 0.5)
  ) %>%
  tab_footnote("*** p<0.001, ** p<0.01, * p<0.05. Sorted by |r|.")
print(tbl3_gt)

## ── 3B.3 Pearson correlation coefficients with 95% confidence intervals --------

corr_ci_tbl <- map_dfr(setdiff(names(corr_df), "life_exp"), function(v) {
  
  test <- cor.test(corr_df$life_exp, corr_df[[v]], method = "pearson")
  
  tibble(
    Variable = v,
    r = as.numeric(test$estimate),
    conf_low = test$conf.int[1],
    conf_high = test$conf.int[2],
    p_value = test$p.value
  )
}) %>%
  mutate(
    Variable = dplyr::recode(Variable,
                             schooling     = "Schooling",
                             income_comp   = "Income Comp. Index",
                             log_gdp       = "log(GDP)",
                             tot_expend    = "Health Expenditure",
                             adult_mort    = "Adult Mortality",
                             infant_deaths = "Infant Deaths",
                             hiv_aids      = "HIV/AIDS",
                             bmi           = "BMI",
                             alcohol       = "Alcohol",
                             hepatitis_b   = "Hepatitis B",
                             polio         = "Polio",
                             diphtheria    = "Diphtheria"
    ),
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  ) %>%
  arrange(r)

corr_ci_tbl

p_corr_forest <- ggplot(corr_ci_tbl,
                        aes(x = r, y = reorder(Variable, r))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high),
                 height = 0.2, colour = "grey40") +
  geom_point(aes(colour = r > 0), size = 3) +
  geom_text(
    aes(label = paste0(sprintf("%.2f", r), " ", sig)),
    hjust = ifelse(corr_ci_tbl$r > 0, -0.15, 1.15),
    size = 3
  ) +
  scale_colour_manual(
    values = c("TRUE" = "#2166AC", "FALSE" = "#D6604D"),
    guide = "none"
  ) +
  xlim(-1.1, 1.1) +
  labs(
    title = "Pearson Correlations with Life Expectancy",
    subtitle = "Points = Pearson r; lines = 95% confidence intervals",
    x = "Pearson r",
    y = NULL
  ) +
  theme_pub

print(p_corr_forest)


## ── 3B.4 Pearson vs Spearman comparison ---------------------------------------

pearson_spearman_comparison <- tbl3 %>%
  mutate(
    difference = abs(r_pearson - rho_spearman)
  ) %>%
  arrange(desc(difference))

pearson_spearman_comparison

p_pearson_spearman <- ggplot(
  pearson_spearman_comparison,
  aes(x = r_pearson, y = rho_spearman, label = Variable)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, colour = "#2166AC") +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
    title = "Pearson vs Spearman Correlations",
    subtitle = "Large deviations from the dashed line suggest non-linearity or outlier sensitivity",
    x = "Pearson r",
    y = "Spearman rho"
  ) +
  theme_pub

print(p_pearson_spearman)



## ── Figure 4: Pearson correlation heatmap -------------------------------------
p_corr <- ggcorrplot(pearson_mat,
                     method   = "square", type = "lower",
                     lab      = TRUE, lab_size = 2.5,
                     colors   = c("#D6604D", "white", "#2166AC"),
                     title    = "Figure 4. Pearson Correlation Matrix",
                     ggtheme  = theme_pub)
print(p_corr)

## ── Figure 5: Correlation bar chart with life expectancy ----------------------
p_corr_bar <- tbl3 %>%
  ggplot(aes(
    x = r_pearson,                    # 1. Map continuous values directly to X
    y = reorder(Variable, r_pearson), # 2. Map discrete variables directly to Y
    fill = r_pearson > 0
  )) +
  geom_col(colour = "white") +
  geom_text(
    aes(
      label = paste0(sprintf("%.3f", r_pearson), " ", sig_pearson),
      hjust = ifelse(r_pearson > 0, -0.1, 1.1) # 3. Moved inside aes() and removed tbl3$
    ), 
    size = 3
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#2166AC", "FALSE" = "#D6604D"),
    guide = "none"
  ) +
  labs(
    title = "Figure 5. Pearson r with Life Expectancy",
    y = NULL,         # Y is now your variable names, so we remove the label here
    x = "Pearson r"   # X is now your correlation value
  ) +
  xlim(c(-1.1, 1.1))  
print(p_corr_bar)


## ── Pairplot / scatterplot matrix ----------------------------------------

pairplot_df <- df %>%
  select(life_exp, schooling, income_comp, log_gdp, adult_mort, hiv_aids, bmi) %>%
  drop_na()

GGally::ggpairs(
  pairplot_df,
  title = "Pairplot of Key Predictors and Life Expectancy"
)



# 4. GROUP COMPARISONS -----------

## ── Table 4: Developed vs Developing binary comparison -----------------------
# Primary: Mann-Whitney (non-normal); Welch t and Cohen's d also reported
t_res  <- t.test(life_exp ~ status, data = df, var.equal = FALSE)
mwu_res <- wilcox.test(life_exp ~ status, data = df, conf.int = TRUE)
cd_res  <- cohens_d(life_exp ~ status, data = df)

compare_vars <- c("life_exp", "schooling", "income_comp", "log_gdp", "gdp",
                  "adult_mort", "hiv_aids", "tot_expend", "bmi", "thinness_1_19")

df_cal <- df %>%
  mutate(across(all_of(compare_vars), as.numeric))

### Improved table ---------
# 1. Isolate the target data and strip haven labels immediately
long_df <- df_cal %>%
  select(status, all_of(compare_vars)) %>%
  filter(status %in% c("Developed", "Developing")) %>%
  pivot_longer(cols = -status, names_to = "Variable", values_to = "value") %>%
  drop_na(value) %>%
  mutate(
    value = as.numeric(value),     # This safely destroys haven labels globally
    status = as.character(status)
  )

# 2. Let rstatix handle the calculations using group_by
stats_tests <- long_df %>%
  group_by(Variable) %>%
  wilcox_test(value ~ status) %>%
  select(Variable, W = statistic, p_value = p)

cohens_eff <- long_df %>%
  group_by(Variable) %>%
  cohens_d(value ~ status) %>%
  select(Variable, Cohens_d = effsize)

# 3. Calculate Medians & IQRs, then join everything together
tbl4 <- long_df %>%
  group_by(Variable, status) %>%
  summarise(
    Med = round(median(value), 2),
    IQR = paste0("(", round(quantile(value, .25), 1), "-", round(quantile(value, .75), 1), ")"),
    .groups = "drop"
  ) %>%
  # Pivot wider to get Developed and Developing as columns
  pivot_wider(names_from = status, values_from = c(Med, IQR)) %>%
  rename(
    Dev_Median   = Med_Developed,
    Dev_IQR      = IQR_Developed,
    Devng_Median = Med_Developing,
    Devng_IQR    = IQR_Developing
  ) %>%
  
  # 4. Join our tests and format
  left_join(stats_tests, by = "Variable") %>%
  left_join(cohens_eff, by = "Variable") %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***", 
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",  
      TRUE            ~ ""
    ),
    Variable = recode_values(Variable,
                             "life_exp"    ~ "Life Expectancy (years)", 
                             "schooling"   ~ "Schooling (years)",
                             "income_comp" ~ "Income Comp. Index",      
                             "gdp"     ~ "GDP per capita",
                             "log_gdp" ~ "log(GDP per capita)",
                             "adult_mort"  ~ "Adult Mortality",         
                             "hiv_aids"    ~ "HIV/AIDS Deaths",
                             "tot_expend"  ~ "Health Expenditure (% GDP)",
                             "bmi" ~ "BMI",
                             "thinness_1_19" ~ "Prevalence of thinness (ages 10-19)",
                             default = Variable
    )
  ) %>%
  select(Variable, Dev_Median, Dev_IQR, Devng_Median, Devng_IQR, W, p_value, Cohens_d, sig)

tbl4_gt <- tbl4 %>%
  select(Variable, Dev_Median, Dev_IQR, Devng_Median, Devng_IQR, W, p_value, sig, Cohens_d) %>%
  gt() %>%
  tab_header(title = "Table 4. Comparison of Key Variables by Development Status") %>%
  tab_spanner(label = "Developed",   columns = c(Dev_Median, Dev_IQR)) %>%
  tab_spanner(label = "Developing",  columns = c(Devng_Median, Devng_IQR)) %>%
  tab_spanner(label = "Test",        columns = c(W, p_value, sig)) %>%
  cols_label(Dev_Median = "Median", Dev_IQR = "IQR",
             Devng_Median = "Median", Devng_IQR = "IQR",
             W = "W", p_value = "p", sig = "",
             Cohens_d = "Cohen's d") %>%
  fmt_number(columns = c(Dev_Median, Devng_Median, Cohens_d), decimals = 2) %>%
  fmt_number(columns = p_value, decimals = 4) %>%
  tab_footnote("Mann-Whitney U test. *** p<0.001, ** p<0.01, * p<0.05. Cohen's d from Welch t-test.")
print(tbl4_gt)

## One 16-year average per country analysed ---------------------------------
df_country_avg <- df_med %>%
  # Remove rows with missing status or life expectancy
  filter(!is.na(status) & status != "", !is.na(life_exp)) %>%
  mutate(status = str_trim(as.character(status))) %>%
  
  # Group by country and status, then calculate the mean for each
  group_by(country, status) %>%
  summarise(
    mean_life_exp = mean(life_exp, na.rm = TRUE),
    .groups = "drop"
  )

# Verify the dimensions (This should ideally output close to 193 rows)
cat("Number of unique countries in new dataset:", nrow(df_country_avg), "\n\n")

# RUN THE STANDARD T-TEST (Parametric)
cat("--- Standard Independent T-Test (Welch's) ---\n")
test_normal <- t.test(mean_life_exp ~ status, data = df_country_avg)
print(test_normal)


### RUN THE WILCOXON RANK-SUM TEST (Non-Parametric) -------------------------
cat("\n--- Wilcoxon Rank-Sum Test (Mann-Whitney U) ---\n")
test_non_normal <- wilcox.test(mean_life_exp ~ status, data = df_country_avg)
print(test_non_normal)



# THE VISUALIZATION: Box Plot with Overlaid Data Points
p_test_box <- ggplot(df_country_avg, aes(x = status, y = mean_life_exp, fill = status)) +
  
  # Draw the boxplots. We hide the default outliers because the jitter layer handles them.
  geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
  
  # THE UPGRADE: Plot the actual 193 countries as dots over the boxes
  #geom_jitter(width = 0.15, alpha = 0.6, size = 2, color = "gray20") + # confusing for no good reason
  
  scale_fill_manual(values = c("Developed" = "#2C3E50", "Developing" = "#E74C3C")) +
  
  labs(
    title = "Baseline Life Expectancy by Development Status",
    subtitle = "Each point represents a country's 16-year historical average (N = 193)",
    x = NULL, # X-axis label removed since the categories are self-explanatory
    y = "Mean Life Expectancy (Years)"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Render the plot
print(p_test_box)

### THE STATISTICAL RESULTS TABLE ---------------------------------
# Run the tests silently to extract their mathematical values
t_res <- t.test(mean_life_exp ~ status, data = df_country_avg)
w_res <- wilcox.test(mean_life_exp ~ status, data = df_country_avg)

# Build a clean dataframe of the results
df_stats <- tibble(
  `Statistical Test` = c("Welch Two Sample t-test", "Wilcoxon Rank-Sum Test"),
  `Statistic` = c(paste0("t = ", round(t_res$statistic, 2)), 
                  paste0("W = ", w_res$statistic)),
  `Degrees of Freedom` = c(round(t_res$parameter, 1), "—"),
  `P-Value` = c("< 0.001", "< 0.001")
)

# Format the dataframe into a stunning presentation table using gt
tbl_stats <- df_stats %>%
  gt() %>%
  tab_header(
    title = md("**Hypothesis Testing Results**"),
    subtitle = "Comparison of historical global life expectancy baselines"
  ) %>%
  # Center align the math, left align the text
  cols_align(align = "center", columns = c(`Statistic`, `Degrees of Freedom`, `P-Value`)) %>%
  cols_align(align = "left", columns = `Statistical Test`) %>%
  # Add some styling to make it pop on a slide
  tab_options(
    table.width = pct(90),
    heading.title.font.weight = "bold",
    column_labels.font.weight = "bold",
    table.border.top.color = "black",
    table.border.bottom.color = "black"
  )

# Render the table
tbl_stats


### EXTRACT DESCRIPTIVE STATISTICS (Mean & SD) -------------------------------
df_country_summary <- df_country_avg %>%
  group_by(status) %>%
  summarise(
    `Total Countries (N)` = n(),
    `Mean Life Expectancy` = mean(mean_life_exp, na.rm = TRUE),
    `Standard Deviation` = sd(mean_life_exp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Round to 1 decimal place so it looks perfect on a slide
  mutate(across(where(is.numeric), ~round(., 1)))

# Print cleanly to the console
print(df_country_summary)

## ── Table 5: Categorical associations ----------------------------------------
# we can make some categorical variables for this assessment
chi_tbl    <- table(df_cal$status, df_cal$schooling_cat)
chi_res    <- chisq.test(chi_tbl)
cv_school  <- effectsize::cramers_v(chi_tbl)


tbl5 <- tibble(
  Test       = c("Chi-square"),
  Comparison = c("Status x Schooling Category"),
  Statistic  = c(round(chi_res$statistic, 3)),
  df         = c(chi_res$parameter),
  p_value    = c(round(chi_res$p.value, 4)),
  Effect     = c(paste0("Cramer's V = ", round(cv_school$Cramers_v, 3)))
)

tbl5_gt <- tbl5 %>%
  gt() %>%
  tab_header(title = "Table 5. Categorical Association Tests") %>%
  fmt_number(columns = c(Statistic), decimals = 3) %>%
  fmt_number(columns = p_value, decimals = 4) %>%
  sub_missing(missing_text = "—") %>%
  tab_footnote("Schooling categories low, medium and high based on number of years <9, 10-17 >18.")
print(tbl5_gt)

## ── 4.1  One-way KW + Dunn: Life Expectancy ~ Schooling category --------------
cat("\n-- Kruskal-Wallis: Life Expectancy ~ Schooling Category --\n")
kw_school <- kruskal.test(life_exp ~ schooling_cat, data = df_cal)
print(kw_school)
cat("\nDunn post-hoc (Bonferroni):\n")
dunn_school <- dunn.test(df_cal$life_exp, df_cal$schooling_cat, method = "bonferroni",
                         alph = 0.05, list = FALSE)

cat("\n-- Kruskal-Wallis: Life Expectancy ~ Income Category --\n")
kw_income <- kruskal.test(life_exp ~ income_cat, data = df_cal)
print(kw_income)
cat("\nDunn post-hoc (Bonferroni):\n")
dunn_income <- dunn.test(df_cal$life_exp, df_cal$income_cat, method = "bonferroni",
                         alph = 0.05, list = FALSE)

## ── 4.2  Two-way ANOVA: Status x Schooling ------------------------------------
cat("\n-- Two-way ANOVA: Life Expectancy ~ Status x Schooling Category --\n")
anova_2way  <- aov(life_exp ~ status * schooling_cat, data = df_cal)
print(summary(anova_2way))
eta_2way <- effectsize::eta_squared(anova_2way)
eta_2way
cat("\nPartial eta-squared:\n"); print(eta_2way)

emm_2way <- emmeans(anova_2way, pairwise ~ status | schooling_cat,
                    adjust = "tukey")
cat("\nEmmeans pairwise (Status within Schooling):\n")
print(emm_2way$contrasts)

##── Figure 7: Boxplot by schooling category -----------------------------------
p_box_school <- ggplot(df_cal, aes(schooling_cat, life_exp, fill = schooling_cat)) +
  geom_violin(alpha = 0.3, colour = NA) +
  geom_boxplot(width = 0.28, outlier.shape = 21, outlier.size = 1.2,
               outlier.alpha = 0.4) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Figure 6. Life Expectancy by Schooling Level",
       x = "Schooling Tertile", y = "Life Expectancy (years)")
print(p_box_school)

## ── Figure 8: Boxplot by income category --------------------------------------
p_box_income <- ggplot(df_cal, aes(income_cat, life_exp, fill = income_cat)) +
  geom_violin(alpha = 0.3, colour = NA) +
  geom_boxplot(width = 0.28, outlier.shape = 21, outlier.size = 1.2,
               outlier.alpha = 0.4) +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  labs(title = "Figure 7. Life Expectancy by Income Composition Category",
       x = "Income Tertile", y = "Life Expectancy (years)")
print(p_box_income)

## ── Figure 9: Two-way interaction plot ----------------------------------------
p_twoway <- ggplot(df_cal, aes(schooling_cat, life_exp, fill = status, colour = status)) +
  geom_boxplot(alpha = 0.45, position = position_dodge(0.8),
               outlier.shape = 21, outlier.size = 1.2, outlier.alpha = 0.4) +
  scale_fill_manual(values = palette_status) +
  scale_colour_manual(values = palette_status) +
  labs(title    = "Figure 8. Life Expectancy by Schooling x Development Status",
       x        = "Schooling Tertile",
       y        = "Life Expectancy (years)",
       fill     = "Status", colour = "Status")
print(p_twoway)

## ── Figure 10: Mosaic-style stacked bar ----------------------------------------
p_mosaic <- df_cal %>%
  count(status, life_exp_cat) %>%
  group_by(status) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(status, pct, fill = life_exp_cat)) +
  geom_col(colour = "white", width = 0.6) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_viridis_d(option = "C") +
  labs(title = "Figure 9. Life Expectancy Category Distribution by Status",
       x = "Status", y = "Proportion", fill = "Life Exp. Category")
print(p_mosaic)

## ── 4.3  Repeated-measures ANOVA (2000, 2005, 2010, 2015) --------------------
cat("\n-- Repeated-Measures ANOVA: Life Expectancy over Time --\n")
rm_data <- df_cal %>%
  filter(year %in% c(2000, 2005, 2010, 2015)) %>%
  group_by(country) %>%
  filter(n() == 4) %>%
  ungroup() %>%
  mutate(year_f = factor(year), country = factor(country))

cat("Countries with all 4 time points:", length(unique(rm_data$country)), "\n")
rm_model <- aov(life_exp ~ year_f + Error(country / year_f), data = rm_data)
print(summary(rm_model))



# 5. TEMPORAL ANALYSIS ---------


yearly_summary <- df_cal %>%
  group_by(year, status) %>%
  summarise(
    avg_life_exp   = mean(life_exp, na.rm = TRUE),
    avg_schooling  = mean(schooling, na.rm = TRUE),
    avg_income     = mean(income_comp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

## ── Figure 11: Life expectancy gap over time ----------------------------------
p_yearly_trend <- ggplot(yearly_summary,
                         aes(year, avg_life_exp, colour = status, group = status)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_colour_manual(values = palette_status) +
  scale_x_continuous(breaks = seq(2000, 2015, 5)) +
  labs(title    = "Figure 10. Average Life Expectancy Over Time by Status",
       x        = "Year", y = "Mean Life Expectancy (years)",
       colour   = "Status")
print(p_yearly_trend)

## ── Table 6: Year-by-year Mann-Whitney gap analysis ---------------------------
test_vars <- c("life_exp", "schooling", "polio")

run_yearly_tests <- function(data, variable) {
  data %>%
    group_by(year) %>%
    summarise(
      Mean_Developed  = round(mean(.data[[variable]][status == "Developed"], na.rm = TRUE), 2),
      Mean_Developing = round(mean(.data[[variable]][status == "Developing"], na.rm = TRUE), 2),
      Gap               = round(abs(Mean_Developed - Mean_Developing), 2),
      p_value           = wilcox.test(.data[[variable]] ~ status, exact = FALSE)$p.value,
      .groups           = "drop"
    ) %>%
    mutate(
      Variable    = variable,
      p_value     = signif(p_value, 3),
      Significant = ifelse(p_value < 0.05, "Yes", "No")
    ) %>%
    select(Variable, year, Mean_Developed, Mean_Developing, Gap, p_value, Significant)
}

yearly_hypothesis_results <- bind_rows(lapply(test_vars, function(v) run_yearly_tests(df, v)))

tbl6_gt <- yearly_hypothesis_results %>%
  mutate(Variable = dplyr::recode(Variable,
    life_exp      = "Life Expectancy",
    schooling = "Schooling",
    polio         = "Polio"
  )) %>%
  gt(groupname_col = "Variable") %>%
  tab_header(title = "Table 6. Mann-Whitney U Tests by Year: Developed vs Developing") %>%
  cols_label(year = "Year", Mean_Developed = "Mean (Developed)",
             Mean_Developing = "Mean (Developing)", Gap = "Absolute Gap",
             p_value = "p-value", Significant = "Sig.") %>%
  fmt_number(columns = c(Mean_Developed, Mean_Developing, Gap), decimals = 2) %>%
  fmt_number(columns = p_value, decimals = 4) %>%
  tab_style(
    style = cell_text(color = "red", weight = "bold"),
    locations = cells_body(columns = Significant, rows = Significant == "Yes")
  ) %>%
  tab_footnote("Wilcoxon rank-sum test per year. Significant = p < 0.05.")
print(tbl6_gt)

## ── Figure 12: Gap significance over time -------------------------------------
p_yearly_sig <- ggplot(yearly_hypothesis_results, aes(year, Gap)) +
  geom_line(colour = "grey50", linewidth = 0.9) +
  geom_point(aes(colour = Significant, size = Significant)) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1,
             labeller = as_labeller(c(life_exp = "Life Expectancy",
                                      infant_deaths = "Infant Deaths",
                                      polio = "Polio"))) +
  scale_colour_manual(values = c("Yes" = "#D6604D", "No" = "#2166AC")) +
  scale_size_manual(values  = c("Yes" = 3, "No" = 2)) +
  scale_x_continuous(breaks = seq(2000, 2015, 5)) +
  labs(title    = "Figure 11. Gap Between Developed and Developing Countries Over Time",
       subtitle = "Red = statistically significant gap (p < 0.05)",
       x = "Year", y = "Absolute Difference in Medians",
       colour = "Significant", size = "Significant")
print(p_yearly_sig)



# 6. SCATTERPLOTS WITH YEAR SUBGROUP --------

# Each predictor (Schooling, Income Composition, log-GDP) is shown in three
# formats: (A) colour gradient by year, (B) faceted by year bracket,
# (C) animated GIF by year.


scatter_data <- df_cal %>%
  mutate(year_f = factor(year))

scatter_spec <- list(
  list(var = "schooling",   xlab = "Years of Schooling"),
  list(var = "income_comp", xlab = "Income Composition Index"),
  list(var = "log_gdp",     xlab = "log(GDP per capita)")
)

fig_counter <- 12   # continue figure numbering

for (spec in scatter_spec) {

  v    <- spec$var
  xlab <- spec$xlab

  # ── (A) Colour gradient by year -------------------------------------------
  p_grad <- ggplot(scatter_data,
                   aes(x = .data[[v]], y = life_exp,
                       colour = year, shape = status)) +
    geom_point(alpha = 0.45, size = 1.6) +
    geom_smooth(aes(group = status, colour = NULL, fill = status),
                method = "lm", se = TRUE, alpha = 0.15, colour = "grey30",
                linewidth = 0.8) +
    scale_colour_viridis_c(option = "C", name = "Year") +
    scale_fill_manual(values = palette_status) +
    scale_shape_manual(values = c("Developed" = 16, "Developing" = 17)) +
    labs(title    = paste0("Figure ", fig_counter, "A. ", xlab, " vs Life Expectancy (year gradient)"),
         x        = xlab, y = "Life Expectancy (years)",
         shape    = "Status", fill = "Status")
  print(p_grad)

  # ── (B) Faceted by year bracket --------------------------------------------
  p_facet <- ggplot(scatter_data %>% filter(!is.na(year_bracket)),
                    aes(x = .data[[v]], y = life_exp, colour = status)) +
    geom_point(alpha = 0.35, size = 1.3) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.9) +
    facet_wrap(~ year_bracket, ncol = 3) +
    scale_colour_manual(values = palette_status) +
    labs(title    = paste0("Figure ", fig_counter, "B. ", xlab, " vs Life Expectancy by Year Bracket"),
         x        = xlab, y = "Life Expectancy (years)",
         colour   = "Status")
  print(p_facet)

  # ── (C) Animated GIF -------------------------------------------------------
  p_anim <- ggplot(scatter_data,
                   aes(x = .data[[v]], y = life_exp,
                       colour = status, group = status)) +
    geom_point(alpha = 0.55, size = 2) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) +
    scale_colour_manual(values = palette_status) +
    labs(title    = paste0(xlab, " vs Life Expectancy  |  Year: {closest_state}"),
         subtitle = "Colour = development status",
         x        = xlab, y = "Life Expectancy (years)",
         colour   = "Status") +
    transition_states(year, transition_length = 2, state_length = 1) +
    ease_aes("cubic-in-out")

  gif_file <- paste0("fig_", fig_counter, "C_anim_", v, ".gif")
  anim_save(gif_file,
            animation = animate(p_anim, nframes = 80, fps = 10,
                                width = 700, height = 480,
                                renderer = gifski_renderer()))
  cat("Saved:", gif_file, "\n")

  fig_counter <- fig_counter + 1
}

# fig_counter is now 15 after the loop (12, 13, 14 used)



# 7. MULTIPLE LINEAR REGRESSION ---------


# 8. Logistic regression ---------
## ── 9.1 Create binary outcome -------------------------------------------------

df_logit <- df %>%
  mutate(
    high_life_exp = ifelse(life_exp > 70, 1, 0),
    high_life_exp = factor(high_life_exp, levels = c(0, 1),
                           labels = c("≤70 years", ">70 years"))
  )

table(df_logit$high_life_exp)

## ── 9.2 Fit logistic regression model ----------------------------------------

log_model <- glm(
  high_life_exp ~ schooling + income_comp + log_gdp +
    adult_mort + hiv_aids + tot_expend + status,
  data = df_logit,
  family = binomial(link = "logit")
)

## ── 9.x Model fit: Likelihood Ratio Test and Pseudo R² ------------------------

pacman::p_load(pscl)

# Null model: intercept only
log_null <- glm(
  high_life_exp ~ 1,
  data = df_logit,
  family = binomial(link = "logit")
)

# Likelihood Ratio Test: full model vs null model
lrt_log <- anova(
  log_null,
  log_model,
  test = "Chisq"
)

lrt_log

# Pseudo R² values
pseudo_r2 <- pscl::pR2(log_model)

pseudo_r2

summary(log_model)

## ── 9.3 Odds ratios with 95% CI ----------------------------------------------

log_results <- broom::tidy(
  log_model,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = dplyr::recode(term,
                         schooling        = "Schooling (years)",
                         income_comp      = "Income Composition Index",
                         log_gdp          = "log(GDP per capita)",
                         adult_mort       = "Adult Mortality",
                         hiv_aids         = "HIV/AIDS",
                         tot_expend       = "Health Expenditure",
                         statusDeveloping = "Status: Developing"
    ),
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  ) %>%
  transmute(
    Predictor = term,
    OR = round(estimate, 3),
    CI_low = round(conf.low, 3),
    CI_high = round(conf.high, 3),
    p_value = round(p.value, 4),
    sig
  )

log_results

log_results_gt <- log_results %>%
  gt() %>%
  tab_header(
    title = "Logistic Regression: Predictors of High Life Expectancy",
    subtitle = "Outcome: Life expectancy >70 years"
  ) %>%
  cols_label(
    Predictor = "Predictor",
    OR = "Odds Ratio",
    CI_low = "95% CI low",
    CI_high = "95% CI high",
    p_value = "p",
    sig = ""
  ) %>%
  tab_footnote("OR > 1 indicates higher odds of life expectancy >70 years. OR < 1 indicates lower odds.")

print(log_results_gt)

## ── 9.4 Predicted probabilities ----------------------------------------------

df_logit <- df_logit %>%
  mutate(
    predicted_prob = predict(log_model, type = "response")
  )

summary(df_logit$predicted_prob)

## ── 9.5 Example patient/country profiles -------------------------------------

profiles <- tibble(
  schooling = c(6, 10, 14),
  income_comp = c(0.45, 0.65, 0.85),
  log_gdp = c(
    quantile(df$log_gdp, 0.25),
    median(df$log_gdp),
    quantile(df$log_gdp, 0.75)
  ),
  adult_mort = c(300, 180, 80),
  hiv_aids = c(5, 1, 0.1),
  tot_expend = c(4, 6, 8),
  status = factor(
    c("Developing", "Developing", "Developed"),
    levels = levels(df$status)
  )
)

profiles$predicted_probability <- predict(
  log_model,
  newdata = profiles,
  type = "response"
)

profiles <- profiles %>%
  mutate(
    predicted_probability_percent = round(predicted_probability * 100, 1)
  )

profiles

## ── 9.6 Forest plot of odds ratios -------------------------------------------

p_log_or <- log_results %>%
  ggplot(aes(x = reorder(Predictor, OR), y = OR)) +
  geom_point(size = 3, colour = "#2166AC") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  coord_flip() +
  scale_y_log10() +
  labs(
    title = "Odds Ratios for High Life Expectancy",
    subtitle = "Outcome: Life expectancy >70 years",
    x = NULL,
    y = "Odds Ratio (log scale)"
  ) +
  theme_pub

print(p_log_or)

## ── 9.7 Model evaluation: ROC and AUC -----------------------------------------

pacman::p_load(pROC)

roc_obj <- pROC::roc(
  response = df_logit$high_life_exp,
  predictor = df_logit$predicted_prob,
  levels = c("≤70 years", ">70 years")
)

auc_value <- pROC::auc(roc_obj)
auc_value

plot(
  roc_obj,
  main = paste("ROC Curve: Logistic Regression, AUC =", round(auc_value, 3)),
  col = "#2166AC",
  lwd = 2
)
abline(a = 0, b = 1, lty = 2, col = "grey50")

## ── Table 12: Logistic regression results ------------------------------------

tbl12 <- broom::tidy(
  log_model,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Predictor = dplyr::recode(term,
                              schooling        = "Schooling (years)",
                              income_comp      = "Income Composition Index",
                              log_gdp          = "log(GDP per capita)",
                              adult_mort       = "Adult Mortality",
                              hiv_aids         = "HIV/AIDS Deaths",
                              tot_expend       = "Health Expenditure (% GDP)",
                              statusDeveloping = "Status: Developing"
    ),
    
    OR        = round(estimate, 3),
    CI        = paste0("(", round(conf.low, 3),
                       " ; ",
                       round(conf.high, 3), ")"),
    
    p_value   = round(p.value, 4),
    
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  ) %>%
  select(Predictor, OR, CI, p_value, sig)

tbl12_gt <- tbl12 %>%
  gt() %>%
  
  tab_header(
    title = "Table 12. Logistic Regression Predicting High Life Expectancy",
    subtitle = "Outcome: Life expectancy >70 years"
  ) %>%
  
  cols_label(
    Predictor = "Predictor",
    OR        = "Odds Ratio",
    CI        = "95% CI",
    p_value   = "p-value",
    sig       = ""
  ) %>%
  
  fmt_number(
    columns = c(OR, p_value),
    decimals = 3
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = p_value < 0.05
    )
  ) %>%
  
  tab_footnote(
    footnote = "OR > 1 indicates higher odds of life expectancy >70 years. OR < 1 indicates lower odds. *** p<0.001, ** p<0.01, * p<0.05."
  )

print(tbl12_gt)




## VISUALISATIONS FOR LOGISTIC REGRESSION--------
broom::glance(log_model)
# --- Figure 24. Predicted probability of high life expectancy

pred_grid <- expand.grid(
  schooling = seq(
    min(df$schooling, na.rm = TRUE),
    max(df$schooling, na.rm = TRUE),
    length.out = 100
  ),
  
  income_comp = median(df$income_comp, na.rm = TRUE),
  log_gdp     = median(df$log_gdp, na.rm = TRUE),
  adult_mort  = median(df$adult_mort, na.rm = TRUE),
  hiv_aids    = median(df$hiv_aids, na.rm = TRUE),
  tot_expend  = median(df$tot_expend, na.rm = TRUE),
  
  status = levels(df$status)
)

# Predict probabilities
pred_grid$predicted_prob <- predict(
  log_model,
  newdata = pred_grid,
  type = "response"
)

# Visualization
p_logit_prob <- ggplot(
  pred_grid,
  aes(
    x = schooling,
    y = predicted_prob,
    colour = status,
    fill = status
  )
) +
  geom_line(linewidth = 1.3) +
  
  scale_colour_manual(values = palette_status) +
  scale_fill_manual(values = palette_status) +
  
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  
  labs(
    title = "Figure 24. Predicted Probability of High Life Expectancy",
    subtitle = "Outcome: probability of life expectancy >70 years",
    x = "Years of Schooling",
    y = "Predicted Probability",
    colour = "Development Status",
    fill = "Development Status"
  ) +
  
  theme_pub

print(p_logit_prob)


# Figure 25. ROC curve and AUC for logistic regression
# Predicted probabilities
df_logit$predicted_prob <- predict(
  log_model,
  type = "response"
)

# ROC object
roc_obj <- pROC::roc(
  response = df_logit$high_life_exp,
  predictor = df_logit$predicted_prob,
  levels = c("≤70 years", ">70 years")
)

# AUC
auc_value <- pROC::auc(roc_obj)
auc_value

# ROC plot
p_roc <- ggroc(roc_obj, linewidth = 1.2, colour = "#2166AC") +
  geom_abline(
    intercept = 1,
    slope = 1,
    linetype = "dashed",
    colour = "grey50"
  ) +
  labs(
    title = "Figure 25. ROC Curve for Logistic Regression",
    subtitle = paste0("Outcome: Life expectancy >70 years | AUC = ", round(auc_value, 3)),
    x = "Specificity",
    y = "Sensitivity"
  ) +
  theme_pub

print(p_roc)

# ── Table 13: Confusion matrix — Logistic Regression
# Life expectancy >70 years

pacman::p_load(caret)

# Predicted probabilities from logistic model
df_logit$predicted_prob <- predict(
  log_model,
  type = "response"
)

# Classification threshold
df_logit$predicted_class <- ifelse(
  df_logit$predicted_prob >= 0.5,
  "1",
  "0"
)

# Convert to factors
df_logit$predicted_class <- factor(
  df_logit$predicted_class,
  levels = c("0", "1")
)

df_logit$high_le_binary <- ifelse(
  df_logit$life_exp > 70,
  "1",
  "0"
)

df_logit$high_le_binary <- factor(
  df_logit$high_le_binary,
  levels = c("0", "1")
)

# Confusion matrix
conf_matrix <- caret::confusionMatrix(
  df_logit$predicted_class,
  df_logit$high_le_binary,
  positive = "1"
)

conf_matrix

# Accuracy
accuracy <- conf_matrix$overall["Accuracy"]
accuracy


# confusion matrix table

## 2x2 confusion matrix
conf_table <- table(
  Observed  = df_logit$high_le_binary,
  Predicted = df_logit$predicted_class
)

## Convert to dataframe
conf_table_df <- as.data.frame.matrix(conf_table)

# Performance metrics
accuracy    <- round(conf_matrix$overall["Accuracy"], 3)
sensitivity <- round(conf_matrix$byClass["Sensitivity"], 3)
specificity <- round(conf_matrix$byClass["Specificity"], 3)

## Pretty table
conf_table_gt <- conf_table_df %>%
  tibble::rownames_to_column("Observed") %>%
  gt() %>%
  
  tab_header(
    title = "2x2 Confusion Matrix",
    subtitle = "Observed vs Predicted High Life Expectancy"
  ) %>%
  
  cols_label(
    Observed = "Observed",
    `0` = "Predicted 0",
    `1` = "Predicted 1"
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = c(1, 2),
      columns = c(`0`, `1`)
    )
  ) %>%
  
  tab_source_note(
    source_note = paste0(
      "Accuracy = ", accuracy,
      " | Sensitivity = ", sensitivity,
      " | Specificity = ", specificity
    )
  )

print(conf_table_gt)

