# ANALYSIS: SOCIOECONOMIC DETERMINANTS OF LIFE EXPECTANCY --------
# Research Question: What is the influence of socioeconomic factors and
#                    education on life expectancy?
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

# ── Global ggplot theme -------------------------------------------------------
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


# ── 1.1  Load & clean ---------------------------------------------------------
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

# ── 1.2  Missing data summary -------------------------------------------------
missing_summary <- df_raw %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(100 * n_missing / nrow(df_raw), 1)) %>%
  filter(n_missing > 0) %>%
  arrange(desc(n_missing))

cat("\nMissing data (variables with any NA):\n")
print(missing_summary)

# ── 1.3  Drop rows missing outcome or key predictors; impute remainder --------
df <- df_raw %>%
  filter(!is.na(life_exp), !is.na(schooling),
         !is.na(income_comp), !is.na(gdp)) %>%
  mutate(across(where(is.numeric),
                ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

cat("\nCleaned dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
cat("Status distribution:\n"); print(table(df$status))

# ── 1.4  Derived variables ----------------------------------------------------
df <- df %>%
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

socio_edu_vars <- c("life_exp", "schooling", "income_comp", "gdp",
                    "log_gdp", "tot_expend", "adult_mort",
                    "infant_deaths", "hiv_aids", "bmi", "alcohol")

# ── Table 1: Sample characteristics (gtsummary) -------------------------------
# One row per country-year; for a cleaner Table 1, use one row per country
# (most recent year) so values are not repeated across years
df_t1 <- df %>%
  group_by(country) %>%
  slice_max(year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(status, life_exp, schooling, income_comp, log_gdp,
         tot_expend, adult_mort, infant_deaths, hiv_aids, bmi, alcohol)

tbl1 <- df_t1 %>%
  tbl_summary(
    by = status,
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits  = all_continuous() ~ 2,
    label   = list(
      life_exp      ~ "Life Expectancy (years)",
      schooling     ~ "Schooling (years)",
      income_comp   ~ "Income Composition Index",
      log_gdp       ~ "log(GDP per capita)",
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
      all_continuous()  ~ "wilcox.test",
      all_categorical() ~ "chisq.test"
    ),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Table 1.** Sample characteristics by country development status (most recent year per country)")

print(tbl1)

# ── Table 2: Normality tests --------------------------------------------------
key_vars <- c("life_exp", "schooling", "income_comp", "log_gdp", "adult_mort")

tbl2 <- map(key_vars, function(v) {
  
  # 1. Remove NAs. If you don't do this, mean(), sd(), and shapiro.test() 
  # will fail or return NA if there is even a single missing value.
  x <- na.omit(df[[v]]) 
  
  # 2. Perform Shapiro-Wilk test (safely capped at 5000 to avoid errors)
  sw <- shapiro.test(sample(x, min(5000, length(x))))
  
  tibble(
    Variable  = v,
    N         = length(x), # This is the valid N (excluding NAs)
    Mean      = round(mean(x), 2),
    SD        = round(sd(x), 2),
    Skewness  = round(moments::skewness(x), 3),
    Kurtosis  = round(moments::kurtosis(x), 3),
    SW_W      = round(sw$statistic, 4),
    SW_p      = round(sw$p.value, 4),
    Normal    = ifelse(sw$p.value > 0.05, "Yes", "No")
  )
  
}) %>% 
  list_rbind() %>%   # 3. Modern replacement for map_dfr()
  mutate(Variable = c("Life Expectancy", "Schooling", "Income Comp. Index",
                      "log(GDP)", "Adult Mortality"))

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
print(tbl2_gt)

# ── Figure 1: Distribution panel ----------------------------------------------
p_hist <- df %>%
  select(all_of(key_vars)) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = dplyr::recode(Variable,
    life_exp    = "Life Expectancy",
    schooling   = "Schooling (years)",
    income_comp = "Income Comp. Index",
    log_gdp     = "log(GDP per capita)",
    adult_mort  = "Adult Mortality"
  )) %>%
  ggplot(aes(Value)) +
  geom_histogram(bins = 35, fill = "#2166AC", colour = "white", alpha = 0.85) +
  facet_wrap(~ Variable, scales = "free", ncol = 3) +
  labs(title = "Figure 1. Distribution of Key Variables",
       x = NULL, y = "Frequency")
print(p_hist)

# ── Figure 2: Violin-boxplot by status ----------------------------------------
p_violin_status <- ggplot(df, aes(status, life_exp, fill = status)) +
  geom_violin(alpha = 0.35, colour = NA) +
  geom_boxplot(width = 0.22, outlier.shape = 21, outlier.size = 1.5,
               outlier.alpha = 0.5) +
  scale_fill_manual(values = palette_status) +
  labs(title = "Figure 2. Life Expectancy by Country Development Status",
       x = "Development Status", y = "Life Expectancy (years)") +
  guides(fill = "none")
print(p_violin_status)

# ── Figure 3: Countries by status (bar) ---------------------------------------
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


# =============================================================================
# 3. CORRELATION ANALYSIS - Alena
# =============================================================================

corr_vars <- c("life_exp", "schooling", "income_comp", "log_gdp",
               "tot_expend", "adult_mort", "infant_deaths",
               "hiv_aids", "bmi", "alcohol", "hepatitis_b", "polio", "diphtheria")

corr_df      <- df %>% select(any_of(corr_vars)) %>% drop_na()
pearson_mat  <- cor(corr_df, method = "pearson",  use = "complete.obs")
spearman_mat <- cor(corr_df, method = "spearman", use = "complete.obs")

# ── Table 3: Pearson r and Spearman rho with life expectancy ------------------
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

# ── 3B.3 Pearson correlation coefficients with 95% confidence intervals --------

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


# ── 3B.4 Pearson vs Spearman comparison ---------------------------------------

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



# ── Figure 4: Pearson correlation heatmap -------------------------------------
p_corr <- ggcorrplot(pearson_mat,
                     method   = "square", type = "lower",
                     lab      = TRUE, lab_size = 2.5,
                     colors   = c("#D6604D", "white", "#2166AC"),
                     title    = "Figure 4. Pearson Correlation Matrix",
                     ggtheme  = theme_pub)
print(p_corr)

# ── Figure 5: Correlation bar chart with life expectancy ----------------------
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


<<<<<<< HEAD

# ── Pairplot / scatterplot matrix ----------------------------------------

pairplot_df <- df %>%
  select(life_exp, schooling, income_comp, log_gdp, adult_mort, hiv_aids, bmi) %>%
  drop_na()

GGally::ggpairs(
  pairplot_df,
  title = "Pairplot of Key Predictors and Life Expectancy"
)



# =============================================================================
# 4. GROUP COMPARISONS
# =============================================================================
=======
# 4. GROUP COMPARISONS ------------
>>>>>>> 0d68e05 (better subheads)

# ── Table 4: Developed vs Developing binary comparison -----------------------
# Primary: Mann-Whitney (non-normal); Welch t and Cohen's d also reported
t_res  <- t.test(life_exp ~ status, data = df, var.equal = FALSE)
mwu_res <- wilcox.test(life_exp ~ status, data = df, conf.int = TRUE)
cd_res  <- cohens_d(life_exp ~ status, data = df)

compare_vars <- c("life_exp", "schooling", "income_comp", "log_gdp",
                  "adult_mort", "hiv_aids", "tot_expend")
df <- df %>%
  mutate(across(all_of(compare_vars), as.numeric))

tbl4 <- map(compare_vars, function(v) {
  
  # 1. Extract and force to base vector
  dev   <- na.omit(df[[v]][df$status == "Developed"])
  devng <- na.omit(df[[v]][df$status == "Developing"])
  
  # 2. Wilcoxon test
  mwu <- wilcox.test(dev, devng, conf.int = FALSE)
  
  # 3. The Fix: Strip ALL hidden attributes when building temp_df
  # as.numeric(as.vector()) completely destroys any haven labels or matrix shapes
  temp_df <- tibble(
    value  = as.numeric(as.vector(c(dev, devng))), 
    status = as.character(c(rep("Developed", length(dev)), rep("Developing", length(devng))))
  )
  
  # 4. Calculate Cohen's d
  cd <- rstatix::cohens_d(temp_df, value ~ status)
  
  tibble(
    Variable     = v,
    Dev_Median   = round(median(dev), 2),
    Dev_IQR      = paste0("(", round(quantile(dev, .25), 1), "-",
                          round(quantile(dev, .75), 1), ")"),
    Devng_Median = round(median(devng), 2),
    Devng_IQR    = paste0("(", round(quantile(devng, .25), 1), "-",
                          round(quantile(devng, .75), 1), ")"),
    W            = round(mwu$statistic, 1),
    p_value      = round(mwu$p.value, 4),
    Cohens_d     = round(abs(cd$effsize), 3)
  )
  
}) %>% 
  bind_rows() %>%
  mutate(
    sig = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                    p_value < 0.05  ~ "*",   TRUE           ~ ""),
    Variable = dplyr::recode(Variable,
                             life_exp    = "Life Expectancy (years)", schooling   = "Schooling (years)",
                             income_comp = "Income Comp. Index",      log_gdp     = "log(GDP per capita)",
                             adult_mort  = "Adult Mortality",         hiv_aids    = "HIV/AIDS Deaths",
                             tot_expend  = "Health Expenditure (% GDP)"
    )
  )

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

# ── Table 5: Categorical associations ----------------------------------------
chi_tbl    <- table(df$status, df$schooling_cat)
chi_res    <- chisq.test(chi_tbl)
cv_school  <- effectsize::cramers_v(chi_tbl)

fisher_tbl <- table(df$status, df$life_exp_cat)
fish_res   <- fisher.test(fisher_tbl, simulate.p.value = TRUE, B = 10000)

tbl5 <- tibble(
  Test       = c("Chi-square", "Fisher's Exact"),
  Comparison = c("Status x Schooling Category", "Status x Life Expectancy Category"),
  Statistic  = c(round(chi_res$statistic, 3), NA_real_),
  df         = c(chi_res$parameter, NA_real_),
  p_value    = c(round(chi_res$p.value, 4), round(fish_res$p.value, 4)),
  Effect     = c(paste0("Cramer's V = ", round(cv_school$Cramers_v, 3)), "—")
)

tbl5_gt <- tbl5 %>%
  gt() %>%
  tab_header(title = "Table 5. Categorical Association Tests") %>%
  fmt_number(columns = c(Statistic), decimals = 3) %>%
  fmt_number(columns = p_value, decimals = 4) %>%
  sub_missing(missing_text = "—") %>%
  tab_footnote("Fisher's Exact p-value via Monte Carlo simulation (B=10,000).")
print(tbl5_gt)

# ── 4.1  One-way KW + Dunn: Life Expectancy ~ Schooling category --------------
cat("\n-- Kruskal-Wallis: Life Expectancy ~ Schooling Category --\n")
kw_school <- kruskal.test(life_exp ~ schooling_cat, data = df)
print(kw_school)
cat("\nDunn post-hoc (Bonferroni):\n")
dunn_school <- dunn.test(df$life_exp, df$schooling_cat, method = "bonferroni",
                         alph = 0.05, list = FALSE)

cat("\n-- Kruskal-Wallis: Life Expectancy ~ Income Category --\n")
kw_income <- kruskal.test(life_exp ~ income_cat, data = df)
print(kw_income)
cat("\nDunn post-hoc (Bonferroni):\n")
dunn_income <- dunn.test(df$life_exp, df$income_cat, method = "bonferroni",
                         alph = 0.05, list = FALSE)

# ── 4.2  Two-way ANOVA: Status x Schooling ------------------------------------
cat("\n-- Two-way ANOVA: Life Expectancy ~ Status x Schooling Category --\n")
anova_2way  <- aov(life_exp ~ status * schooling_cat, data = df)
print(summary(anova_2way))
eta_2way <- effectsize::eta_squared(anova_2way)
eta_2way
cat("\nPartial eta-squared:\n"); print(eta_2way)

emm_2way <- emmeans(anova_2way, pairwise ~ status | schooling_cat,
                    adjust = "tukey")
cat("\nEmmeans pairwise (Status within Schooling):\n")
print(emm_2way$contrasts)

# ── Figure 7: Boxplot by schooling category -----------------------------------
p_box_school <- ggplot(df, aes(schooling_cat, life_exp, fill = schooling_cat)) +
  geom_violin(alpha = 0.3, colour = NA) +
  geom_boxplot(width = 0.28, outlier.shape = 21, outlier.size = 1.2,
               outlier.alpha = 0.4) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Figure 6. Life Expectancy by Schooling Level",
       x = "Schooling Tertile", y = "Life Expectancy (years)")
print(p_box_school)

# ── Figure 8: Boxplot by income category --------------------------------------
p_box_income <- ggplot(df, aes(income_cat, life_exp, fill = income_cat)) +
  geom_violin(alpha = 0.3, colour = NA) +
  geom_boxplot(width = 0.28, outlier.shape = 21, outlier.size = 1.2,
               outlier.alpha = 0.4) +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  labs(title = "Figure 7. Life Expectancy by Income Composition Category",
       x = "Income Tertile", y = "Life Expectancy (years)")
print(p_box_income)

# ── Figure 9: Two-way interaction plot ----------------------------------------
p_twoway <- ggplot(df, aes(schooling_cat, life_exp, fill = status, colour = status)) +
  geom_boxplot(alpha = 0.45, position = position_dodge(0.8),
               outlier.shape = 21, outlier.size = 1.2, outlier.alpha = 0.4) +
  scale_fill_manual(values = palette_status) +
  scale_colour_manual(values = palette_status) +
  labs(title    = "Figure 8. Life Expectancy by Schooling x Development Status",
       x        = "Schooling Tertile",
       y        = "Life Expectancy (years)",
       fill     = "Status", colour = "Status")
print(p_twoway)

# ── Figure 10: Mosaic-style stacked bar ----------------------------------------
p_mosaic <- df %>%
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

# ── 4.3  Repeated-measures ANOVA (2000, 2005, 2010, 2015) --------------------
cat("\n-- Repeated-Measures ANOVA: Life Expectancy over Time --\n")
rm_data <- df %>%
  filter(year %in% c(2000, 2005, 2010, 2015)) %>%
  group_by(country) %>%
  filter(n() == 4) %>%
  ungroup() %>%
  mutate(year_f = factor(year), country = factor(country))

cat("Countries with all 4 time points:", length(unique(rm_data$country)), "\n")
rm_model <- aov(life_exp ~ year_f + Error(country / year_f), data = rm_data)
print(summary(rm_model))



# 5. TEMPORAL ANALYSIS ---------


yearly_summary <- df %>%
  group_by(year, status) %>%
  summarise(
    avg_life_exp   = mean(life_exp, na.rm = TRUE),
    avg_schooling  = mean(schooling, na.rm = TRUE),
    avg_income     = mean(income_comp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

# ── Figure 11: Life expectancy gap over time ----------------------------------
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

# ── Table 6: Year-by-year Mann-Whitney gap analysis ---------------------------
test_vars <- c("life_exp", "infant_deaths", "polio")

run_yearly_tests <- function(data, variable) {
  data %>%
    group_by(year) %>%
    summarise(
      Median_Developed  = round(median(.data[[variable]][status == "Developed"], na.rm = TRUE), 2),
      Median_Developing = round(median(.data[[variable]][status == "Developing"], na.rm = TRUE), 2),
      Gap               = round(abs(Median_Developed - Median_Developing), 2),
      p_value           = wilcox.test(.data[[variable]] ~ status, exact = FALSE)$p.value,
      .groups           = "drop"
    ) %>%
    mutate(
      Variable    = variable,
      p_value     = signif(p_value, 3),
      Significant = ifelse(p_value < 0.05, "Yes", "No")
    ) %>%
    select(Variable, year, Median_Developed, Median_Developing, Gap, p_value, Significant)
}

yearly_hypothesis_results <- bind_rows(lapply(test_vars, function(v) run_yearly_tests(df, v)))

tbl6_gt <- yearly_hypothesis_results %>%
  mutate(Variable = dplyr::recode(Variable,
    life_exp      = "Life Expectancy",
    infant_deaths = "Infant Deaths",
    polio         = "Polio"
  )) %>%
  gt(groupname_col = "Variable") %>%
  tab_header(title = "Table 6. Mann-Whitney U Tests by Year: Developed vs Developing") %>%
  cols_label(year = "Year", Median_Developed = "Median (Developed)",
             Median_Developing = "Median (Developing)", Gap = "Absolute Gap",
             p_value = "p-value", Significant = "Sig.") %>%
  fmt_number(columns = c(Median_Developed, Median_Developing, Gap), decimals = 2) %>%
  fmt_number(columns = p_value, decimals = 4) %>%
  tab_style(
    style = cell_text(color = "red", weight = "bold"),
    locations = cells_body(columns = Significant, rows = Significant == "Yes")
  ) %>%
  tab_footnote("Wilcoxon rank-sum test per year. Significant = p < 0.05.")
print(tbl6_gt)

# ── Figure 12: Gap significance over time -------------------------------------
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


scatter_data <- df %>%
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

# ── Model building (hierarchical) --------------------------------------------
m0 <- lm(life_exp ~ 1, data = df)
m1 <- lm(life_exp ~ log_gdp + income_comp + tot_expend + status, data = df)
m2 <- lm(life_exp ~ log_gdp + income_comp + tot_expend + status + schooling, data = df)
m3 <- lm(life_exp ~ log_gdp + income_comp + tot_expend + status +
            schooling + adult_mort + infant_deaths + hiv_aids +
            bmi + alcohol + hepatitis_b + polio + diphtheria, data = df)
m4 <- lm(life_exp ~ income_comp + schooling + adult_mort + hiv_aids +
            bmi + status, data = df)

# ── Table 7: Model fit statistics --------------------------------------------
fchange_m2 <- anova(m1, m2)
fchange_m3 <- anova(m2, m3)

tbl7 <- tibble(
  Model      = c("M1: Socioeconomic", "M2: + Education", "M3: Full", "M4: Parsimonious"),
  R2         = round(c(summary(m1)$r.squared,  summary(m2)$r.squared,
                       summary(m3)$r.squared,  summary(m4)$r.squared), 3),
  Adj_R2     = round(c(summary(m1)$adj.r.squared, summary(m2)$adj.r.squared,
                       summary(m3)$adj.r.squared, summary(m4)$adj.r.squared), 3),
  AIC        = round(c(AIC(m1), AIC(m2), AIC(m3), AIC(m4)), 1),
  BIC        = round(c(BIC(m1), BIC(m2), BIC(m3), BIC(m4)), 1),
  RMSE       = round(c(sqrt(mean(m1$residuals^2)), sqrt(mean(m2$residuals^2)),
                       sqrt(mean(m3$residuals^2)), sqrt(mean(m4$residuals^2))), 3),
  F_change_p = c(NA_real_, round(fchange_m2$`Pr(>F)`[2], 4),
                 round(fchange_m3$`Pr(>F)`[2], 4), NA_real_)
)

tbl7_gt <- tbl7 %>%
  gt() %>%
  tab_header(title = "Table 7. Linear Regression Model Fit Comparison") %>%
  cols_label(Adj_R2 = "Adj. R\u00b2", F_change_p = "F-change p") %>%
  fmt_number(columns = c(R2, Adj_R2, RMSE), decimals = 3) %>%
  fmt_number(columns = c(AIC, BIC), decimals = 1) %>%
  fmt_number(columns = F_change_p, decimals = 4) %>%
  sub_missing(missing_text = "—") %>%
  tab_style(
    style = cell_fill(color = "#EBF5FB"),
    locations = cells_body(rows = Model == "M4: Parsimonious")
  ) %>%
  tab_footnote("Highlighted row (M4) selected as parsimonious final model.")
print(tbl7_gt)

# ── Table 8: Regression coefficients M3 + M4 side-by-side --------------------
coef_m3 <- tidy(m3, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  transmute(term,
            b_m3   = round(estimate, 3),
            ci_m3  = paste0("(", round(conf.low, 3), ", ", round(conf.high, 3), ")"),
            p_m3   = round(p.value, 4),
            sig_m3 = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                               p.value < 0.05  ~ "*",   TRUE           ~ ""))

coef_m4 <- tidy(m4, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  transmute(term,
            b_m4   = round(estimate, 3),
            ci_m4  = paste0("(", round(conf.low, 3), ", ", round(conf.high, 3), ")"),
            p_m4   = round(p.value, 4),
            sig_m4 = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                               p.value < 0.05  ~ "*",   TRUE           ~ ""))

std_m4 <- parameters::model_parameters(m4, standardize = "refit") %>%
  as.data.frame() %>%
  filter(Parameter != "(Intercept)") %>%
  transmute(term = Parameter, beta_std = round(Coefficient, 3))

tbl8 <- coef_m3 %>%
  full_join(coef_m4, by = "term") %>%
  left_join(std_m4,  by = "term") %>%
  mutate(term = dplyr::recode(term,
    income_comp        = "Income Comp. Index",
    schooling          = "Schooling (years)",
    adult_mort         = "Adult Mortality",
    hiv_aids           = "HIV/AIDS Deaths",
    bmi                = "BMI",
    statusDeveloping   = "Status: Developing",
    log_gdp            = "log(GDP)",
    tot_expend         = "Health Expenditure",
    infant_deaths      = "Infant Deaths",
    alcohol            = "Alcohol",
    hepatitis_b        = "Hepatitis B",
    polio              = "Polio",
    diphtheria         = "Diphtheria"
  ))

tbl8_gt <- tbl8 %>%
  gt() %>%
  tab_header(title    = "Table 8. Multiple Linear Regression: Predictors of Life Expectancy",
             subtitle = "Unstandardised coefficients with 95% CI") %>%
  tab_spanner(label = "M3: Full model",        columns = c(b_m3, ci_m3, p_m3, sig_m3)) %>%
  tab_spanner(label = "M4: Parsimonious",      columns = c(b_m4, ci_m4, p_m4, sig_m4)) %>%
  tab_spanner(label = "Std. (M4)",             columns = beta_std) %>%
  cols_label(term = "Predictor",
             b_m3 = "b", ci_m3 = "95% CI", p_m3 = "p", sig_m3 = "",
             b_m4 = "b", ci_m4 = "95% CI", p_m4 = "p", sig_m4 = "",
             beta_std = "\u03b2") %>%
  fmt_number(columns = c(b_m3, b_m4, beta_std), decimals = 3) %>%
  fmt_number(columns = c(p_m3, p_m4), decimals = 4) %>%
  sub_missing(missing_text = "—") %>%
  tab_footnote("*** p<0.001, ** p<0.01, * p<0.05. beta = standardised coefficient.")
print(tbl8_gt)

# ── Table 9: Model assumption checks -----------------------------------------
bp_res  <- bptest(m3)
dw_res  <- dwtest(m3)
vif_res <- vif(m3)

tbl9_vif <- tibble(
  Predictor = names(vif_res),
  VIF       = round(vif_res, 3),
  Flag      = ifelse(vif_res > 5, "High", "OK")
) %>%
  gt() %>%
  tab_header(title = "Table 9a. Variance Inflation Factors (Model 3)") %>%
  tab_style(
    style = cell_text(color = "red", weight = "bold"),
    locations = cells_body(columns = Flag, rows = Flag == "High")
  ) %>%
  tab_footnote("VIF > 5 flagged as potentially problematic multicollinearity.")
print(tbl9_vif)

tbl9_tests <- tibble(
  Test        = c("Breusch-Pagan", "Durbin-Watson"),
  Statistic   = c(round(bp_res$statistic, 4), round(dw_res$statistic, 4)),
  p_value     = c(round(bp_res$p.value, 4), round(dw_res$p.value, 4)),
  Null        = c("Homoscedasticity", "No autocorrelation"),
  Conclusion  = c(
    ifelse(bp_res$p.value < 0.05, "Heteroscedasticity present", "Homoscedastic"),
    ifelse(dw_res$p.value < 0.05, "Autocorrelation present", "No autocorrelation")
  )
) %>%
  gt() %>%
  tab_header(title = "Table 9b. Regression Assumption Tests (Model 3)") %>%
  fmt_number(columns = c(Statistic, p_value), decimals = 4)
print(tbl9_tests)

# ── Figures: regression diagnostics ------------------------------------------
par(mfrow = c(2, 2))
plot(m3, main = "Figure 15. Regression Diagnostics (Model 3)")
par(mfrow = c(1, 1))

pred_df <- df %>%
  mutate(predicted = fitted(m3), resid = residuals(m3))

p_avp <- ggplot(pred_df, aes(predicted, life_exp, colour = status)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_manual(values = palette_status) +
  labs(title = "Figure 16. Actual vs Predicted Life Expectancy (Model 3)",
       x = "Predicted", y = "Actual", colour = "Status")
print(p_avp)

p_resid <- ggplot(pred_df, aes(predicted, resid)) +
  geom_point(alpha = 0.3, size = 1.4, colour = "#2166AC") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, colour = "#D6604D") +
  labs(title = "Figure 17. Residuals vs Fitted (Model 3)",
       x = "Fitted Values", y = "Residuals")
print(p_resid)

std_coefs <- parameters::model_parameters(m4, standardize = "refit") %>%
  as.data.frame() %>%
  filter(Parameter != "(Intercept)") %>%
  mutate(across(where(is.numeric), ~ round(., 3)),
         Parameter = dplyr::recode(Parameter,
           income_comp      = "Income Comp. Index", schooling  = "Schooling",
           adult_mort       = "Adult Mortality",    hiv_aids   = "HIV/AIDS",
           bmi              = "BMI",                statusDeveloping = "Status: Developing"
         ))

p_std_coef <- std_coefs %>%
  ggplot(aes(reorder(Parameter, Coefficient), Coefficient,
             fill = Coefficient > 0)) +
  geom_col(colour = "white") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D6604D"),
                    guide = "none") +
  labs(title    = "Figure 18. Standardised Coefficients (Model 4)",
       subtitle = "Beta with 95% CI",
       x = NULL, y = "Standardised beta")
print(p_std_coef)

# ── 7A. Exploratory model plots: linear vs non-linear trends -------------------

key_model_pairs <- list(
  list(x = "schooling",   xlab = "Schooling (years)"),
  list(x = "income_comp", xlab = "Income Composition Index"),
  list(x = "log_gdp",     xlab = "log(GDP per capita)"),
  list(x = "hiv_aids",    xlab = "HIV/AIDS deaths")
)

for (pair in key_model_pairs) {
  
  p <- ggplot(df, aes(x = .data[[pair$x]], y = life_exp)) +
    geom_point(alpha = 0.35, colour = "#2166AC") +
    geom_smooth(method = "lm", se = TRUE, colour = "#D6604D") +
    geom_smooth(method = "loess", se = FALSE, colour = "black", linetype = "dashed") +
    labs(
      title = paste("Life Expectancy vs", pair$xlab),
      subtitle = "Red = linear trend; dashed black = LOESS trend",
      x = pair$xlab,
      y = "Life Expectancy (years)"
    ) +
    theme_pub
  
  print(p)
}



# ── 7B. Partial correlations: adjusted associations ---------------------------

partial_schooling_gdp <- ppcor::pcor.test(
  x = df$schooling,
  y = df$life_exp,
  z = df$log_gdp
)

partial_income_gdp <- ppcor::pcor.test(
  x = df$income_comp,
  y = df$life_exp,
  z = df$log_gdp
)

partial_results <- tibble(
  Relationship = c(
    "Schooling vs Life Expectancy adjusted for log(GDP)",
    "Income Composition vs Life Expectancy adjusted for log(GDP)"
  ),
  Partial_r = c(
    partial_schooling_gdp$estimate,
    partial_income_gdp$estimate
  ),
  p_value = c(
    partial_schooling_gdp$p.value,
    partial_income_gdp$p.value
  )
) %>%
  mutate(
    Partial_r = round(Partial_r, 3),
    p_value = round(p_value, 4),
    Significant = ifelse(p_value < 0.05, "Yes", "No")
  )

partial_results



# 8. INTERACTION / MODERATION ANALYSIS ---------

# ── 8.1  Schooling x Status --------------------------------------------------
m_add1 <- lm(life_exp ~ schooling + status + income_comp + adult_mort + hiv_aids, data = df)
m_int1 <- lm(life_exp ~ schooling * status + income_comp + adult_mort + hiv_aids, data = df)
lrt1   <- anova(m_add1, m_int1)
cat("\n-- LRT: Schooling x Status --\n"); print(lrt1)

emm_int1 <- emmeans(m_int1, ~ schooling | status,
                    at = list(schooling = seq(0, 20, by = 2)))

p_emm1 <- as.data.frame(emm_int1) %>%
  ggplot(aes(schooling, emmean, colour = status, fill = status)) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = palette_status) +
  scale_fill_manual(values   = palette_status) +
  labs(title    = "Figure 19. Marginal Effect of Schooling on Life Expectancy by Status",
       subtitle = "Schooling x Development Status interaction",
       x = "Years of Schooling", y = "Estimated Life Expectancy (years)",
       colour = "Status", fill = "Status")
print(p_emm1)

# ── 8.2  Income x Status -----------------------------------------------------
m_add2 <- lm(life_exp ~ income_comp + status + schooling + adult_mort + hiv_aids, data = df)
m_int2 <- lm(life_exp ~ income_comp * status + schooling + adult_mort + hiv_aids, data = df)
lrt2   <- anova(m_add2, m_int2)
cat("\n-- LRT: Income x Status --\n"); print(lrt2)

emm_int2 <- emmeans(m_int2, ~ income_comp | status,
                    at = list(income_comp = seq(0.2, 0.9, by = 0.05)))

p_emm2 <- as.data.frame(emm_int2) %>%
  ggplot(aes(income_comp, emmean, colour = status, fill = status)) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = palette_status) +
  scale_fill_manual(values   = palette_status) +
  labs(title    = "Figure 20. Marginal Effect of Income Composition on Life Expectancy by Status",
       subtitle = "Income Composition x Development Status interaction",
       x = "Income Composition Index", y = "Estimated Life Expectancy (years)",
       colour = "Status", fill = "Status")
print(p_emm2)

# ── 8.3  Schooling x Income category -----------------------------------------
m_add3 <- lm(life_exp ~ schooling + income_cat + adult_mort + hiv_aids + status, data = df)
m_int3 <- lm(life_exp ~ schooling * income_cat + adult_mort + hiv_aids + status, data = df)
lrt3   <- anova(m_add3, m_int3)
cat("\n-- LRT: Schooling x Income Category --\n"); print(lrt3)

emm_int3 <- emmeans(m_int3, ~ schooling | income_cat,
                    at = list(schooling = seq(0, 20, by = 2)))

p_emm3 <- as.data.frame(emm_int3) %>%
  ggplot(aes(schooling, emmean, colour = income_cat, fill = income_cat)) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.2) +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette   = "Set1") +
  labs(title    = "Figure 21. Marginal Effect of Schooling by Income Tertile",
       subtitle = "Schooling x Income Category interaction",
       x = "Years of Schooling", y = "Estimated Life Expectancy (years)",
       colour = "Income Category", fill = "Income Category")
print(p_emm3)

# ── Figure 22: Heatmap (Schooling x Income) ----------------------------------
heatmap_df <- expand.grid(
  schooling   = seq(min(df$schooling), max(df$schooling), length.out = 50),
  income_comp = seq(min(df$income_comp), max(df$income_comp), length.out = 50)
) %>%
  mutate(
    status     = "Developing",
    income_cat = cut(income_comp,
                     breaks = quantile(df$income_comp, c(0, 0.25, 0.75, 1)),
                     labels = c("Low", "Middle", "High"),
                     include.lowest = TRUE),
    adult_mort = median(df$adult_mort),
    hiv_aids   = median(df$hiv_aids),
    predicted  = predict(m_int3, newdata = pick(everything()))
  )

p_heatmap <- ggplot(heatmap_df, aes(schooling, income_comp, fill = predicted)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", name = "Predicted\nLife Exp.") +
  labs(title    = "Figure 22. Predicted Life Expectancy: Schooling x Income",
       subtitle = "Adult mortality and HIV/AIDS held at median (Developing countries)",
       x = "Years of Schooling", y = "Income Composition Index")
print(p_heatmap)

# ── Table 10: Logistic regression — high life expectancy (>70) ---------------
df_logit <- df %>% mutate(high_le = as.integer(life_exp > 70))
log_m1   <- glm(high_le ~ schooling + income_comp + log_gdp +
                  adult_mort + hiv_aids + status + tot_expend,
                data = df_logit, family = binomial(link = "logit"))

tbl10 <- tidy(log_m1, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    sig  = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                     p.value < 0.05  ~ "*",   TRUE            ~ ""),
    term = dplyr::recode(term,
      schooling        = "Schooling (years)", income_comp = "Income Comp. Index",
      log_gdp          = "log(GDP)",          adult_mort  = "Adult Mortality",
      hiv_aids         = "HIV/AIDS Deaths",   statusDeveloping = "Status: Developing",
      tot_expend       = "Health Expenditure"
    )
  ) %>%
  select(term, estimate, conf.low, conf.high, p.value, sig) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ round(., 3)),
         p.value = round(p.value, 4))

tbl10_gt <- tbl10 %>%
  gt() %>%
  tab_header(title    = "Table 10. Logistic Regression: Predictors of High Life Expectancy (>70 years)",
             subtitle = "Odds ratios with 95% confidence intervals") %>%
  cols_label(term = "Predictor", estimate = "OR", conf.low = "95% CI (low)",
             conf.high = "95% CI (high)", p.value = "p", sig = "") %>%
  fmt_number(columns = c(estimate, conf.low, conf.high), decimals = 3) %>%
  fmt_number(columns = p.value, decimals = 4) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = p.value < 0.05)
  ) %>%
  tab_footnote("*** p<0.001, ** p<0.01, * p<0.05. Bold rows = statistically significant.")
print(tbl10_gt)

# OR forest plot
p_or <- tbl10 %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate,
             colour = p.value < 0.05)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  coord_flip() +
  scale_y_log10() +
  scale_colour_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#BDBDBD")) +
  labs(title    = "Figure 23. Odds Ratios: High Life Expectancy (>70 years)",
       subtitle = "Log scale; dashed line = OR 1.0",
       x = NULL, y = "Odds Ratio (log scale)", colour = "p < 0.05")
print(p_or)

# ── Table 11: Interaction summary ---------------------------------------------
lrt_summary <- function(lrt_obj, label) {
  tibble(
    Interaction  = label,
    F_statistic  = round(lrt_obj$F[2], 3),
    df_num       = lrt_obj$Df[2],
    df_den       = lrt_obj$Res.Df[2],
    p_value      = round(lrt_obj$`Pr(>F)`[2], 4),
    Significant  = ifelse(lrt_obj$`Pr(>F)`[2] < 0.05, "Yes", "No")
  )
}

tbl11 <- bind_rows(
  lrt_summary(lrt1, "Schooling x Development Status"),
  lrt_summary(lrt2, "Income Composition x Development Status"),
  lrt_summary(lrt3, "Schooling x Income Category")
)

tbl11_gt <- tbl11 %>%
  gt() %>%
  tab_header(title = "Table 11. Likelihood-Ratio Tests for Interaction Terms") %>%
  cols_label(F_statistic = "F", df_num = "df (num)", df_den = "df (den)",
             p_value = "p-value", Significant = "Sig.") %>%
  fmt_number(columns = F_statistic, decimals = 3) %>%
  fmt_number(columns = p_value, decimals = 4) %>%
  tab_style(
    style = cell_text(color = "#D6604D", weight = "bold"),
    locations = cells_body(columns = Significant, rows = Significant == "Yes")
  ) %>%
  tab_footnote("F-test from ANOVA comparison of additive vs interaction model.")
print(tbl11_gt)


# =============================================================================
# 9. Logistic regression
# =============================================================================

# ── 9.1 Create binary outcome -------------------------------------------------

df_logit <- df %>%
  mutate(
    high_life_exp = ifelse(life_exp > 70, 1, 0),
    high_life_exp = factor(high_life_exp, levels = c(0, 1),
                           labels = c("≤70 years", ">70 years"))
  )

table(df_logit$high_life_exp)

# ── 9.2 Fit logistic regression model ----------------------------------------

log_model <- glm(
  high_life_exp ~ schooling + income_comp + log_gdp +
    adult_mort + hiv_aids + tot_expend + status,
  data = df_logit,
  family = binomial(link = "logit")
)

# ── 9.x Model fit: Likelihood Ratio Test and Pseudo R² ------------------------

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

# ── 9.3 Odds ratios with 95% CI ----------------------------------------------

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

# ── 9.4 Predicted probabilities ----------------------------------------------

df_logit <- df_logit %>%
  mutate(
    predicted_prob = predict(log_model, type = "response")
  )

summary(df_logit$predicted_prob)

# ── 9.5 Example patient/country profiles -------------------------------------

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

# ── 9.6 Forest plot of odds ratios -------------------------------------------

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

# ── 9.7 Model evaluation: ROC and AUC -----------------------------------------

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

# ── Table 12: Logistic regression results ------------------------------------

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

broom::glance(log_model)


## VISUALISATIONS FOR LOGISTIC REGRESSION
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
# ============================================================

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