###############################################################################
# Title: Script to Produce Results (Tables and Figures) for our Manuscript
# Manuscript Title: "Iodized salt has large impacts on child health and test scores in Ethiopia"
# Authors: Robel Alemu, Kibrom Tafere, Dawd Gashu, Edward J. M. Joy, Elizabeth H. Bailey, 
#          R. Murray Lark, Martin R. Broadley, William A. Masters
# Date: This version was revised on 10 December 2024
#
# Description:
# This script demonstrates the steps taken to produce key tables and figures
# for the manuscript "Iodized salt has large impacts on child health and test scores in Ethiopia."
#
# The script includes:
# 1. Pre-processing of the dataset, including:
#    - Computing total test scores.
#    - Normalizing test scores into z-scores by centering on the mean and scaling 
#      by the standard deviation of students who were born in the same year, took the exam 
#      in the same year, and resided in the same administrative region.
#    - Creating standardized (z-score) variables for simulated nutrient measures.
#    - Applying filters to restrict the sample to the intended study population 
#      (e.g., rural districts outside of the Oromia special zone surrounding Addis Ababa).
#
# 2. Running linear fixed-effects regressions (felm) to estimate treatment effects
#    for children born after loss of access to iodized salt (main analysis), and a placebo test 
#    for urban districts around Addis Ababa (where food is largely not locally sourced).
#
# 3. Formatting the regression results into a publication-ready table (Table 1).
#    Table 1 shows the dose-response of secondary-school exam scores to local soil iodine 
#    for children born after loss of access to iodized salt in rural districts (treatment effect),
#    and for children in urban districts around Addis Ababa (placebo test).
#
# In addition, we also test for placebo effects using soil selenium instead of soil iodine.
# This falsification test, where selenium replaces iodine as the local nutrient measure, 
# is presented in Supplementary Table S5:
#
# "Table S5. Placebo test of dose-response in exam scores to local soil selenium 
#  for children born after loss of access to iodized salt"
#
# Note: This placebo test uses the same modeling approach described above, but 
# replaces the district-level soil iodine variable with soil selenium. If the 
# observed relationships with iodine were purely artifactual (due to model specification), 
# we would expect to find similar significant effects with selenium. The absence of 
# effects with selenium supports the validity of the iodine results.
#
# Note:
# - This script uses placeholder directories for the input and output data.
# - You should customize the paths to point to your actual data directories.
# - Ensure that the datasets referenced (e.g., cohort1_1yearwindow, 
#   cohort2_2yearwindow, cohort3_3yearwindow) are loaded prior to running the models.
#
###############################################################################

# ---------------------------- Setup -------------------------------------------
# Clear the global environment
rm(list = ls())

# Install and load necessary libraries (if needed).
install.packages(c("dplyr", "tidyverse", "broom", "data.table", 
                   "lme4", "lfe", "haven", 
                   "lmtest", "fixest", "broom", "openxlsx", 
                   "sandwich", "doParallel", "car", "fmsb", 
                   "boot", "ggplot2", "stringr", "devtools"), 
                 dependencies = TRUE)

library(haven)
library(lmtest)
library(sandwich)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(boot)
library(conflicted)
library(fixest)
library(lfe)
library(broom)
library(openxlsx)
library(car)
library(fmsb)
library(data.table)

# ---------------------------- Directories -------------------------------------
# Set your working directory and define paths to your data and output directories.
# IMPORTANT: Update these paths as needed. These are placeholder directories.
data_dir <- "path/to/data_directory/"
output_dir <- "path/to/output_directory/"

# ---------------------------- Data Loading ------------------------------------
# Load the main dataset.
# Replace with your actual data file name and directory.
TestScore_data <- haven::read_dta(file.path(data_dir, "TestScore_only_04May2021.dta"))
print(colnames(TestScore_data))

# ---------------------------- Data Pre-processing ------------------------------
#
# STEP 1: Compute total test scores in two domains: social and natural.
TestScore_data <- TestScore_data %>%
  mutate(
    total_social = english + aptitude + civics + economics + socialmath + geography + history,
    total_natural = english + aptitude + civics + naturalmath + physics + biology + chemistry
  )

# STEP 2: Normalize test scores into z-scores by centering on the mean and scaling by 
# the standard deviation of students who were born in the same year, took the exam 
# in the same year, and resided in the same administrative region.
TestScore_data <- TestScore_data %>%
  group_by(yeargc, birth_yearGC, region_code) %>%
  mutate(across(
    c(english, aptitude, civics, total_3Sub, physics, chemistry, biology, 
      geography, history, economics, naturalmath, socialmath, total_social, 
      total_natural),
    ~ ( . - mean(., na.rm = TRUE) ) / sd(., na.rm = TRUE),
    .names = "{.col}_std"
  )) %>%
  ungroup()

# STEP 3: Convert Simulated_OrgIodine and Simulated_SolubleSe to z-scores.
TestScore_data$Simulated_OrgIodine_ZScore <- scale(TestScore_data$Simulated_OrgIodine)
TestScore_data$Simulated_SolubleSe_ZScore <- scale(TestScore_data$Simulated_SolubleSe)

# STEP 4: Identify observations in the Oromia Special Zone Surrounding Finfinne (near Addis).
TestScore_data <- TestScore_data %>%
  mutate(urban_nearAddis = ifelse(zone_name == "Oromia Special Zone Surrounding Finfinne", 1, 0))

summary(TestScore_data$urban_nearAddis)

# STEP 5: Restrict the sample to rural districts excluding the special zone near Addis
TestScore_data <- dplyr::filter(TestScore_data, urban_nearAddis == 0)

# Summarize the data by region to ensure filters were applied correctly.
region_summary <- TestScore_data %>%
  group_by(region_name) %>%
  summarize(count = n())

print(region_summary)
summary(as.factor(TestScore_data$urban_nearAddis))

# At this point, the data (TestScore_data) is filtered and normalized according 
# to the study's criteria for the rural treatment sample. 
#
# For the placebo test (urban districts around Addis Ababa), you would implement 
# similar steps but filter differently (e.g., keep urban_nearAddis == 1) before 
# running analogous regressions.

# For the selenium placebo test (Supplementary Table S5), you would follow the 
# same procedure but replace the iodine variable (Simulated_OrgIodine_ZScore) 
# with the selenium variable (Simulated_SolubleSe_ZScore) in the regression models.

# --------------------------------------------------------------------------------
# Load or define your cohort datasets here. For example:
#
# cohort1_1yearwindow <- ...
# cohort2_2yearwindow <- ...
# cohort3_3yearwindow <- ...
#
# These datasets should be subsets of TestScore_data or another source, containing:
# total_3Sub_std, Post, Simulated_OrgIodine_ZScore (or Simulated_SolubleSe_ZScore for placebo), 
# female, natural_stream, woreda_code, schoolcode, yeargc, birth_yearGC, 
# Simulated_OrgIodine, total_3Sub.
#
# Adjust these data loading steps as necessary.

# ---------------------------- Custom Function -----------------------------------
tidy_felm <- function(model) {
  coefs <- summary(model, robust = TRUE)$coefficients
  tibble(
    term = rownames(coefs),
    estimate = coefs[, "Estimate"],
    std.error = coefs[, "Cluster s.e."],
    statistic = coefs[, "t value"],
    p.value = coefs[, "Pr(>|t|)"]
  )
}

# ---------------------------- Regression Models --------------------------------
# Below is an example of running the FELM regressions for the rural sample with iodine.
# For the placebo test using selenium (Supplementary Table S5), replace 
# "Simulated_OrgIodine_ZScore" with "Simulated_SolubleSe_ZScore" in the model formula.
#
# For the urban placebo test (using iodine or selenium), run the same regressions on 
# the urban sample (TestScore_data filtered for urban_nearAddis == 1).

model_1 <- felm(
  total_3Sub_std ~ Post + Post:Simulated_OrgIodine_ZScore + female + natural_stream 
  | factor(woreda_code) + factor(schoolcode) + factor(yeargc) 
  | 0 | woreda_code, 
  data = cohort1_1yearwindow
)

model_2 <- felm(
  total_3Sub_std ~ Post + Post:Simulated_OrgIodine_ZScore + female + natural_stream 
  | factor(woreda_code) + factor(schoolcode) + factor(yeargc) 
  | 0 | woreda_code, 
  data = cohort2_2yearwindow
)

model_3 <- felm(
  total_3Sub_std ~ Post + Post:Simulated_OrgIodine_ZScore + female + natural_stream 
  | factor(woreda_code) + factor(schoolcode) + factor(yeargc) 
  | 0 | woreda_code, 
  data = cohort3_3yearwindow
)

results_1 <- tidy_felm(model_1) %>% mutate(Model = "Born 1997-1998")
results_2 <- tidy_felm(model_2) %>% mutate(Model = "Born 1996-1999")
results_3 <- tidy_felm(model_3) %>% mutate(Model = "Born 1995-2000")

all_results <- bind_rows(results_1, results_2, results_3)

# ---------------------------- Formatting Results --------------------------------
all_results$stars <- ifelse(
  all_results$p.value < 0.01, "***", 
  ifelse(all_results$p.value < 0.05, "**", 
         ifelse(all_results$p.value < 0.1, "*", "")))

all_results$coeff_se <- paste0(
  round(all_results$estimate, 3), all_results$stars, "\n(",
  round(all_results$std.error, 3), ")"
)

desired_coeffs <- c("Post", "Post:Simulated_OrgIodine_ZScore", "female", 
                    "natural_stream", "(Intercept)")
all_results$term <- factor(all_results$term, levels = desired_coeffs)

# ---------------------------- Summary Statistics --------------------------------
cohort1_pre_war <- dplyr::filter(cohort1_1yearwindow, birth_yearGC == 1997)
cohort2_pre_war <- dplyr::filter(cohort2_2yearwindow, birth_yearGC %in% 1996:1997)
cohort3_pre_war <- dplyr::filter(cohort3_3yearwindow, birth_yearGC %in% 1995:1997)

sd_pre_war_total_3Sub <- sd(
  c(cohort1_pre_war$total_3Sub, cohort2_pre_war$total_3Sub, cohort3_pre_war$total_3Sub),
  na.rm = TRUE
)

median_total_3Sub <- c(
  median(cohort1_1yearwindow$total_3Sub, na.rm = TRUE),
  median(cohort2_2yearwindow$total_3Sub, na.rm = TRUE),
  median(cohort3_3yearwindow$total_3Sub, na.rm = TRUE)
)

mean_total_3Sub <- c(
  mean(cohort1_1yearwindow$total_3Sub, na.rm = TRUE),
  mean(cohort2_2yearwindow$total_3Sub, na.rm = TRUE),
  mean(cohort3_3yearwindow$total_3Sub, na.rm = TRUE)
)

mean_Simulated_OrgIodine <- c(
  mean(cohort1_1yearwindow$Simulated_OrgIodine, na.rm = TRUE),
  mean(cohort2_2yearwindow$Simulated_OrgIodine, na.rm = TRUE),
  mean(cohort3_3yearwindow$Simulated_OrgIodine, na.rm = TRUE)
)

sd_Simulated_OrgIodine <- c(
  sd(cohort1_1yearwindow$Simulated_OrgIodine, na.rm = TRUE),
  sd(cohort2_2yearwindow$Simulated_OrgIodine, na.rm = TRUE),
  sd(cohort3_3yearwindow$Simulated_OrgIodine, na.rm = TRUE)
)

median_total_3Sub_row <- tibble(
  term = "Median total_3Sub (SD pre-war)",
  `Born 1997-1998` = paste0(round(median_total_3Sub[1], 3), " (", round(sd_pre_war_total_3Sub, 3), ")"),
  `Born 1996-1999` = paste0(round(median_total_3Sub[2], 3), " (", round(sd_pre_war_total_3Sub, 3), ")"),
  `Born 1995-2000` = paste0(round(median_total_3Sub[3], 3), " (", round(sd_pre_war_total_3Sub, 3), ")")
)

mean_total_3Sub_row <- tibble(
  term = "Average total_3Sub (SD pre-war)",
  `Born 1997-1998` = paste0(round(mean_total_3Sub[1], 3), " (", round(sd_pre_war_total_3Sub, 3), ")"),
  `Born 1996-1999` = paste0(round(mean_total_3Sub[2], 3), " (", round(sd_pre_war_total_3Sub, 3), ")"),
  `Born 1995-2000` = paste0(round(mean_total_3Sub[3], 3), " (", round(sd_pre_war_total_3Sub, 3), ")")
)

mean_Simulated_OrgIodine_row <- tibble(
  term = "Average Simulated_OrgIodine",
  `Born 1997-1998` = paste0(round(mean_Simulated_OrgIodine[1], 3), " (", round(sd_Simulated_OrgIodine[1], 3), ")"),
  `Born 1996-1999` = paste0(round(mean_Simulated_OrgIodine[2], 3), " (", round(sd_Simulated_OrgIodine[2], 3), ")"),
  `Born 1995-2000` = paste0(round(mean_Simulated_OrgIodine[3], 3), " (", round(sd_Simulated_OrgIodine[3], 3), ")")
)

formatted_table <- all_results %>%
  dplyr::filter(term %in% desired_coeffs) %>%
  select(term, Model, coeff_se) %>%
  spread(Model, coeff_se) %>%
  arrange(match(term, desired_coeffs)) %>%
  select(term, `Born 1997-1998`, `Born 1996-1999`, `Born 1995-2000`)

formatted_table <- dplyr::bind_rows(
  dplyr::filter(formatted_table, term != "(Intercept)"),
  dplyr::filter(formatted_table, term == "(Intercept)"),
  median_total_3Sub_row,
  mean_total_3Sub_row,
  mean_Simulated_OrgIodine_row
)

# ---------------------------- Model Statistics ---------------------------------
extract_model_stats <- function(model, cohort_data, model_name) {
  n <- nrow(cohort_data)
  r2 <- summary(model)$r.squared
  tibble(
    term = c("N", "R2"),
    Model = model_name,
    Value = as.character(c(n, round(r2, 3)))
  )
}

stats_1 <- extract_model_stats(model_1, cohort1_1yearwindow, "Born 1997-1998")
stats_2 <- extract_model_stats(model_2, cohort2_2yearwindow, "Born 1996-1999")
stats_3 <- extract_model_stats(model_3, cohort3_3yearwindow, "Born 1995-2000")

combined_stats <- bind_rows(stats_1, stats_2, stats_3) %>%
  spread(Model, Value)

formatted_table <- dplyr::bind_rows(
  dplyr::filter(formatted_table, term != "N" & term != "R2"),
  combined_stats
)

# ---------------------------- Save the Table -----------------------------------
# Save the formatted table to an Excel file
output_file <- file.path(output_dir, "Standard_doseResponseDiD_TotalScore_Table1.xlsx")
write.xlsx(formatted_table, output_file, overwrite = TRUE)

# Print the final table
print(formatted_table)

###############################################################################
# Table Title:
# Table 1. Dose-response of secondary-school exam scores to local soil iodine 
#           for children born after loss of access to iodized salt
#
# Notes: 
# The first 3 columns show results for the rural study sample (treatment effect),
# where food is mostly locally sourced. The next 3 columns (in the actual final 
# manuscript table) would show results for the urban districts around Addis Ababa 
# (placebo test), where food largely comes from elsewhere, thereby serving as a 
# placebo. 
#
# To generate the placebo results, repeat the steps for data filtering and 
# regression on the urban sample (e.g., TestScore_data filtered to urban_nearAddis == 1).
#
# For Supplementary Table S5, the same steps are used but the iodine variable is replaced 
# with selenium (Simulated_SolubleSe_ZScore) to perform a falsification test. If the 
# iodine-related results were artifactual, similar effects would appear with selenium. 
# The absence of significant selenium effects supports the specificity and validity 
# of our iodine findings.
###############################################################################


###############################################################################
# Figure 3: Event-Study of Dose-Response by Year of Birth
#
# Manuscript: "Iodized salt has large impacts on child health and test scores in Ethiopia"
#
# Description:
# This section of the script implements an event-study design to show how the 
# relationship between local soil iodine (and selenium as a placebo) and children's
# exam scores evolves for each birth cohort relative to the reference year (1997).
#
# Key idea:
# - Instead of just comparing children born before and after 1998 (as in the main DID), 
#   we estimate separate coefficients for each birth year from 1994 to 2000.
# - 1997 is chosen as the reference year because it is the last year before the war-induced 
#   iodized salt disruption in 1998.
# - For each birth year, we estimate how a 1 SD increase in soil iodine relates to 
#   test scores, subtracting the reference year's coefficient so that all effects 
#   are shown relative to 1997.
#
# In this example, we use "total_3Sub_std" as the outcome, but for other outcomes 
# (like individual subject scores), you would only need to change the dependent 
# variable in the regression formulas.
###############################################################################

# Event-study DID models

# 1. Define the birth years of interest
year_range <- 1994:2000

# 2. Create dummy variables for each birth year
#    For each year in year_range, we create a binary indicator that equals 1 if the child 
#    was born in that year, and 0 otherwise.
for (year in year_range) {
  dummy_name <- paste0("birth_year_", year)
  TestScore_data[[dummy_name]] <- ifelse(TestScore_data$birth_yearGC == year, 1, 0)
}

# 3. Initialize a data frame to store regression results
results <- data.frame(
  birth_year = integer(),
  iodine_coef = double(),
  selenium_coef = double(),
  iodine_se = double(),
  selenium_se = double()
)

# 4. Run regressions for each birth year and extract coefficients
#    For each year, we run two regressions:
#    - One with soil iodine (Simulated_OrgIodine_ZScore) interacted with the birth-year dummy
#      to measure the dose-response for iodine.
#    - One with soil selenium (Simulated_SolubleSe_ZScore) for a placebo check.
#
#    We include controls (female, natural_stream, fixed effects for school and year)
#    and cluster standard errors at the district (woreda_code) level.

for (year in year_range) {
  dummy_name <- paste0("birth_year_", year)
  
  # Iodine Model
  # The formula includes an interaction between the birth-year dummy and Simulated_OrgIodine_ZScore.
  # total_3Sub_std ~ birth_year_dummy * iodine + controls + fixed effects
  formula_string_iodine <- sprintf(
    "total_3Sub_std ~ %s * Simulated_OrgIodine_ZScore + female + natural_stream + factor(schoolcode) + factor(yeargc)",
    dummy_name
  )
  iodine_model <- lm(as.formula(formula_string_iodine), data = TestScore_data)
  iodine_se <- sqrt(diag(vcovHC(iodine_model, type = "HC1", cluster = "group", group = TestScore_data$woreda_code)))
  
  # Selenium Model (Placebo)
  # Same structure but with Simulated_SolubleSe_ZScore instead of iodine.
  formula_string_selenium <- sprintf(
    "total_3Sub_std ~ %s * Simulated_SolubleSe_ZScore + female + natural_stream + factor(schoolcode) + factor(yeargc)",
    dummy_name
  )
  selenium_model <- lm(as.formula(formula_string_selenium), data = TestScore_data)
  selenium_se <- sqrt(diag(vcovHC(selenium_model, type = "HC1", cluster = "group", group = TestScore_data$woreda_code)))
  
  # Extract the interaction coefficients and standard errors
  results <- rbind(results, data.frame(
    birth_year = year,
    iodine_coef = coefficients(iodine_model)[paste0(dummy_name, ":Simulated_OrgIodine_ZScore")],
    selenium_coef = coefficients(selenium_model)[paste0(dummy_name, ":Simulated_SolubleSe_ZScore")],
    iodine_se = iodine_se[paste0(dummy_name, ":Simulated_OrgIodine_ZScore")],
    selenium_se = selenium_se[paste0(dummy_name, ":Simulated_SolubleSe_ZScore")]
  ))
}

# 5. Adjust coefficients relative to the 1997 reference year
#    We subtract the 1997 coefficient from all other years so that 1997 becomes the zero point.
coef_1997_iodine <- results[results$birth_year == 1997, "iodine_coef"]
coef_1997_selenium <- results[results$birth_year == 1997, "selenium_coef"]

results$iodine_coef <- results$iodine_coef - coef_1997_iodine
results$selenium_coef <- results$selenium_coef - coef_1997_selenium

# 6. Reshape data for plotting
results_long <- results %>%
  pivot_longer(cols = c(iodine_coef, selenium_coef, iodine_se, selenium_se),
               names_to = "type",
               values_to = "value") %>%
  mutate(measure = ifelse(grepl("coef", type), "coef", "se"),
         treatment_type = case_when(
           grepl("iodine", type) ~ "iodine",
           grepl("selenium", type) ~ "selenium"
         )) %>%
  select(birth_year, treatment_type, measure, value)

results_merged <- merge(
  results_long[results_long$measure == "coef", ],
  results_long[results_long$measure == "se", ],
  by = c("birth_year", "treatment_type")
)

results_merged$treatment_type <- factor(results_merged$treatment_type, levels = c("iodine", "selenium"))

# 7. Plot the event-study coefficients
#    Each point represents the coefficient for a given birth year (relative to 1997),
#    with error bars showing 95% confidence intervals. Red = iodine, Blue = selenium.

plot <- ggplot(data = results_merged, 
               aes(x = birth_year, y = value.x, color = treatment_type)) +
  geom_point(size = 4) +
  geom_line(aes(group = treatment_type), linetype = "dashed", size = 1) +
  geom_errorbar(aes(ymin = value.x - 1.96 * value.y, ymax = value.x + 1.96 * value.y), 
                width = 0.2, size = 1) +
  labs(y = "Test score (SD), Difference from children born in 1997",
       x = "Year of Birth") +
  geom_vline(xintercept = 1997, linetype = "dashed", color = "grey", size = 1) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 1) +
  scale_x_continuous(breaks = 1994:2000) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey92", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        legend.position = c(0.05, 0.95),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_color_manual(values = c("iodine" = "red", "selenium" = "blue"),
                     labels = c("Local Soil Iodine Concentration (SD)",
                                "Local Soil Selenium Concentration (SD)"))

print(plot)

# Saving the figure if needed
# ggsave(filename = "EventStudy_Figure3.png", plot = plot, width = 14, height = 8, units = "in")
###############################################################################

# Figure 3: Event study of dose-response in exam scores to 1 SD higher local soil nutrients by year of birth
# Note:
# The figure displays coefficients from the event-study regressions, showing how the 
# impact of a 1 SD increase in soil iodine concentration on test scores differs by birth year, 
# relative to the reference year 1997. The red points and bars are for iodine (the actual 
# treatment), and the blue points and bars are for selenium (the placebo test). If the 
# placebo results remain around zero, it supports the specificity of the iodine disruption effect.

# The event-study specification allows us to check for pre-trends (birth years before 1997) 
# and to see whether there's a distinct break starting in 1998, the year after iodized salt 
# was disrupted. This approach confirms that the dose-response effect of soil iodine on test 
# performance emerges post-1997, consistent with reduced protection from iodized salt. 

# For other outcomes (e.g., different subject test scores), the same code can be applied, 
# simply replacing 'total_3Sub_std' with the chosen dependent variable in the regression formula.
###############################################################################

###############################################################################
# Dose-response DID estimator by subject test score and birth-year cohort
#
# Figures S2 and S3 in the supplementary information:
#
# - Figure S2: Coefficients on 'Post' term by subject and cohort window (Table S2).
#   Shows how test scores differ for those born after vs. before iodized salt loss,
#   independent of soil iodine levels.
#
# - Figure S3: Coefficients on 'Post × Soil Iodine' term by subject and cohort window (Table S3).
#   Shows the dose-response to local soil iodine for those born after vs. before the cutoff.
#
# Steps:
# 1. Define test scores and labels.
# 2. Use previously defined cohorts: 1-year, 2-year, 3-year windows around 1998.
# 3. Run regressions for each subject and cohort, extracting 'Post' and 'Post×Iodine' coefficients.
# 4. Compute confidence intervals and arrange data for plotting.
# 5. Plot the results for each term, displaying estimates and CIs by subject and cohort.
###############################################################################

# (1) Define standardized test score variables and create formatted labels
test_scores <- c("english_std", "aptitude_std", "civics_std", "naturalmath_std", 
                 "socialmath_std", "physics_std", "biology_std", "chemistry_std", 
                 "geography_std", "history_std", "economics_std")

formatted_labels <- gsub("_std", "", test_scores)
formatted_labels <- sapply(strsplit(formatted_labels, "_"), function(x) { 
  paste0(toupper(substring(x,1,1)), substring(x,2))
})

# The data subsets (cohort1_1yearwindow, cohort2_2yearwindow, cohort3_3yearwindow)
# and the 'run_did_regression()' function are assumed to be already defined elsewhere.
# These subsets compare cohorts born:
# - 1-year window: 1997 vs 1998
# - 2-year window: 1996-97 vs 1998-99
# - 3-year window: 1995-97 vs 1998-2000

# (2) Initialize data frames to store 'Post' and 'Post×Iodine' results
results_post <- data.frame()
results_interaction <- data.frame()

# (3) Loop over each test score and each cohort set. For each combination:
#     - Extract the 'Post' coefficient (Figure S2)
#     - Extract the 'Post×SoilIodine' coefficient (Figure S3)
for (score in test_scores) {
  for (cohort in 1:3) {
    # Select the relevant cohort data (already created)
    cohort_data <- get(paste0("cohort", cohort, "_", cohort, "yearwindow"))
    
    # Run DID regression and store 'Post' coefficient and SE
    stats_post <- run_did_regression(cohort_data, score, "Post")
    results_post <- rbind(results_post, cbind(
      TestScore = score, 
      Cohort = paste0("Cohort ", cohort), 
      Coef = stats_post[1], 
      SE = stats_post[2], 
      N = stats_post[3]
    ))
    
    # Run DID regression and store 'Post×Simulated_OrgIodine_ZScore' coefficient and SE
    stats_interaction <- run_did_regression(cohort_data, score, "Post:Simulated_OrgIodine_ZScore")
    results_interaction <- rbind(results_interaction, cbind(
      TestScore = score, 
      Cohort = paste0("Cohort ", cohort), 
      Coef = stats_interaction[1], 
      SE = stats_interaction[2], 
      N = stats_interaction[3]
    ))
  }
}

# (4) Convert results to data frames and ensure numeric format
results_post_df <- as.data.frame(results_post)
results_post_df$Coef <- as.numeric(as.character(results_post_df$Coef))
results_post_df$SE <- as.numeric(as.character(results_post_df$SE))

results_interaction_df <- as.data.frame(results_interaction)
results_interaction_df$Coef <- as.numeric(as.character(results_interaction_df$Coef))
results_interaction_df$SE <- as.numeric(as.character(results_interaction_df$SE))

# Compute 95% confidence intervals (Coef ± 1.96*SE)
results_post_df$Lower <- results_post_df$Coef - 1.96 * results_post_df$SE
results_post_df$Upper <- results_post_df$Coef + 1.96 * results_post_df$SE

results_interaction_df$Lower <- results_interaction_df$Coef - 1.96 * results_interaction_df$SE
results_interaction_df$Upper <- results_interaction_df$Coef + 1.96 * results_interaction_df$SE

# Assign numeric positions along x-axis for each subject, 
# and slightly shift positions by cohort so all three cohorts appear side-by-side
central_positions <- seq_along(test_scores)

results_post_df$Position <- with(results_post_df, 
                                 central_positions[match(TestScore, test_scores)] + 
                                   (as.numeric(as.factor(Cohort)) - 2) / 4)

results_interaction_df$Position <- with(results_interaction_df, 
                                        central_positions[match(TestScore, test_scores)] + 
                                          (as.numeric(as.factor(Cohort)) - 2) / 4)

# Relabel Cohorts for the legend to show a clearer interpretation
cohort_levels <- c("Cohort 1", "Cohort 2", "Cohort 3")
cohort_labels <- c("Born in 1998 (vs 1997)",
                   "Born 1998-99 (vs. 1996-97)",
                   "Born 1998-2000 (vs 1995-97)")

results_post_df$Cohort <- factor(results_post_df$Cohort, levels = cohort_levels, labels = cohort_labels)
results_interaction_df$Cohort <- factor(results_interaction_df$Cohort, levels = cohort_levels, labels = cohort_labels)

library(ggplot2)

# (5) Plotting Figure S2 (Post coefficients):
# Shows differences in test scores for those born just after (vs. before) iodized salt loss,
# independent of soil iodine. Different colors = different birth-year windows.
plot_post <- ggplot(results_post_df, aes(x = Position, y = Coef, color = Cohort)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = central_positions, labels = formatted_labels) +
  labs(
    x = "",
    y = "Test score (SD), difference for post-salt-loss cohorts",
    color = "Cohorts"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white")
  )

# Plotting Figure S3 (Post × Soil Iodine interaction):
# Shows how 1 SD higher soil iodine affects post-cutoff cohorts vs. pre-cutoff cohorts.
plot_interaction <- ggplot(results_interaction_df, aes(x = Position, y = Coef, color = Cohort)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = central_positions, labels = formatted_labels) +
  labs(
    x = "",
    y = "Test score (SD) per 1 SD higher soil iodine (Post vs. Pre)",
    color = "Cohorts"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white")
  )

# Print the plots to the R console
print(plot_post)
print(plot_interaction)

# If needed, save the figures:
# ggsave(filename = file.path(output_dir, "FigureS2_Post_bySubjectCohort.png"), plot = plot_post, width = 16, height = 10, units = "in")
# ggsave(filename = file.path(output_dir, "FigureS3_Interaction_bySubjectCohort.png"), plot = plot_interaction, width = 18, height = 10, units = "in")