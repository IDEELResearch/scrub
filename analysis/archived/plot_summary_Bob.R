library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)

# Load rds file
data <- readRDS("analysis/data-derived/final_data.rds")

# Extract year from the date
data <- data %>%
  mutate(year = year(as.Date(collection_day, format = "%Y-%m-%d")))

# Count unique study_IDs per year and country
study_counts <- data %>%
  distinct(study_ID, year, database) %>%
  count(year, database) 
write.csv(study_counts, "analysis/Bob_plots_data/study_counts.csv")
# Count unique survey_IDs per year and country
survey_counts <- data %>%
  distinct(survey_ID, year, database) %>%
  count(year, database)
write.csv(survey_counts, "analysis/Bob_plots_data/survey_counts.csv")
# Count sample size per year and country
sample_counts <- data %>%
  mutate(total_num = as.numeric(total_num)) %>%                    # Ensure numeric
  group_by(survey_ID, year, database) %>%                            # Group by survey & country
  summarise(max_sample = max(total_num, na.rm = TRUE), .groups = "drop") %>%  # Get max per survey
  group_by(year, database) %>%                                       # Now group by country
  summarise(sample_size = sum(max_sample, na.rm = TRUE), .groups = "drop")    # Sum across surveys
write.csv(sample_counts, "analysis/Bob_plots_data/sample_counts.csv")

# Plot Barplot of study IDs per year  -------------------------------------
study_count_barplot <- ggplot(study_counts, aes(x = year, y = n, fill = database)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_x_continuous(breaks = seq(min(study_counts$year), max(study_counts$year), by = 5)) +
  labs(x = "Year", y = "Number of Study IDs", fill = "Database") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()
ggsave("analysis/Bob_plots_data/study_count_barplot.png", study_count_barplot)

# Plot Barplot of survey IDs per year  ------------------------------------
survey_count_barplot <- ggplot(survey_counts, aes(x = year, y = n, fill = database)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_x_continuous(breaks = seq(min(survey_counts$year), max(survey_counts$year), by = 5)) +
  labs(x = "Year", y = "Number of Survey IDs", fill = "Database") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()
ggsave("analysis/Bob_plots_data/survey_count_barplot.png", survey_count_barplot)

# Plot Barplot of sample size per year for each country -------------------
sample_size_barplot <- ggplot(sample_counts, aes(x = year, y = sample_size, fill = database)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_x_continuous(breaks = seq(min(sample_counts$year), max(sample_counts$year), by = 5)) +
  labs(
    x = "Year",
    y = "Sample Size",
    fill = "Database"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()
ggsave("analysis/Bob_plots_data/sample_size_barplot.png", sample_size_barplot)
