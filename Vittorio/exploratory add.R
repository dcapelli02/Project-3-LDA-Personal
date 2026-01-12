library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

dropout_check <- alz_long %>%
  group_by(id) %>%
  arrange(time) %>%
  summarise(
    first_missing = ifelse(any(!observed),
                           min(time[!observed]),
                           NA),
    non_monotone = ifelse(
      is.na(first_missing),
      FALSE,
      any(observed & time > first_missing)
    ),
    .groups = "drop"
  )
table(dropout_check$non_monotone)


# Pivot wide -> long
cdrsb_bin_long <- alz %>%
  select(matches("cdrsb_bin\\d+$")) %>%
  mutate(id = row_number()) %>%  # aggiunge un id
  pivot_longer(
    cols = starts_with("cdrsb_bin"),
    names_to = "visit",
    values_to = "value"
  ) %>%
  mutate(
    visit = as.numeric(gsub("cdrsb_bin", "", visit))
  )

# Trovo il primo dropout (prima NA o ultima osservazione disponibile)
# Qui assumiamo che NA indica dropout
cdrsb_bin_long <- cdrsb_bin_long %>%
  group_by(id) %>%
  mutate(dropout_time = ifelse(all(!is.na(value)), max(visit), min(visit[is.na(value)]))) %>%
  ungroup()

# Calcolo proporzione di value = 1 per visit e gruppo dropout
cdrsb_prop_group <- cdrsb_bin_long %>%
  group_by(dropout_time, visit) %>%
  summarise(prop_positive = mean(value == 1, na.rm = TRUE), .groups = "drop")

# Grafico
ggplot(cdrsb_prop_group, aes(x = visit, y = prop_positive, color = factor(dropout_time))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Visit",
    y = "Proportion with CDRSB > 10",
    color = "Dropout Time",
    title = "Evolution of Binary CDRSB Outcome by Dropout Time"
  ) +
  theme_minimal()


# Calcolo proporzioni
alz_long_lag <- alz_long %>%
  group_by(id) %>%
  mutate(
    is_missing_next = lead(is.na(cdrsb)) 
  ) %>%
  filter(time < max(time)) 

cdrsb_bin_prop <- alz_long_lag %>%
  group_by(time, is_missing_next) %>%
  summarise(prop_positive = mean(cdrsb_bin == 1, na.rm = TRUE), .groups = "drop")

# Grafico
ggplot(cdrsb_bin_prop, aes(x = as.factor(time), y = prop_positive, color = is_missing_next, group = is_missing_next)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = c("skyblue", "salmon")) +
  labs(
    x = "Year",
    y = "Proportion with CDRSB > 10",
    color = "Missing at t+1?"
  ) +
  theme_minimal()
