
# Script: Three Parks Savanna Fire-effects Plot Network: The effect of fire
#         treatment on soil properties
# Author: Luke Florence
# Date:   22 January 2024

# Note: Models have previously undergone fit tests and have been assessed for 
# appropriate structure and transformations. Models fit on a Gaussian
# distribution are best fit with "block" as a random effect and Poisson
# distribution with "block/sample" as random effects. Nevertheless, because I
# have undertaken a novel approach to SEM model selection (see notes in section
# (9) Structural equation models) I include the model selection process within
# the corresponding section.

# Required packages and functions:
require(multcompView)
require(lme4)
require(emmeans)
require(ggeffects)
require(performance)
require(parameters)
require(DHARMa)
require(openxlsx)
require(patchwork)
require(viridis)
require(ggpubr)
require(piecewiseSEM)
require(tidyverse)
source("code/functions.R")

# Read in the data:
data <- read.csv("data/data.csv", stringsAsFactors = TRUE) %>%
  # Order levels
  mutate(treatment = factor(treatment, levels = c(
    "U", "E5", "E4", "E3", "E2","L2", "E1"))
  ) %>%
  glimpse(.)

# Ensure 'treatment' is a factor
data$treatment <- factor(data$treatment)

# Check the unique levels of 'treatment'
unique_levels <- levels(data$treatment)

# Set sum contrasts
contrasts(data$treatment) <- contr.treatment(unique_levels)

### (1) Fungal gene abundance ##################################################

model_fungi <- glmer(
  fungi_abundance ~ treatment + (1 | block/sample_id),
  data = data,
  family = poisson
  )

# Model validation
simulateResiduals(model_fungi) %>%
  plot(.)

# Model coefficients
summary_fungi <- summarise_results(
  model_fungi, test_type = "z"
)

# Estimated marginal means
emmeans_fungi <- emmeans(
  model_fungi,
  list(pairwise ~ treatment),
  adjust ="tukey",
  type = "response",
  biased.adj = TRUE
  )

# Means for plots
estimate_fungi <- emmeans_fungi$`emmeans of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment, rate, SE,  df, lower_CI = asymp.LCL, upper_CI = asymp.UCL
  ) %>%
  mutate(
    across(c(-treatment, -df), ~round(., 0)),
    df = as.character(df)
  ) %>%
  print(.)

# Pairwise comparison
pairwise_fungi <- emmeans_fungi$`pairwise differences of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment = 1, ratio, SE,  df, "z_ratio" = z.ratio, "p_value" = p.value) %>%
  mutate(
    across(c(-treatment, df, -p_value), ~round(., 2)),
    p_value = round(p_value, 3),
    df = as.character(df)
  ) %>%
  print(.)

### (2) Bacterial gene abundance ###############################################

model_bacteria <- glmer(
  bacteria_abundance ~ treatment + (1 | block/sample_id),
  data = data,
  family = poisson
  )

# Model validation
simulateResiduals(model_bacteria) %>%
  plot(.)

# Model coefficients
summary_bacteria <- summarise_results(
  model_bacteria, test_type = "z"
  )

# Estimated marginal means
emmeans_bacteria <- emmeans(
  model_bacteria,
  list(pairwise ~ treatment),
  adjust ="tukey",
  type = "response",
  biased.adj = TRUE
  )

# Means for plots
estimate_bacteria <- emmeans_bacteria$`emmeans of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment, rate, SE, df, lower_CI = asymp.LCL, upper_CI = asymp.UCL
  ) %>%
  mutate(
    across(c(-treatment, -df), ~round(., 0)),
    df = as.character(df)
  ) %>%
  print(.)

# Pairwise comparison
pairwise_bacteria <- emmeans_bacteria$`pairwise differences of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment = 1, ratio, SE,  df, "z_ratio" = z.ratio, "p_value" = p.value
  ) %>%
  mutate(
    across(c(-treatment, df, -p_value), ~round(., 2)),
    p_value = round(p_value, 3),
    df = as.character(df)
  ) %>%
  print(.)

### (3) Organic carbon #########################################################

model_carbon <- lmer(
  organic_carbon ~ treatment + (1 | block),
  data = data
  )

# Model validation
simulateResiduals(model_carbon) %>%
  plot(.)

# Model coefficients
summary_carbon <- summarise_results(model_carbon, test_type = "t")

# Estimated marginal means
emmeans_carbon <- emmeans(
  model_carbon,
  list(pairwise ~ treatment),
  adjust ="tukey"
  )

# Means for plots
estimate_carbon <- emmeans_carbon$`emmeans of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment, mean = emmean, SE, df, lower_CI = lower.CL, upper_CI = upper.CL
  ) %>%
  mutate(
    across(c(-treatment, -df), ~round(., 2)),
    df = round(df, 0)
  ) %>%
  print(.)

# Pairwise comparison
pairwise_carbon <- emmeans_carbon$`pairwise differences of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment = 1, estimate, SE, df, "t_ratio" = t.ratio, "p_value" = p.value
  ) %>%
  mutate(
    across(c(-treatment, -df, -p_value), ~round(., 2)),
    p_value = round(p_value, 3),
    df = round(df, 0)
  ) %>%
  print(.)

### (4) Fungal to bacterial ratio ##############################################

model_ratio <- lmer(
  fungi_bacteria_ratio ~ treatment + (1 | block),
  data = data
  )

# Model validation
simulateResiduals(model_ratio) %>%
  plot(.)

# Model coefficients
summary_ratio <- summarise_results(model_ratio, test_type = "t")

# Estimated marginal means
emmeans_ratio <- emmeans(
  model_ratio,
  list(pairwise ~ treatment),
  adjust ="tukey"
  )

# Means for plots
estimate_ratio <- emmeans_ratio$`emmeans of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment, mean = emmean, SE, df, lower_CI = lower.CL, upper_CI = upper.CL
  ) %>%
  mutate(
    across(c(-treatment, -df), ~round(., 2)),
    df = round(df, 0)
  ) %>%
  print(.)

# Pairwise comparison
pairwise_ratio <- emmeans_ratio$`pairwise differences of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment = 1, estimate, SE, df, "t_ratio" = t.ratio, "p_value" = p.value
  ) %>%
  mutate(
    across(c(-treatment, -df, -p_value), ~round(., 2)),
    p_value = round(p_value, 3),
    df = round(df, 0)
  ) %>%
  print(.)

### (5) Nitrogen ###############################################################

# Here I had to log transform nitrogen and drop a single upper outlier to get 
# the model to fit

model_nitrogen <- data %>%
  filter(nitrogen < 0.6) %>%
  lmer(
    log(nitrogen) ~ treatment + (1 | block),
    data = .
    )

# Model validation
simulateResiduals(model_nitrogen) %>%
  plot(.)

# Model coefficients
summary_nitrogen <- summarise_results(
  model_nitrogen, test_type = "t"
  )

# Estimated marginal means
emmeans_nitrogen <- emmeans(
  model_nitrogen,
  list(pairwise ~ treatment),
  adjust ="tukey",
  type = "response",
  biased.adj = TRUE
  )

# Means for plots
estimate_nitrogen <- emmeans_nitrogen$`emmeans of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment, response, SE, df, lower_CI = lower.CL, upper_CI = upper.CL
  ) %>%
  mutate(
    across(c(-treatment, -df), ~round(., 2)),
    df = round(df, 0)
  ) %>%
  print(.)

# Pairwise comparison
pairwise_nitrogen <- emmeans_nitrogen$`pairwise differences of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment = 1, ratio, SE, df, "t_ratio" = t.ratio, "p_value" = p.value
  ) %>%
  mutate(
    across(c(-treatment, -df, -p_value), ~round(., 2)),
    p_value = round(p_value, 3),
    df = round(df, 0)
  ) %>%
  print(.)

### (6) Carbon-to-nitrogen ratio ###############################################

model_CN <- lmer(
  carbon_nitrogen_ratio ~ treatment + (1 | block),
  data = data
  )

# Model validation
simulateResiduals(model_CN) %>%
  plot(.)

# Model coefficients
summary_CN <- summarise_results(model_CN, test_type = "t")

# Estimated marginal means
emmeans_CN <- emmeans(
  model_CN,
  list(pairwise ~ treatment),
  adjust ="tukey"
  )

# Means for plots
estimate_CN <- emmeans_CN$`emmeans of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment, mean = emmean, SE, df, lower_CI = lower.CL, upper_CI = upper.CL
  ) %>%
  mutate(
    across(c(-treatment, -df), ~round(., 2)),
    df = round(df, 0)
  ) %>%
  print(.)

# Pairwise comparison
pairwise_CN <- emmeans_CN$`pairwise differences of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment = 1, mean = estimate, SE, df, "t_ratio" = t.ratio,
    "p_value" = p.value
  ) %>%
  mutate(
    across(c(-treatment, -df, -p_value), ~round(., 2)),
    p_value = round(p_value, 3),
    df = round(df, 0)
  ) %>%
  print(.)

# Recalculate means for plots: 
# For some reason the confidence intervals are not being calculated correctly
# using the emmeans function directly, so I recalculate emmeans here using
# ggeffects
estimate_CN <- ggeffects::ggemmeans(model_CN, terms = "treatment") %>%
  as_tibble(.) %>%
  select(
    treatment = x, mean = predicted, SE = std.error, lower_CI = conf.low,
    upper_CI = conf.high
  ) %>%
  mutate(
    across(c(-treatment), ~round(., 2)),
    df = estimate_CN$df
  ) %>%
  select(
    treatment, mean, SE, df, lower_CI, upper_CI
  ) %>%
  print(.)

### (7) pH #####################################################################

model_pH <- lmer(
  pH ~ treatment + (1 | block),
  data = data
  )

# Model validation
simulateResiduals(model_pH) %>%
  plot(.)

# Model coefficients
summary_pH <- summarise_results(model_pH, test_type = "t")

# Estimated marginal means
emmeans_pH <- emmeans(
  model_pH,
  list(pairwise ~ treatment),
  adjust ="tukey"
  )

# Means for plots
estimate_pH <- emmeans_pH$`emmeans of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment, mean = emmean, SE, df, lower_CI = lower.CL, upper_CI = upper.CL
  ) %>%
  mutate(
    across(c(-treatment), ~round(., 2))
  ) %>%
  print(.)

# Pairwise comparison
pairwise_pH <- emmeans_pH$`pairwise differences of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment = 1, mean = estimate, SE, df, "t_ratio" = t.ratio,
    "p_value" = p.value
  ) %>%
  mutate(
    across(c(-treatment, -df, -p_value), ~round(., 2)),
    p_value = round(p_value, 3),
    df = round(df, 0)
  ) %>%
  print(.)

### (8) Total grass cover ####################################################

model_total_grass <- lmer(
  total_grass_cover ~ treatment + (1 | block),
  data = data %>% drop_na()
)

# Model validation
simulateResiduals(model_total_grass) %>%
  plot(.)

# Model coefficients
summary_total_grass <- summarise_results(
  model_total_grass, test_type = "t"
)

# Estimated marginal means
emmeans_total_grass <- emmeans(
  model_total_grass,
  list(pairwise ~ treatment),
  adjust ="tukey",
  type = "response",
  biased.adj = TRUE
)

# Means for plots
estimate_total_grass <- emmeans_total_grass$`emmeans of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment, mean = emmean, SE,  df, lower_CI = lower.CL, upper_CI = upper.CL
  ) %>%
  mutate(
    across(c(-treatment, -df), ~round(., 0)),
    df = as.character(df)
  ) %>%
  print(.)

# Pairwise comparison
pairwise_total_grass <- emmeans_total_grass$`pairwise differences of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment = 1, estimate, SE,  df, "t_ratio" = t.ratio, "p_value" = p.value) %>%
  mutate(
    across(c(-treatment, df, -p_value), ~round(., 2)),
    p_value = round(p_value, 3),
    df = as.character(df)
  ) %>%
  print(.)

### (9) Modern grass cover ####################################################

model_modern_grass <- lmer(
  sqrt(modern_grass_cover) ~ treatment + (1 | block),
  data = data
)

# Model validation
simulateResiduals(model_modern_grass) %>%
  plot(.)

# Model coefficients
summary_modern_grass <- summarise_results(
  model_modern_grass, test_type = "t"
)

# Estimated marginal means
emmeans_modern_grass <- emmeans(
  model_modern_grass,
  list(pairwise ~ treatment),
  adjust ="tukey",
  type = "response",
  biased.adj = TRUE
)

# Means for plots
estimate_modern_grass <- emmeans_modern_grass$`emmeans of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment, mean = response, SE,  df, lower_CI = lower.CL, upper_CI = upper.CL
  ) %>%
  mutate(
    across(c(-treatment, -df), ~round(., 0)),
    df = as.character(df)
  ) %>%
  print(.)

# Pairwise comparison
pairwise_modern_grass <- emmeans_modern_grass$`pairwise differences of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment = 1, estimate, SE,  df, "t_ratio" = t.ratio, "p_value" = p.value) %>%
  mutate(
    across(c(-treatment, df, -p_value), ~round(., 2)),
    p_value = round(p_value, 3),
    df = as.character(df)
  ) %>%
  print(.)

### (10) Ancient grass cover ####################################################

model_ancient_grass <- lmer(
  sqrt(ancient_grass_cover) ~ treatment + (1 | block),
  data = data
)

# Model validation
simulateResiduals(model_ancient_grass) %>%
  plot(.)

# Model coefficients
summary_ancient_grass <- summarise_results(
  model_ancient_grass, test_type = "t"
)

# Estimated marginal means
emmeans_ancient_grass <- emmeans(
  model_ancient_grass,
  list(pairwise ~ treatment),
  adjust ="tukey",
  type = "response",
  biased.adj = TRUE
)

# Means for plots
estimate_ancient_grass <- emmeans_ancient_grass$`emmeans of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment, mean = response, SE,  df, lower_CI = lower.CL, upper_CI = upper.CL
  ) %>%
  mutate(
    across(c(-treatment, -df), ~round(., 0)),
    df = as.character(df)
  ) %>%
  print(.)

# Pairwise comparison
pairwise_ancient_grass <- emmeans_ancient_grass$`pairwise differences of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment = 1, estimate, SE,  df, "t_ratio" = t.ratio, "p_value" = p.value) %>%
  mutate(
    across(c(-treatment, df, -p_value), ~round(., 2)),
    p_value = round(p_value, 3),
    df = as.character(df)
  ) %>%
  print(.)

### (11) Litter cover ##########################################################

model_litter <- lmer(
  litter_cover ~ treatment + (1 | block),
  data = data
)

# Model validation
simulateResiduals(model_litter) %>%
  plot(.)

# Model coefficients
summary_litter <- summarise_results(
  model_litter, test_type = "t"
)

# Estimated marginal means
emmeans_litter <- emmeans(
  model_litter,
  list(pairwise ~ treatment),
  adjust ="tukey",
  type = "response",
  biased.adj = TRUE
)

# Means for plots
estimate_litter <- emmeans_litter$`emmeans of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment, mean = emmean, SE,  df, lower_CI = lower.CL, upper_CI = upper.CL
  ) %>%
  mutate(
    across(c(-treatment, -df), ~round(., 0)),
    df = as.character(df)
  ) %>%
  print(.)

# Pairwise comparison
pairwise_litter <- emmeans_litter$`pairwise differences of treatment` %>%
  as_tibble(.) %>%
  select(
    treatment = 1, estimate, SE,  df, "t_ratio" = t.ratio, "p_value" = p.value) %>%
  mutate(
    across(c(-treatment, df, -p_value), ~round(., 2)),
    p_value = round(p_value, 3),
    df = as.character(df)
  ) %>%
  print(.)

### (12) Save test statistics ###################################################

# List of response variables
response_variables <- c(
  "fungi", "bacteria", "carbon","ratio", "nitrogen", "CN", "pH", "total_grass",
  "modern_grass", "ancient_grass", "litter"
)

save_results_to_excel(response_variables, "output/test_statistics.xlsx")

### (8) Main plots #############################################################

# Palette
pal <- rev(viridis::magma(6))

#### (8a) Fungi plots ####

# I'm not sure how to overlay observed values in ggplot without this weird
# workaround, where all metrics called from the main data frame need to be
# present in subsequent layers. I therefore make false metrics that are 
# repetition of the value I want to overlay, that is the observed values in the
# raw data
observed_values_fungi <- data %>%
  select(
    treatment, rate = fungi_abundance
    ) %>%
  mutate(
    rate = rate / 1000000,
    lower_CI = rate,
    upper_CI = rate,
    SE = rate
  )

# Observed values and marginal means plot: I divide abundance by a million to
# make the results to increase interretability
means_fungi <- ggplot(
  data = estimate_fungi,
  aes(x = treatment,
      ymin = lower_CI / 1000000,
      lower = (rate - SE) / 1000000,
      middle = rate / 1000000,
      upper = (rate + SE) / 1000000,
      ymax = upper_CI / 1000000,
      fill = treatment)
  )  +
  geom_boxplot(stat = 'identity', alpha = 0.6) +
  geom_point(
    observed_values_fungi,
    mapping = aes(x = treatment, y = rate, fill = treatment),
    position = position_jitter(seed = 1986, width = 0.25),
    shape = 21, size = 1, stroke = 0.25
  ) +
  scale_fill_manual(
    values = pal
  ) +
  scale_y_continuous(
    limits = c(0, (max(observed_values_fungi$rate) + 0.01))
    ) +
  xlab(NULL) +
  ylab("Fungal Abundance") +
  MyTheme() +
  theme(
     axis.text.x = element_blank()
     ) +
  annotate(
    "text", x = -0.5, y = (max(observed_values_fungi$rate) + 0.01), 
    label = "b", size = 5, fontface = "bold"
    ) +
  coord_cartesian(
    # This allows me to annotate outside panel borders, but it also emboldens
    # the panel borders so include this line in all plots regardless of
    # annotation
    xlim = c(1, 6), clip = "off"
    )

# Effect size means
mean_effect_fungi <- parameters(
  model_fungi,
  effects = "fixed"
) %>%
  as_tibble() %>%
  select(
    Coefficient, CI_low, CI_high
  ) %>%
  mutate(
    treatment = c("U", "E5", "E3", "E2", "L2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "E2", "L2", "E1"))
  )

# Format pairwise comparisons data frame
pairwise_plot_fungi <- pairwise_fungi %>%
  select(treatment, p_value) %>%
  filter(p_value <= 0.05) %>%
  mutate(
    group1 = word(treatment, 1, sep = " / "),
    group2 = word(treatment, 2, sep = " / ")
  ) %>%
  select(-treatment)

# Create effect size plot
effects_fungi <- ggplot(
  mean_effect_fungi,
  aes(x = treatment, y = Coefficient)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_errorbar(
    mapping = aes(x = treatment, ymax = CI_high, ymin = CI_low),
    color = "black",
    width = 0,
    linewidth = 0.6
  ) +
  geom_point(
    shape = 16, size = 3,
    colour = alpha(c(
      "#FEAF77FF", "#F1605DFF", "#B63679FF",
      "#721F81FF", "#2D1160FF", "#000004FF"),
      0.7)
  ) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme(
    axis.text.x = element_blank()
  ) +
  xlab(NULL) +
  ylab("Mean Difference in Effect Size") +
  MyTheme() +
  coord_cartesian(
    xlim = c(1, 6), clip = "off"
    ) +
  # Add pairwise comparisons
  stat_pvalue_manual(
    data = pairwise_plot_fungi, 
    y.position = 0.70,
    step.increase = 0.0075,
    tip.length = 0.002,
    bracket.size = 0.2,
    label.size = 2.5,
    label = "p_value"
  )

#### (8b) Bacteria plots ####
observed_values_bacteria <- data %>%
  select(
    treatment, rate = bacteria_abundance
    ) %>%
  mutate(
    rate = rate / 1000000,
    lower_CI = rate,
    upper_CI = rate,
    SE = rate
  )

# Observed values and marginal means plot
means_bacteria <- ggplot(
  data = estimate_bacteria,
  aes(x = treatment,
      ymin = lower_CI / 1000000,
      lower = (rate - SE) / 1000000,
      middle = rate / 1000000,
      upper = (rate + SE) / 1000000,
      ymax = upper_CI / 1000000,
      fill = treatment)
  )  +
  geom_boxplot(stat = 'identity', alpha = 0.6) +
  geom_point(
    observed_values_bacteria,
    mapping = aes(x = treatment, y = rate, fill = treatment),
    position = position_jitter(seed = 86, width = 0.25),
    shape = 21, size = 1, stroke = 0.25
  ) +
  scale_fill_manual(
    values = pal
  ) +
  scale_y_continuous(limits = c(0, 15)) +
  xlab(NULL) +
  ylab("Bacterial Abundance") +
  MyTheme() +
  theme(
    axis.text.x = element_blank()
    ) + 
  annotate(
    "text", x = -0.5, y = 15, 
    label = "a", size = 5, fontface = "bold"
    ) +
  coord_cartesian(
    xlim = c(1, 6), clip = "off"
    )

# Effect size means
mean_effect_bacteria <- parameters(
  model_bacteria,
  effects = "fixed"
) %>%
  as_tibble() %>%
  select(
    Coefficient, CI_low, CI_high
  ) %>%
  mutate(
    treatment = c("U", "E5", "E3", "E2", "L2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "E2", "L2", "E1"))
  )

# Format pairwise comarisons data frame
pairwise_plot_bacteria <- pairwise_bacteria %>%
  select(treatment, p_value) %>%
  filter(p_value <= 0.05) %>%
  mutate(
    group1 = word(treatment, 1, sep = " / "),
    group2 = word(treatment, 2, sep = " / ")
  ) %>%
  select(-treatment)

# Create effect size plot
effects_bacteria <- ggplot(
  mean_effect_bacteria,
  aes(x = treatment, y = Coefficient)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_errorbar(
    mapping = aes(x = treatment, ymax = CI_high, ymin = CI_low),
    color = "black",
    width = 0,
    linewidth = 0.6
  ) +
  geom_point(
    shape = 16, size = 3,
    colour = alpha(c(
      "#FEAF77FF", "#F1605DFF", "#B63679FF",
      "#721F81FF", "#2D1160FF", "#000004FF"),
      0.7)
  ) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme(
    axis.text.x = element_blank()
  ) +
  xlab(NULL) +
  ylab("Mean Difference in Effect Size") +
  MyTheme() +
  coord_cartesian(
    xlim = c(1, 6), clip = "off"
    ) +
  # Add pairwise comparisons
  stat_pvalue_manual(
    data = pairwise_plot_bacteria, 
    y.position = 0.70,
    step.increase = 0.0075,
    tip.length = 0.002,
    bracket.size = 0.2,
    label.size = 2.5,
    label = "p_value"
  )

#### (8c) Carbon plots ####
observed_values_carbon <- data %>%
  select(
    treatment, mean = organic_carbon
    ) %>%
  mutate(
    lower_CI = mean,
    upper_CI = mean,
    SE = mean)

# Observed values and marginal means plot
means_carbon <- ggplot(
  data = estimate_carbon,
  aes(x = treatment,
      ymin = lower_CI,
      lower = mean - SE,
      middle = mean,
      upper = mean + SE,
      ymax = upper_CI,
      fill = treatment)
  )  +
  geom_boxplot(stat = 'identity', alpha = 0.6) +
  geom_point(
    observed_values_carbon,
    mapping = aes(x = treatment, y = mean, fill = treatment),
    position = position_jitter(seed = 1986, width = 0.25),
    shape = 21, size = 1, stroke = 0.25
  ) +
  scale_fill_manual(
    values = pal
  ) +
  scale_y_continuous(limits = c(0, 3)) +
  xlab("Fire Treatment") +
  ylab("Organic Carbon (%)") +
  MyTheme() +
  annotate(
    "text", x = -0.5, y = 3, 
    label = "d", size = 5, fontface = "bold"
    ) +
  coord_cartesian(
    xlim = c(1, 6), clip = "off"
    )

# Effect size means
mean_effect_carbon <- parameters(
  model_carbon,
  effects = "fixed"
) %>%
  as_tibble() %>%
  select(
    Coefficient, CI_low, CI_high
  ) %>%
  mutate(
    treatment = c("U", "E5", "E3", "E2", "L2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "E2", "L2", "E1"))
  )

# Create effect size plot
effects_carbon <- ggplot(
  mean_effect_carbon,
  aes(x = treatment, y = Coefficient, fill = treatment)) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_errorbar(
    mapping = aes(x = treatment, ymax = CI_high, ymin = CI_low),
    color = "black",
    width = 0,
    linewidth = 0.7
  ) +
  geom_point(
    shape = 16, size = 3,
    colour = alpha(c(
      "#FEAF77FF", "#F1605DFF", "#B63679FF",
      "#721F81FF", "#2D1160FF", "#000004FF"),
      0.7)
  ) +
  scale_y_continuous(limits = c(-1, 1)) +
  xlab("Fire Treatment") +
  ylab("Mean Difference in Effect Size") +
  MyTheme() +
  coord_cartesian(
    xlim = c(1, 6), clip = "off"
    )

#### (8d) Fungal to bacterial ratio  ####
observed_values_ratio <- data %>%
  select(
    treatment, mean = fungi_bacteria_ratio
  ) %>%
  mutate(
    lower_CI = mean,
    upper_CI = mean,
    SE = mean)

# Observed values and marginal means plot
means_ratio <- ggplot(
  data = estimate_ratio,
  aes(x = treatment,
      ymin = lower_CI,
      lower = mean - SE,
      middle = mean,
      upper = mean + SE,
      ymax = upper_CI,
      fill = treatment)
)  +
  geom_boxplot(stat = 'identity', alpha = 0.6) +
  geom_point(
    observed_values_ratio,
    mapping = aes(x = treatment, y = mean, fill = treatment),
    position = position_jitter(seed = 1986, width = 0.25),
    shape = 21, size = 1, stroke = 0.25
  ) +
  scale_fill_manual(
    values = pal
  ) +
  scale_y_continuous(limits = c(0, 3)) +
  xlab(NULL) +
  ylab("Fungal-to-Bacterial Ratio") +
  MyTheme() +
  theme(
    axis.text.x = element_blank()
  ) + 
  annotate(
    "text", x = -0.5, y = 3, 
    label = "c", size = 5, fontface = "bold"
  ) +
  coord_cartesian(
    xlim = c(1, 6), clip = "off"
  )

# Effect size means
mean_effect_ratio <- parameters(
  model_ratio,
  effects = "fixed"
) %>%
  as_tibble() %>%
  select(
    Coefficient, CI_low, CI_high
  ) %>%
  mutate(
    treatment = c("U", "E5", "E3", "E2", "L2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "E2", "L2", "E1"))
  )

# Create effect size plot
effects_ratio <- ggplot(
  mean_effect_ratio,
  aes(x = treatment, y = Coefficient, fill = treatment)) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_errorbar(
    mapping = aes(x = treatment, ymax = CI_high, ymin = CI_low),
    color = "black",
    width = 0,
    linewidth = 0.7
  ) +
  geom_point(
    shape = 16, size = 3,
    colour = alpha(c(
      "#FEAF77FF", "#F1605DFF", "#B63679FF",
      "#721F81FF", "#2D1160FF", "#000004FF"),
      0.7)
  ) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme(
    axis.text.x = element_blank()
  ) +
  xlab(NULL) +
  ylab("Mean Difference in Effect Size") +
  MyTheme() +
  coord_cartesian(
    xlim = c(1, 6), clip = "off"
  )


#### (8e) Join and save Figure 2 ####
main_plots <- 
  means_bacteria + effects_bacteria + means_fungi + effects_fungi +
  means_ratio + effects_ratio + means_carbon + effects_carbon + 
  plot_layout(ncol = 2)
ggsave("output/Figure1.jpeg", width = 4.75, height = 8.4)
ggsave("output/Figure1.tiff", width = 4.75, height = 8.4)
ggsave("output/Figure1.pdf", width = 4.75, height = 8.4)

#### (8f) Additional effects plots: Main text ####

# Will make effects plots here, which will be used to aid the interpretation of
# the SEMs below

# Carbon-to-nitrogen: Format mean effect size data frame
mean_effect_CN <- parameters(
  model_CN,
  effects = "fixed"
) %>%
  as_tibble() %>%
  select(
    Coefficient, CI_low, CI_high
  ) %>%
  mutate(
    treatment = c("U", "E5", "E3", "E2", "L2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "E2", "L2", "E1"))
  )

# Carbon-to-nitrogen: Format pairwise comarisons data frame
pairwise_plot_CN <- pairwise_CN %>%
  select(treatment, p_value) %>%
  filter(p_value <= 0.05) %>%
  mutate(
    group1 = word(treatment, 1, sep = " - "),
    group2 = word(treatment, 2, sep = " - "),
    p_value = case_when(
      treatment == "U - E5" ~ "< 0.001",
      TRUE ~ as.character(p_value)
    )
  )

# Carbon-to-nitrogen: Create effect size plot
effects_CN <- ggplot(
  mean_effect_CN,
  aes(x = treatment, y = Coefficient, fill = treatment)) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_errorbar(
    mapping = aes(x = treatment, ymax = CI_high, ymin = CI_low),
    color = "black",
    width = 0,
    linewidth = 0.7
  ) +
  geom_point(
    shape = 16, size = 3,
    colour = alpha(c(
      "#FEAF77FF", "#F1605DFF", "#B63679FF",
      "#721F81FF", "#2D1160FF", "#000004FF"),
      0.7)
  ) +
  scale_y_continuous(
    limits = c(-10, 6.05),
    breaks = c(-10, -5, 0, 5)
    ) +
  xlab("Fire Treatment") +
  ylab(NULL) +
  MyTheme() +
  annotate(
    "text", x = 0.7, y = 6.05 * 0.975, 
    label = "e", size = 5, fontface = "bold"
  ) +
  # Add pairwise comparisons
  stat_pvalue_manual(
    data = pairwise_plot_CN, 
    y.position = 0.70,
    step.increase = 0.032,
    tip.length = 0.005,
    bracket.size = 0.2,
    label.size = 2.5,
    label = "p_value"
  ) +
  coord_cartesian(
    # This allows me to annotate outside panel borders, but it also emboldens
    # the panel borders so include this line in all plots regardless of
    # annotation
    xlim = c(1, 6), clip = "off"
  )

# Ancient grass cover: Format mean effect size data frame
mean_effect_ancient_grass <- parameters(
  model_ancient_grass,
  effects = "fixed"
) %>%
  as_tibble() %>%
  select(
    Coefficient, CI_low, CI_high
  ) %>%
  mutate(
    treatment = c("U", "E5", "E3", "E2", "L2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "E2", "L2", "E1"))
  )

# Ancient grass cover: Format pairwise comarisons data frame
pairwise_plot_ancient_grass <- pairwise_ancient_grass %>%
  select(treatment, p_value) %>%
  filter(p_value <= 0.05) %>%
  mutate(
    group1 = word(treatment, 1, sep = " - "),
    group2 = word(treatment, 2, sep = " - "),
    p_value = case_when(
      treatment == "E3 - E1" ~ "< 0.001",
      TRUE ~ as.character(p_value)
    )
  )

# Ancient grass cover: Create effect size plot
effects_ancient_grass <- ggplot(
  mean_effect_ancient_grass,
  aes(x = treatment, y = Coefficient, fill = treatment)) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_errorbar(
    mapping = aes(x = treatment, ymax = CI_high, ymin = CI_low),
    color = "black",
    width = 0,
    linewidth = 0.7
  ) +
  geom_point(
    shape = 16, size = 3,
    colour = alpha(c(
      "#FEAF77FF", "#F1605DFF", "#B63679FF",
      "#721F81FF", "#2D1160FF", "#000004FF"),
      0.7)
  ) +
  scale_y_continuous(
    limits = c(-4, 4.15),
    breaks = c(-4, -2, 0, 2, 4)
    ) +
  xlab(NULL) +
  ylab("Mean Difference in Effect Size") +
  MyTheme() +
  annotate(
    "text", x = 0.7, y = 4.15 * 0.9775, 
    label = "b", size = 5, fontface = "bold"
  ) +
  # Add pairwise comparisons
  stat_pvalue_manual(
    data = pairwise_plot_ancient_grass, 
    y.position = 2.5,
    step.increase = 0.0482,
    tip.length = 0.005,
    bracket.size = 0.2,
    label.size = 2.5,
    label = "p_value"
  ) +
  coord_cartesian(
    # This allows me to annotate outside panel borders, but it also emboldens
    # the panel borders so include this line in all plots regardless of
    # annotation
    xlim = c(1, 6), clip = "off"
  )

# Modern grass cover: Format mean effect size data frame
mean_effect_modern_grass <- parameters(
  model_modern_grass,
  effects = "fixed"
) %>%
  as_tibble() %>%
  select(
    Coefficient, CI_low, CI_high
  ) %>%
  mutate(
    treatment = c("U", "E5", "E3", "E2", "L2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "E2", "L2", "E1"))
  )

# Modern grass cover: Format pairwise comarisons data frame
pairwise_plot_modern_grass <- pairwise_modern_grass %>%
  select(treatment, p_value) %>%
  filter(p_value <= 0.05) %>%
  mutate(
    group1 = word(treatment, 1, sep = " - "),
    group2 = word(treatment, 2, sep = " - "),
    p_value = case_when(
      treatment == "U - E2" ~ "< 0.001",
      treatment == "U - L2" ~ "< 0.001",
      treatment == "U - E1" ~ "< 0.001",
      treatment == "E5 - L2" ~ "< 0.001",
      TRUE ~ as.character(p_value)
    )
  )

# Modern grass cover: Create effect size plot
effects_modern_grass <- ggplot(
  mean_effect_modern_grass %>%
    mutate(across(where(is.numeric), ~ if_else(row_number() == 1, . * 100, .))) %>%
    ungroup(),
  aes(x = treatment, y = Coefficient, fill = treatment)) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_errorbar(
    mapping = aes(x = treatment, ymax = CI_high, ymin = CI_low),
    color = "black",
    width = 0,
    linewidth = 0.7
  ) +
  geom_point(
    shape = 16, size = 3,
    colour = alpha(c(
      "#FEAF77FF", "#F1605DFF", "#B63679FF",
      "#721F81FF", "#2D1160FF", "#000004FF"),
      0.7)
  ) +
  scale_y_continuous(
    limits = c(-0.5, 12.75),
    breaks = c(0, 2, 4, 6, 8, 10, 12)
  ) +
  xlab(NULL) +
  ylab(NULL) +
  MyTheme() +
  annotate(
    "text", x = 0.7, y = 12.75 * 0.99, 
    label = "c", size = 5, fontface = "bold"
  ) +
  # Add pairwise comparisons
  stat_pvalue_manual(
    data = pairwise_plot_modern_grass, 
    y.position = 5.75,
    step.increase = 0.001825,
    tip.length = 0.0002,
    bracket.size = 0.2,
    label.size = 2.5,
    label = "p_value"
  ) +
  coord_cartesian(
    # This allows me to annotate outside panel borders, but it also emboldens
    # the panel borders so include this line in all plots regardless of
    # annotation
    xlim = c(1, 6), clip = "off"
  )

# Litter cover: Format mean effect size data frame
mean_effect_litter <- parameters(
  model_litter,
  effects = "fixed"
) %>%
  as_tibble() %>%
  select(
    Coefficient, CI_low, CI_high
  ) %>%
  mutate(
    treatment = c("U", "E5", "E3", "E2", "L2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "E2", "L2", "E1"))
  )

# Litter cover: Format pairwise comarisons data frame
pairwise_plot_litter <- pairwise_litter %>%
  select(treatment, p_value) %>%
  filter(p_value <= 0.05) %>%
  mutate(
    group1 = word(treatment, 1, sep = " - "),
    group2 = word(treatment, 2, sep = " - "),
    p_value = case_when(
      treatment == "U - L2" ~ "< 0.001",
      treatment == "E5 - L2" ~ "< 0.001",
      TRUE ~ as.character(p_value)
    )
  )

# Modern grass cover: Create effect size plot
effects_litter <- ggplot(
  mean_effect_litter,
  aes(x = treatment, y = Coefficient, fill = treatment)) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_errorbar(
    mapping = aes(x = treatment, ymax = CI_high, ymin = CI_low),
    color = "black",
    width = 0,
    linewidth = 0.7
  ) +
  geom_point(
    shape = 16, size = 3,
    colour = alpha(c(
      "#FEAF77FF", "#F1605DFF", "#B63679FF",
      "#721F81FF", "#2D1160FF", "#000004FF"),
      0.7)
  ) +
  scale_y_continuous(
    limits = c(-50, 39),
    breaks = c(-40, -20, 0, 20)
  ) +
  xlab("Fire Treatment") +
  ylab("Mean Difference in Effect Size") +
  MyTheme() +
  annotate(
    "text", x = 0.7, y = 39 * 0.975, 
    label = "d", size = 5, fontface = "bold"
  ) +
  # Add pairwise comparisons
  stat_pvalue_manual(
    data = pairwise_plot_litter, 
    y.position = 18,
    step.increase = 0.048,
    tip.length = 0.005,
    bracket.size = 0.2,
    label.size = 2.5,
    label = "p_value"
  ) +
  coord_cartesian(
    # This allows me to annotate outside panel borders, but it also emboldens
    # the panel borders so include this line in all plots regardless of
    # annotation
    xlim = c(1, 6), clip = "off"
  )

# Join and save components b, c, d and e of Figure 2
SEM_plots <-
  effects_ancient_grass + effects_modern_grass + effects_litter +
  effects_CN +
  plot_layout(ncol = 2)
ggsave("output/Figure2_bcde.jpeg", width = 4.4, height = 4.3)
ggsave("output/Figure2_bcde.tiff", width = 4.4, height = 4.3)
ggsave("output/Figure2_bcde.pdf", width = 4.4, height = 4.3)

#### (8g) Additional effects plots: Main text ####

# Total grass cover: Format mean effect size data frame
mean_effect_total_grass <- parameters(
  model_total_grass,
  effects = "fixed"
) %>%
  as_tibble() %>%
  select(
    Coefficient, CI_low, CI_high
  ) %>%
  mutate(
    treatment = c("U", "E5", "E3", "E2", "L2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "E2", "L2", "E1"))
  )

# Total grass cover: Format pairwise comarisons data frame
pairwise_plot_total_grass <- pairwise_total_grass %>%
  select(treatment, p_value) %>%
  filter(p_value <= 0.05) %>%
  mutate(
    group1 = word(treatment, 1, sep = " - "),
    group2 = word(treatment, 2, sep = " - "),
    p_value = case_when(
      treatment == "U - L2" ~ "< 0.001",
      TRUE ~ as.character(p_value)
    )
  )

# Total grass cover: Create effect size plot
effects_total_grass <- ggplot(
  mean_effect_total_grass %>%
    mutate(across(where(is.numeric), ~ if_else(row_number() == 1, . * 100, .))) %>%
    ungroup(),
  aes(x = treatment, y = Coefficient, fill = treatment)) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_errorbar(
    mapping = aes(x = treatment, ymax = CI_high, ymin = CI_low),
    color = "black",
    width = 0,
    linewidth = 0.7
  ) +
  geom_point(
    shape = 16, size = 3,
    colour = alpha(c(
      "#FEAF77FF", "#F1605DFF", "#B63679FF",
      "#721F81FF", "#2D1160FF", "#000004FF"),
      0.7)
  ) +
  scale_y_continuous(
    limits = c(-2, 110),
    breaks = c(0, 20, 40, 60, 80, 100)
  ) +
  xlab("Fire Treatment") +
  ylab("Mean Difference in Effect Size") +
  MyTheme() +
  annotate(
    "text", x = 0.7, y = 110 * 0.985, 
    label = "b", size = 5, fontface = "bold"
  ) +
  # Add pairwise comparisons
  stat_pvalue_manual(
    data = pairwise_plot_total_grass, 
    y.position = 80,
    step.increase = 0.001,
    tip.length = 0.0001,
    bracket.size = 0.2,
    label.size = 2.5,
    label = "p_value"
  ) +
  coord_cartesian(
    # This allows me to annotate outside panel borders, but it also emboldens
    # the panel borders so include this line in all plots regardless of
    # annotation
    xlim = c(1, 6), clip = "off"
  )
effects_total_grass

# Relabel litter plot
effects_litter_2 <- ggplot(
  mean_effect_litter,
  aes(x = treatment, y = Coefficient, fill = treatment)) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_errorbar(
    mapping = aes(x = treatment, ymax = CI_high, ymin = CI_low),
    color = "black",
    width = 0,
    linewidth = 0.7
  ) +
  geom_point(
    shape = 16, size = 3,
    colour = alpha(c(
      "#FEAF77FF", "#F1605DFF", "#B63679FF",
      "#721F81FF", "#2D1160FF", "#000004FF"),
      0.7)
  ) +
  scale_y_continuous(
    limits = c(-50, 39),
    breaks = c(-40, -20, 0, 20)
  ) +
  xlab("Fire Treatment") +
  ylab(NULL) +
  MyTheme() +
  annotate(
    "text", x = 0.7, y = 39 * 0.975, 
    label = "c", size = 5, fontface = "bold"
  ) +
  # Add pairwise comparisons
  stat_pvalue_manual(
    data = pairwise_plot_litter, 
    y.position = 18,
    step.increase = 0.048,
    tip.length = 0.005,
    bracket.size = 0.2,
    label.size = 2.5,
    label = "p_value"
  ) +
  coord_cartesian(
    # This allows me to annotate outside panel borders, but it also emboldens
    # the panel borders so include this line in all plots regardless of
    # annotation
    xlim = c(1, 6), clip = "off"
  )

# Relabel CN plot
effects_CN_2 <- ggplot(
  mean_effect_CN,
  aes(x = treatment, y = Coefficient, fill = treatment)) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_errorbar(
    mapping = aes(x = treatment, ymax = CI_high, ymin = CI_low),
    color = "black",
    width = 0,
    linewidth = 0.7
  ) +
  geom_point(
    shape = 16, size = 3,
    colour = alpha(c(
      "#FEAF77FF", "#F1605DFF", "#B63679FF",
      "#721F81FF", "#2D1160FF", "#000004FF"),
      0.7)
  ) +
  scale_y_continuous(
    limits = c(-10, 6.05),
    breaks = c(-10, -5, 0, 5)
  ) +
  xlab("Fire Treatment") +
  ylab(NULL) +
  MyTheme() +
  annotate(
    "text", x = 0.7, y = 6.05 * 0.975, 
    label = "d", size = 5, fontface = "bold"
  ) +
  # Add pairwise comparisons
  stat_pvalue_manual(
    data = pairwise_plot_CN, 
    y.position = 0.70,
    step.increase = 0.032,
    tip.length = 0.005,
    bracket.size = 0.2,
    label.size = 2.5,
    label = "p_value"
  ) +
  coord_cartesian(
    # This allows me to annotate outside panel borders, but it also emboldens
    # the panel borders so include this line in all plots regardless of
    # annotation
    xlim = c(1, 6), clip = "off"
  )

# Join and save components b and c of Figure 2
SEM_plots <-
  effects_total_grass + effects_litter_2 +
  effects_CN_2 +
  plot_layout(nrow = 1)
ggsave("output/Sup_figure_bcd.jpeg", width = 6.6, height = 2.4)
ggsave("output/Sup_figure_bcd.tiff", width = 6.4, height = 2.4)
ggsave("output/Sup_figure_bcd.pdf", width = 6.4, height = 2.4)

# Clean up the environment
#rm(list = ls())

### (9) Structural equation models #############################################

# Structural equation models are typically complex by design. They are intended
# to contain a given number of paths based on a priori model. Model selection
# typically does follow stepwise processes and traditional selection based
# on information criterion. Log-likelihood χ2 and Fishers C are the most common
# methods to assess SEM fit. However, here we follow a SEM model selection 
# approach outliined by "Garrido et al. (2022). A model selection approach to 
# structural equation modelling: A critical evaluation and a road map for 
# ecologists. Methods in Ecology and Evolution, 13(1), 42-53". This is because
# our sample size is small and including all variable would likely overfit the
# model.

# Read in the data:
data <- read.csv("data/data.csv", stringsAsFactors = TRUE) %>%
  # Remove the nitrogen outlier and NAs
  filter(nitrogen < 0.6) %>%
  mutate(nitrogen = log(nitrogen)) %>%
  mutate(modern_grass_cover = sqrt(modern_grass_cover)) %>%
  mutate(ancient_grass_cover = sqrt(ancient_grass_cover)) %>%
  drop_na() %>%
  # Order levels
  mutate(treatment = factor(treatment, levels = c(
    "U", "E5", "E4", "E3", "L2", "E2", "E1"))
  ) %>%
  # Scale columns 8 to 15
  mutate_at(vars(6:14), ~ as.vector(scale(.))) %>%
  glimpse(.)

# Ensure 'treatment' is a factor
data$treatment <- factor(data$treatment)

# Check the unique levels of 'treatment'
unique_levels <- levels(data$treatment)

# Set sum contrasts
contrasts(data$treatment) <- contr.treatment(unique_levels)

#### (9a) Ancient and modern grass ####

# A priori SEM
sem_ancient_modern_full_model <- psem(
  glmer(
    bacteria_abundance ~ treatment + modern_grass_cover + ancient_grass_cover +
      litter_cover + organic_carbon + nitrogen + carbon_nitrogen_ratio + pH + 
    (1 | block/sample_id),
  family = 'poisson',
  data
  ),
  
  glmer(
    fungi_abundance ~ treatment + modern_grass_cover + ancient_grass_cover +
      litter_cover + organic_carbon + nitrogen + carbon_nitrogen_ratio + pH + 
      (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    fungi_bacteria_ratio ~ treatment + modern_grass_cover + ancient_grass_cover +
      litter_cover + organic_carbon + nitrogen + carbon_nitrogen_ratio + pH + 
      (1 | block),
    data
  ),
  
  lmer(
    modern_grass_cover ~ treatment + litter_cover + organic_carbon + nitrogen +
      carbon_nitrogen_ratio + pH + (1 | block),
    data
  ),
  
  lmer(
    ancient_grass_cover ~ treatment + litter_cover + organic_carbon + nitrogen +
      carbon_nitrogen_ratio + pH + (1 | block),
    data
  ),
  
  lmer(
    litter_cover ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    organic_carbon ~ treatment + litter_cover + modern_grass_cover + 
      ancient_grass_cover + (1 | block),
    data
  ),
  
  lmer(
    nitrogen ~ treatment + litter_cover + modern_grass_cover + 
      ancient_grass_cover +  (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + litter_cover + modern_grass_cover + 
      ancient_grass_cover + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + litter_cover + modern_grass_cover + 
      ancient_grass_cover + (1 | block),
    data
  )
)

# Reduced model after backward stepwise selection
sem_ancient_modern_reduced_model <- psem(
  
  glmer(
    bacteria_abundance ~ treatment + modern_grass_cover + ancient_grass_cover +
      litter_cover + organic_carbon + nitrogen +
      (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  glmer(
    fungi_abundance ~ modern_grass_cover + ancient_grass_cover +
      litter_cover + carbon_nitrogen_ratio + nitrogen +
      (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    fungi_bacteria_ratio ~
      organic_carbon + nitrogen + carbon_nitrogen_ratio +
      (1 | block),
    data
  ),
  
  lmer(
    modern_grass_cover ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    ancient_grass_cover ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    litter_cover ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    pH ~
      litter_cover + (1 | block),
    data
  )
)

#### (9b) Total grass ####

# A priori SEM
sem_total_full_model <- psem(
  glmer(
    bacteria_abundance ~ treatment + total_grass_cover +
      litter_cover + organic_carbon + nitrogen + carbon_nitrogen_ratio + pH + 
      (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  glmer(
    fungi_abundance ~ treatment + total_grass_cover +
      litter_cover + organic_carbon + nitrogen + carbon_nitrogen_ratio + pH + 
      (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    fungi_bacteria_ratio ~ treatment + total_grass_cover +
      litter_cover + organic_carbon + nitrogen + carbon_nitrogen_ratio + pH + 
      (1 | block),
    data
  ),
  
  lmer(
    total_grass_cover ~ treatment + litter_cover + organic_carbon + nitrogen +
      carbon_nitrogen_ratio + pH + (1 | block),
    data
  ),

  lmer(
    litter_cover ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    organic_carbon ~ treatment + litter_cover + total_grass_cover + 
      (1 | block),
    data
  ),
  
  lmer(
    nitrogen ~ treatment + litter_cover + total_grass_cover + 
      (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + litter_cover + total_grass_cover + 
      (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + litter_cover + total_grass_cover + 
      (1 | block),
    data
  )
)

# Reduced model after backward stepwise selection
sem_total_reduced_model <- psem(
  glmer(
    bacteria_abundance ~ treatment + total_grass_cover +
      litter_cover + organic_carbon + nitrogen +
      (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  glmer(
    fungi_abundance ~ treatment +
      litter_cover + nitrogen + carbon_nitrogen_ratio +
      (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    fungi_bacteria_ratio ~
      carbon_nitrogen_ratio +
      (1 | block),
    data
  ),
  
  lmer(
    total_grass_cover ~ treatment +
      (1 | block),
    data
  ),
  
  lmer(
    litter_cover ~ treatment +
      (1 | block),
    data
  ),
  
  lmer(
    organic_carbon ~
      litter_cover + (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + total_grass_cover + 
      (1 | block),
    data
  ),
  
  lmer(
    pH ~ total_grass_cover +
      (1 | block),
    data
  )
)
