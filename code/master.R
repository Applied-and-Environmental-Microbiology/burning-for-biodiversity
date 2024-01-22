
# Script: Three Parks Savanna Fire-effects Plot Network: The effect of fire
#         treatment on soil properties
# Author: Luke Florence
# Date:   21 Janurary 2024

# Note: Models have previously undergone fit tests and have been assessed for 
# appropriate structure and transformations. Models fit on a Gaussian
# distribution are best fit with "block" as a random effect and Poisson
# distribution with "block/sample" as random effects. Nevertheless, because I
# have undertaken a novel approach to SEM model selection (see notes in section
# (9) Structural equation models) I include the model selection process within
# the corresponding section.

# Required packages and functions:
require(lme4)
require(emmeans)
require(performance)
require(parameters)
require(DHARMa)
require(openxlsx)
require(patchwork)
require(ggpubr)
require(tidyverse)
source("code/functions.R")

# Read in the data:
data <- read.csv("data/data.csv", stringsAsFactors = TRUE) %>%
  # Order levels
  mutate(treatment = factor(treatment, levels = c(
    "U", "E5", "E4", "E3", "L2", "E2", "E1"))
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

### (8) Save test statistics ###################################################

# List of response variables
response_variables <- c(
  "fungi", "bacteria", "carbon","ratio", "nitrogen", "CN", "pH"
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
    "text", x = -0.75, y = (max(observed_values_fungi$rate) + 0.01), 
    label = "b", size = 5, fontface = "bold"
    ) +
  coord_cartesian(
    # This allows me to annotate outside panel borders, but it also emboldens
    # the panel borders so Iinclude this line in all plots regardless of
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
    treatment = c("U", "E5", "E3", "L2", "E2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "L2", "E2", "E1"))
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
  xlab("Fire treatment") +
  ylab(NULL) +
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
    "text", x = -0.75, y = 15, 
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
    treatment = c("U", "E5", "E3", "L2", "E2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "L2", "E2", "E1"))
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
  xlab(NULL) +
  ylab("Organic Carbon (%)") +
  MyTheme() +
  theme(
    axis.text.x = element_blank()
    ) + 
  annotate(
    "text", x = -0.75, y = 3, 
    label = "c", size = 5, fontface = "bold"
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
    treatment = c("U", "E5", "E3", "L2", "E2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "L2", "E2", "E1"))
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
  xlab(NULL) +
  ylab(NULL) +
  MyTheme() +
  coord_cartesian(
    xlim = c(1, 6), clip = "off"
    )

#### (8d) Join and save Figure 2 ####
main_plots <- 
  means_bacteria + means_fungi + means_carbon + 
  effects_bacteria + effects_fungi + effects_carbon + 
  plot_layout(nrow = 2)
ggsave("output/Figure1.jpeg", width = 7, height = 4.5)
ggsave("output/Figure1.tiff", width = 7, height = 4.5)
ggsave("output/Figure1.pdf", width = 7, height = 4.5)

#### (8e) Additional effects plots ####

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
    treatment = c("U", "E5", "E3", "L2", "E2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "L2", "E2", "E1"))
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
    limits = c(-10, 6),
    breaks = c(-10, -5, 0, 5)
    ) +
  xlab(NULL) +
  ylab("Mean Difference in Effect Size") +
  MyTheme() +
  annotate(
    "text", x = 0.7, y = 5.75, 
    label = "b", size = 5, fontface = "bold"
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
  )

# pH: Format mean effect size data frame
mean_effect_pH <- parameters(
  model_pH,
  effects = "fixed"
) %>%
  as_tibble() %>%
  select(
    Coefficient, CI_low, CI_high
  ) %>%
  mutate(
    treatment = c("U", "E5", "E3", "L2", "E2", "E1"),
    treatment = factor(treatment, levels = c(
      "U", "E5", "E3", "L2", "E2", "E1"))
  )
# pH: Format pairwise comarisons data frame
pairwise_plot_pH <- pairwise_pH %>%
  select(treatment, p_value) %>%
  filter(p_value <= 0.05) %>%
  mutate(
    group1 = word(treatment, 1, sep = " - "),
    group2 = word(treatment, 2, sep = " - "),
    p_value = case_when(
      treatment == "U - E3" ~ "0.050",
      TRUE ~ as.character(p_value)
  )
)

# pH: Create effect size plot
effects_pH <- ggplot(
  mean_effect_pH,
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
    limits = c(-0.2, 0.515),
    breaks = c(-0.2, 0, 0.2, 0.4)
    ) +
  xlab(NULL) +
  ylab("Mean Difference in Effect Size") +
  MyTheme() +
  annotate(
    "text", x = 0.7, y = 0.515, 
    label = "c", size = 5, fontface = "bold"
  ) +
  # Add pairwise comparisons
  stat_pvalue_manual(
    data = pairwise_plot_pH, 
    y.position = 0.40,
    step.increase = 0.0125,
    tip.length = 0.005,
    bracket.size = 0.2,
    label.size = 2.5,
    label = "p_value"
  )

# Join and save components b and c of Figure 2
SEM_plots <- 
  effects_CN + effects_pH + 
  plot_layout(ncol = 1)
ggsave("output/Figure2.jpeg", width = 2.5, height = 4.5)
ggsave("output/Figure2.tiff", width = 2.5, height = 4.5)
ggsave("output/Figure2.pdf", width = 2.5, height = 4.5)

### (9) Structural equation models #############################################

# Structural equation models are typically complex by design. They are intended
# to contain X number of paths based on a priori model. As such, model selection
# typically does follow stepwise processes and traditional selection based
# on an information criterion. In contrast, log-likelihood χ2 and Fishers C are
# common method to assess SEM fit, and may result in the addition rather than 
# subtraction of paths to capture underlying/missing interactions - but a priori
# knowledge is a key component of SEMs and may allow a researcher to disregard
# such statistics to a certain extent. However, a recent paper has identified
# the need for a more 'traditional' SEM model selection approach in cases
# where, for example, the variance of dependent variables is high, there are
# high levels of multicollinearity, and sample size is relatively small 
# (Garrido et al. (2022). A model selection approach to structural equation
# modelling: A critical evaluation and a road map for ecologists. Methods in
# Ecology and Evolution, 13(1), 42-53), all of which are present in our 
# dataset. To that end, here I will fit fully saturated SEMs and follow a
# backward stepwise approach the removal of the most insignificant variable
# and stopping criteria of ΔAICc ≤ 2. I will then select my model based on AICc
# and normalised model weight a using a ΔAICc ≤ 2 and normalised model weight 
# change ≥ 55% as selection thresholds.

# Read in the data:
data <- read.csv("data/data.csv", stringsAsFactors = TRUE) %>%
  # Remove the nitrogen outlier
  filter(nitrogen < 0.6) %>%
  # Order levels
  mutate(treatment = factor(treatment, levels = c(
    "U", "E5", "E4", "E3", "L2", "E2", "E1"))
  ) %>%
  glimpse(.)

#### (9a) Fungi ####

# A priori SEM
sem_fungi_full <- psem(
  glmer(
  fungi_abundance ~ treatment + organic_carbon + log(nitrogen) + 
    carbon_nitrogen_ratio + pH + (1 | block/sample_id),
  family = 'poisson',
  data
  ),
  
  lmer(
    organic_carbon ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    log(nitrogen) ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  )
  )
# Evaluate Significance values
coefs(sem_fungi_full) %>%
  select(Response, Predictor, Crit.Value, P.Value)
# Remove 'organic_carbon' from the 'fungi_abundance' model
# Extract AICc
aic_fungi_full <- AIC_psem(sem_fungi_full)$AICc

# Reduced SEM 1
sem_fungi_1 <- psem(
  glmer(
    fungi_abundance ~ treatment + log(nitrogen) + 
      carbon_nitrogen_ratio + pH + (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    organic_carbon ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    log(nitrogen) ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  )
)
# Evaluate Significance values
coefs(sem_fungi_1) %>%
  select(Response, Predictor, Crit.Value, P.Value)
# Remove 'pH' from the 'fungi_abundance' model
# Extract AICc
aic_fungi_1 <- AIC_psem(sem_fungi_1)$AICc
# Compute delta AIC
aic_fungi_full - aic_fungi_1

# Reduced model 2
sem_fungi_2 <- psem(
  glmer(
    fungi_abundance ~ treatment + log(nitrogen) + 
      carbon_nitrogen_ratio + (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    organic_carbon ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    log(nitrogen) ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  )
)
# Evaluate Significance values
coefs(sem_fungi_2) %>%
  select(Response, Predictor, Crit.Value, P.Value)
# Remove the 'organic_carbon ~ treatmente' model
# Extract AICc
aic_fungi_2 <- AIC_psem(sem_fungi_2)$AICc
# Compute delta AIC
aic_fungi_1 - aic_fungi_2

# Reduced model 3
sem_fungi_3 <- psem(
  glmer(
    fungi_abundance ~ treatment + log(nitrogen) + 
      carbon_nitrogen_ratio + (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    log(nitrogen) ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  )
)
# Evaluate Significance values
coefs(sem_fungi_3) %>%
  select(Response, Predictor, Crit.Value, P.Value)
# Remove the 'log(nitrogen) ~ treatment' model
# Extract AICc
aic_fungi_3 <- AIC_psem(sem_fungi_3)$AICc
# Compute delta AIC
aic_fungi_2 - aic_fungi_3

# Reduced model 4
sem_fungi_4 <- psem(
  glmer(
    fungi_abundance ~ treatment + log(nitrogen) + 
      carbon_nitrogen_ratio + (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  )
)
# Evaluate Significance values
coefs(sem_fungi_4) %>%
  select(Response, Predictor, Crit.Value, P.Value)
# Remove 'treatment' from the 'fungi_abundance' model
# Extract AICc
aic_fungi_4 <- AIC_psem(sem_fungi_4)$AICc
# Compute delta AIC
aic_fungi_3 - aic_fungi_4

# Reduced model 5
sem_fungi_5 <- psem(
  glmer(
    fungi_abundance ~ log(nitrogen) + 
      carbon_nitrogen_ratio + (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  )
)
# Evaluate Significance values
coefs(sem_fungi_5) %>%
  select(Response, Predictor, Crit.Value, P.Value)
# Remove 'carbon_nitrogen_ratio' from 'fungi_abundance' model
# Extract AICc
aic_fungi_5 <- AIC_psem(sem_fungi_5)$AICc
# Compute delta AIC
aic_fungi_4 - aic_fungi_5

# Reduced model 6
sem_fungi_6 <- psem(
  glmer(
    fungi_abundance ~ log(nitrogen) + 
      (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  )
)
# Evaluate Significance values
coefs(sem_fungi_6) %>%
  select(Response, Predictor, Crit.Value, P.Value)
# Remove the 'pH ~ treatment' model
# Extract AICc
aic_fungi_6 <- AIC_psem(sem_fungi_6)$AICc
# Compute delta AIC
aic_fungi_5 - aic_fungi_6
# !!!Stop!!!

# Create a list of models
sem_models_fungi <- list(
  sem_fungi_full, sem_fungi_1, sem_fungi_2, sem_fungi_3, sem_fungi_4,
  sem_fungi_5
)

# Combine LLChisq stats:

# Step 1: Extract and format LLChisq stats
LLChisq_fungi <- map(
  sem_models_fungi, ~ LLchisq(.x) %>%
    mutate(
      formatted_result = sprintf("%.2f(%d), p = %.3f", Chisq, df, P.Value)
      ) %>%
    pull(formatted_result)
)

# Step 2: Combine the formatted results into a single vector
LLChisq_fungi <- unlist(LLChisq_fungi)

# Combine FisherC stats:

# Step 1: Create an empty list to store results
Fisher_C_fungi_list <- list()

# Step 2: Loop through the list of SEM models
for (model in sem_models_fungi) {
  # Extract FisherC stats
  fisher_c_stats <- summary(model) %>%
    unlist(.) %>%
    as.data.frame(
      ., stringsAsFactors = FALSE
    ) %>%
    rownames_to_column(var = "name") %>%
    filter(
      name %in% c("Cstat.Fisher.C", "Cstat.df", "Cstat.P.Value")
    ) %>%
    column_to_rownames(var = "name") %>%
    mutate(across(everything(), as.numeric)) %>%
    t(.) %>%
    as_tibble(.) %>%
    mutate(
      formatted_result = sprintf(
        "%.2f(%d), p = %.3f",
        Cstat.Fisher.C, Cstat.df, Cstat.P.Value)
    ) %>%
    pull(formatted_result)
  
  # Append the results to the list
  Fisher_C_fungi_list[[length(Fisher_C_fungi_list) + 1]] <- fisher_c_stats
}

# Step 3: Combine the formatted results into a single vector
Fisher_C_fungi <- unlist(Fisher_C_fungi_list)

# Combine model structures:

# Step 1: Create an empty list to store results
model_call_list <- list()

# Step 2: Loop through the list of SEM models
for (model in sem_models_fungi) {
  # Extract call information
  call_info <- summary(model) %>%
    unlist(.) %>%
    as.data.frame(
      ., stringsAsFactors = FALSE
    ) %>%
    slice(2) %>%
    pull()
  
  # Assign column names to the call information
  names(call_info) <- "call"
  
  # Store the call information in the list
  model_call_list[[length(model_call_list) + 1]] <- call_info
}

# Step 3: Reduce to a single data frame
model_call_df <- bind_rows(model_call_list, .id = "model")

# Compute standardised model weights:

# Step 1: Combine AIC values
aic_values <- c(
  aic_fungi_full, aic_fungi_1, aic_fungi_2, aic_fungi_3, aic_fungi_4,
  aic_fungi_5
  )

# Step 2: Identify the model with the lowest AIC
best_model_index <- which.min(aic_values)

# Step 3: Calculate the difference in AIC for each model
delta_aic <- aic_values - aic_values[best_model_index]

# Step 4: Compute unnormalised weights
weights <- exp(-delta_aic / 2)

# Step 5: Normalise weights
normalised_weights <- weights / sum(weights)

# Step6: Assess the results
SEM_fit_fungi <- data.frame(
  model = c("full", 2:length(aic_values)-1), 
  AIC = aic_values, 
  Δ_AIC = delta_aic,
  model_weight = round(normalised_weights, digits = 3)
) %>%
  mutate(
    Δ_model_weight = c(0, diff(model_weight)),
    "LL_χ2(df)" = LLChisq_fungi,
    "Fisher_C(df)" = Fisher_C_fungi,
    model_structure = model_call_df$call
  ) %>%
  print(.)

# Grab coefficients of the best-fit model
SEM_coefs_fungi <- coefs(sem_fungi_5)

#### (9b) Bacteria ##### 

# A priori SEM
sem_bacteria_full <- psem(
  glmer(
    bacteria_abundance ~ treatment + organic_carbon + log(nitrogen) + 
      carbon_nitrogen_ratio + pH + (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    organic_carbon ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    log(nitrogen) ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  )
)
# Evaluate Significance values
coefs(sem_bacteria_full) %>%
  select(Response, Predictor, Crit.Value, P.Value)
# Remove the 'organic_carbon ~ treatment) model
# Extract AICc
aic_bacteria_full <- AIC_psem(sem_bacteria_full)$AICc

# Reduced SEM 1
sem_bacteria_1 <- psem(
  glmer(
    bacteria_abundance ~ treatment + organic_carbon + log(nitrogen) +
      carbon_nitrogen_ratio + pH + (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    log(nitrogen) ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  )
)
# Evaluate Significance values
coefs(sem_bacteria_1) %>%
  select(Response, Predictor, Crit.Value, P.Value)
# Remove 'log(nitrogen)' from the 'bacteria_abundance' model
# Extract AICc
aic_bacteria_1 <- AIC_psem(sem_bacteria_1)$AICc
# Compute delta AIC
aic_bacteria_full - aic_bacteria_1

# Reduced model 2
sem_bacteria_2 <- psem(
  glmer(
    bacteria_abundance ~ treatment + organic_carbon +
      carbon_nitrogen_ratio + pH + (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    log(nitrogen) ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  )
)
# Evaluate Significance values
coefs(sem_bacteria_2) %>%
  select(Response, Predictor, Crit.Value, P.Value)
# Remove 'pH' from the 'bacteria_abundance' model
# Extract AICc
aic_bacteria_2 <- AIC_psem(sem_bacteria_2)$AICc
# Compute delta AIC
aic_bacteria_1 - aic_bacteria_2

# Reduced model 3
sem_bacteria_3 <- psem(
  glmer(
    bacteria_abundance ~ treatment + organic_carbon +
      carbon_nitrogen_ratio + (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    log(nitrogen) ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  )
)
# Evaluate Significance values
coefs(sem_bacteria_3) %>%
  select(Response, Predictor, Crit.Value, P.Value)
# Remove the 'log(nitrogen) ~ treatment' model
# Extract AICc
aic_bacteria_3 <- AIC_psem(sem_bacteria_3)$AICc
# Compute delta AIC
aic_bacteria_2 - aic_bacteria_3

# Reduced model 4
sem_bacteria_4 <- psem(
  glmer(
    bacteria_abundance ~ treatment + organic_carbon +
      carbon_nitrogen_ratio + (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  )
)
# Evaluate Significance values
coefs(sem_bacteria_4) %>%
  select(Response, Predictor, Crit.Value, P.Value)
# Remove 'pH ~ treatment' model
# Extract AICc
aic_bacteria_4 <- AIC_psem(sem_bacteria_4)$AICc
# Compute delta AIC
aic_bacteria_3 - aic_bacteria_4

# Reduced model 5
sem_bacteria_5 <- psem(
  glmer(
    bacteria_abundance ~ treatment + organic_carbon +
      carbon_nitrogen_ratio + (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + (1 | block),
    data
  )
)
# Evaluate Significance values
coefs(sem_bacteria_5) %>%
  select(Response, Predictor, Crit.Value, P.Value)
# Remove 'bacteria_abundance' from 'bacteria_abundance' model
# Extract AICc
aic_bacteria_5 <- AIC_psem(sem_bacteria_5)$AICc
# Compute delta AIC
aic_bacteria_4 - aic_bacteria_5
# !!!Stop!!!

# Create a list of models
sem_models_bacteria <- list(
  sem_bacteria_full, sem_bacteria_1, sem_bacteria_2, sem_bacteria_3, 
  sem_bacteria_4
)

# Combine LLChisq stats

# Step 1: Extract and format LLChisq stats
LLChisq_bacteria <- map(
  sem_models_bacteria, ~ LLchisq(.x) %>%
    mutate(
      formatted_result = sprintf("%.2f(%d), p = %.3f", Chisq, df, P.Value)
    ) %>%
    pull(formatted_result)
)

# Step 2: Combine the formatted results into a single vector
LLChisq_bacteria <- unlist(LLChisq_bacteria)

# Combine FisherC stats:

# Step 1: Create an empty list to store results
Fisher_C_bacteria_list <- list()

# Step 2: Loop through the list of SEM models
for (model in sem_models_bacteria) {
  # Extract FisherC stats
  fisher_c_stats <- summary(model) %>%
    unlist(.) %>%
    as.data.frame(
      ., stringsAsFactors = FALSE
    ) %>%
    rownames_to_column(var = "name") %>%
    filter(
      name %in% c("Cstat.Fisher.C", "Cstat.df", "Cstat.P.Value")
    ) %>%
    column_to_rownames(var = "name") %>%
    mutate(across(everything(), as.numeric)) %>%
    t(.) %>%
    as_tibble(.) %>%
    mutate(
      formatted_result = sprintf(
        "%.2f(%d), p = %.3f",
        Cstat.Fisher.C, Cstat.df, Cstat.P.Value)
    ) %>%
    pull(formatted_result)
  
  # Append the results to the list
  Fisher_C_bacteria_list[[length(Fisher_C_bacteria_list) + 1]] <- fisher_c_stats
}

# Step 3: Combine the formatted results into a single vector
Fisher_C_bacteria <- unlist(Fisher_C_bacteria_list)

# Combine model structures:

# Step 1: Create an empty list to store results
model_call_list <- list()

# Step 2: Loop through the list of SEM models
for (model in sem_models_bacteria) {
  # Extract call information
  call_info <- summary(model) %>%
    unlist(.) %>%
    as.data.frame(
      ., stringsAsFactors = FALSE
    ) %>%
    slice(2) %>%
    pull()
  
  # Assign column names to the call information
  names(call_info) <- "call"
  
  # Store the call information in the list
  model_call_list[[length(model_call_list) + 1]] <- call_info
}

# Step 3: Reduce to a single data frame
model_call_df <- bind_rows(model_call_list, .id = "model")

# Compute standardised model weights:

# Step 1: Combine AIC values
aic_values <- c(
  aic_bacteria_full, aic_bacteria_1, aic_bacteria_2, aic_bacteria_3,
  aic_bacteria_4
)

# Step 2: Identify the model with the lowest AIC
best_model_index <- which.min(aic_values)

# Step 3: Calculate the difference in AIC for each model
delta_aic <- aic_values - aic_values[best_model_index]

# Step 4: Compute unnormalised weights
weights <- exp(-delta_aic / 2)

# Step 5: Normalise weights
normalised_weights <- weights / sum(weights)

# Step6: Assess the results
SEM_fit_bacteria <- data.frame(
  model = c("full", 2:length(aic_values)-1), 
  AIC = aic_values, 
  Δ_AIC = delta_aic,
  model_weight = round(normalised_weights, digits = 3)
) %>%
  mutate(
    Δ_model_weight = c(0, diff(model_weight)),
    "LL_χ2(df)" = LLChisq_bacteria,
    "Fisher_C(df)" = Fisher_C_bacteria,
    model_structure = model_call_df$call
  ) %>%
  print(.)

# Grab coefficients of the best-fit model
SEM_coefs_bacteria <- coefs(sem_bacteria_4)

#### (9c) Save the SEM test statistics ####

# Read in the test statistics workbook
test_statistics_wb <- loadWorkbook(
  "output/test_statistics.xlsx"
  )

# Define sheet names
sheet_names <- c(
  "SEM_fit_fungi", "SEM_fit_bacteria", "SEM_coefs_fungi", "SEM_coefs_bacteria"
  )

# Loop through sheet names
for (sheet_name in sheet_names) {
  # Add worksheet
  addWorksheet(test_statistics_wb, sheet_name)
  
  # Get the corresponding data frame based on sheet_name
  df <- get(sheet_name)
  
  # Write data to the sheet
  writeData(test_statistics_wb, sheet = sheet_name, x = df)
}

# Check the sheets have been writen
read.xlsx(test_statistics_wb, sheet = "SEM_fit_fungi")

# Save the workbook
saveWorkbook(test_statistics_wb, "output/test_statistics.xlsx")
