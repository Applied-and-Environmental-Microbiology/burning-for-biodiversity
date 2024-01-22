
# Here I build I will assess the impact of fire treatment on soil microbial 
# abundance and chemical properties using structural equation models. I will
# first make an a prior model based on these assumptions:
#   soil properties ~ fire treatment
#   microbial abundance ~ fire treatment
#   microbial abundance ~ soil properties
#
# Model fit will be assessed on AIC (balance of goodness of fit and complexity)
# and Log-Likelihood Chi-Squared (χ²) (measure of goodness of fit against the
# null model).

# Required packages and my functions
require(piecewiseSEM)
require(lme4)
require(DHARMa)
require(performance)
require(tidyverse)

# Read in the data:
data <- read.csv("data/data.csv", stringsAsFactors = TRUE) %>%
  # Remove the nitrogen outlier
  filter(nitrogen < 0.6) %>%
  # Order levels
  mutate(treatment = factor(treatment, levels = c(
    "U", "E5", "E4", "E3", "L2", "E2", "E1"))
  ) %>%
  glimpse(.)


sem_fungi <- psem(
  
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
  ))

summary(sem_fungi)

#### (9b) Bacteria ####
model_bacteria <- glmer(
  bacteria_abundance ~ treatment + organic_carbon + log(nitrogen) + 
    carbon_nitrogen_ratio + pH + (1 | block/sample_id),
  family = 'poisson',
  data
)
summary(model_bacteria)

model_bacteria_1 <- glmer(
  bacteria_abundance ~ treatment + organic_carbon + log(nitrogen) + 
    carbon_nitrogen_ratio + (1 | block/sample_id),
  family = 'poisson',
  data
)
summary(model_bacteria_1)

model_bacteria_2 <- glmer(
  bacteria_abundance ~ treatment + organic_carbon +
    carbon_nitrogen_ratio + (1 | block/sample_id),
  family = 'poisson',
  data
)
summary(model_bacteria_2)

model_bacteria_3 <- glmer(
  bacteria_abundance ~ treatment + organic_carbon + 
    (1 | block/sample_id),
  family = 'poisson',
  data
)

compare_performance(
  model_bacteria, model_bacteria_1, model_bacteria_2,model_bacteria_3
)

check_collinearity(model_bacteria_2)

DHARMa::simulateResiduals(model_bacteria_2) %>%
  plot(.)

sem_bacteria <- psem(
  
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
  ))

summary(sem_bacteria)





sem_fungi <- psem(
  
  model_fungi_1 <- glmer(
    fungi_abundance ~ treatment + organic_carbon + log(nitrogen) + 
      carbon_nitrogen_ratio + pH 
    + (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    organic_carbon ~ treatment + log(nitrogen) + carbon_nitrogen_ratio + pH + 
      (1 | block),
    data
  ),
  
  lmer(
    log(nitrogen) ~ treatment + carbon_nitrogen_ratio + pH + (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + pH + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  ))

summary(sem_fungi)
AIC_psem(sem_fungi)




sem_fungi_full <- psem(
  
  glmer(
    fungi_abundance ~ treatment + organic_carbon + log(nitrogen) + 
      carbon_nitrogen_ratio + pH + (1 | block/sample_id),
    family = 'poisson',
    data
  ),
  
  lmer(
    organic_carbon ~ treatment + log(nitrogen) + carbon_nitrogen_ratio + pH + 
      (1 | block),
    data
  ),
  
  lmer(
    log(nitrogen) ~ treatment + carbon_nitrogen_ratio + pH + (1 | block),
    data
  ),
  
  lmer(
    carbon_nitrogen_ratio ~ treatment + pH + (1 | block),
    data
  ),
  
  lmer(
    pH ~ treatment + (1 | block),
    data
  ))


summary(sem_fungi_full)
AIC_psem(sem_fungi_full)

# Given AIC values
aic_values <- c(574.1, 571.6, 569.1, 394.2, 747.8)

# Step 2: Identify the model with the lowest AIC
best_model_index <- which.min(aic_values)

# Step 3: Calculate the difference in AIC for each model
delta_aic <- aic_values - aic_values[best_model_index]

# Step 4: Compute unnormalised weights
weights <- exp(-delta_aic / 2)

# Step 5: Normalise weights
normalised_weights <- weights / sum(weights)

# Display the results
results <- data.frame(Model = 1:length(aic_values), AIC = aic_values,
                      Delta_AIC = delta_aic, Unnormalised_weight = weights,
                      Normalised_weight = normalised_weights)
print(results)
