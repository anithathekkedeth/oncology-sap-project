# Create dataset

n_A <- 120
n_B <- 118

# Treatment A (45% responders → 54)
A <- data.frame(
  id = 1:n_A,
  treatment = "A",
  response = c(rep(1, 54), rep(0, n_A - 54)),
  pdl1 = rep(c("<1%", "1-49%", ">=50%"), length.out = n_A),
  time = rexp(n_A, rate = 0.08),
  event = rbinom(n_A, 1, 0.7)
)

# Treatment B (27.1% responders → 32)
B <- data.frame(
  id = (n_A + 1):(n_A + n_B),
  treatment = "B",
  response = c(rep(1, 32), rep(0, n_B - 32)),
  pdl1 = rep(c("<1%", "1-49%", ">=50%"), length.out = n_B),
  time = rexp(n_B, rate = 0.12),
  event = rbinom(n_B, 1, 0.7)
)

data <- rbind(A, B)
head(data)

#ORR TABLE (Table 14.1.1)
library(dplyr)
library(broom)

data$treatment <- relevel(as.factor(data$treatment), ref = "B")

# Summary
summary_tbl <- data %>%
  group_by(treatment) %>%
  summarise(
    N = n(),
    Responders = sum(response),
    ORR = mean(response)
  )

# Logistic model
fit <- glm(response ~ treatment, data = data, family = binomial)

model_res <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)

OR <- model_res %>% filter(term == "treatmentA")

# Final formatted table
final_table <- summary_tbl %>%
  mutate(
    `Responders (n, %)` = paste0(Responders, " (", round(ORR*100,1), "%)")
  ) %>%
  select(treatment, N, `Responders (n, %)`) %>%
  rename(`Treatment Arm` = treatment)

final_table$`Odds Ratio` <- c(round(OR$estimate, 2), NA)
final_table$`95% CI` <- c(
  paste0("(", round(OR$conf.low,2), ", ", round(OR$conf.high,2), ")"),
  NA
)

print(final_table)

library(knitr)

kable(final_table, caption = "Table 14.1.1: ORR Summary")
# SUBGROUP
library(dplyr)
library(broom)

subgroup <- data %>%
  group_by(pdl1) %>%
  do({
    fit <- glm(response ~ treatment, data = ., family = binomial)
    res <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
    
    res %>%
      filter(term == "treatmentA") %>%
      select(estimate, conf.low, conf.high)
  })

names(subgroup) <- c("Subgroup", "OR", "Lower", "Upper")

subgroup$`95% CI` <- paste0("(", round(subgroup$Lower,2), ", ", round(subgroup$Upper,2), ")")

subgroup <- subgroup %>%
  select(Subgroup, OR, `95% CI`)

subgroup
# KM
library(survival)
library(survminer)


# Factor 
data$treatment <- factor(data$treatment, levels = c("B", "A"))

fit <- survfit(Surv(time, event) ~ treatment, data = data)

ggsurvplot(
  fit,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  legend.title = "Treatment",
  legend.labs = c("SOC", "Immunotherapy A"),
  title = "Figure 14.2.1: Kaplan–Meier Curve (PFS)"
)


#Forest


library(dplyr)
library(broom)
library(ggplot2)


data$treatment <- relevel(as.factor(data$treatment), ref = "B")

subgroup <- data %>%
  group_by(pdl1) %>%
  do({
    fit <- glm(response ~ treatment, data = ., family = binomial)
    tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term == "treatmentA")
  })

# Fix ordering
subgroup$pdl1 <- factor(subgroup$pdl1,
                        levels = c("<1%", "1-49%", ">=50%"))

# Add labels
subgroup$label <- paste0(
  round(subgroup$estimate,2),
  " (",
  round(subgroup$conf.low,2), ", ",
  round(subgroup$conf.high,2), ")"
)

ggplot(subgroup, aes(x = estimate, y = pdl1)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  
  # Add labels
  geom_text(aes(label = label), hjust = -0.2, size = 4) +
  
  # Reference line
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # Better scaling
  coord_cartesian(xlim = c(0, 4)) +
  
  # Labels
  xlab("Odds Ratio (A vs SOC)") +
  ylab("PD-L1 Group") +
  ggtitle("Figure 14.2.2: Subgroup Forest Plot") +
  
  theme_minimal()
# Baysian
library(rstanarm)
library(bayesplot)
library(ggplot2)


data$treatment <- relevel(as.factor(data$treatment), ref = "B")

model <- stan_glm(
  response ~ treatment + pdl1,
  data = data,
  family = binomial,
  prior = normal(0, 2.5),
  prior_intercept = normal(0, 5),
  chains = 4,
  iter = 2000,
  seed = 123
)

posterior <- as.matrix(model)

# Convert to OR scale
posterior_OR <- exp(posterior[, "treatmentA"])

df <- data.frame(OR = posterior_OR)

ggplot(df, aes(x = OR)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  
  # Reference line (no effect)
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # Clinical threshold
  geom_vline(xintercept = 1.2, color = "red") +
  
  # Labels
  labs(
    title = "Figure 14.2.3: Bayesian Posterior Distribution",
    x = "Odds Ratio (A vs SOC)",
    y = "Density"
  ) +
  
  theme_minimal()