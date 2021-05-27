library(ggplot2)
library(data.table)


# Discrete covariate, independent -----------------------------------------

cov_indep = readRDS("data/2021-05-26_covariate_independent.RDS")

cov_indep_long = melt(cov_indep, id.vars = c("individual", "simulation"))

ggplot(cov_indep_long, aes(x = value)) +
  geom_boxplot() +
  facet_wrap(vars(variable)) +
  coord_flip()


# No covariate ------------------------------------------------------------

no_cov = readRDS("data/2021-05-26_no_covariate.RDS")

no_cov_long = melt(no_cov, id.vars = c("individual", "simulation"))

ggplot(no_cov_long, aes(x = value)) +
  geom_boxplot() +
  facet_wrap(vars(variable)) +
  coord_flip()


# Discrete covariate, dependent -------------------------------------------

cov_dep = readRDS("data/2021-05-26_p=10_dependent_covariate.RDS")

cov_dep_long = melt(cov_dep, id.vars = c("individual", "simulation", "p"))

ggplot(cov_dep_long, aes(x = value)) +
  geom_boxplot() +
  facet_wrap(vars(variable)) +
  coord_flip()

