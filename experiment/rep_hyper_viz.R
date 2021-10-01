library(ggplot2)
library(data.table)

# Discrete covariate, independent -----------------------------------------

metrics = readRDS("data/hyper_avg/discrete_independent.RDS")

metrics_long = melt(metrics, 
                    id.vars = c("method", "individual", 
                                "sim_idx", "covariate", "p", "weight1"))
metrics_long[method == "ss_ep2+varbvs", method := paste0(method, "_", weight1)]

plt = ggplot(metrics_long, aes(x = value, fill = method)) +
  geom_boxplot() +
  facet_wrap(vars(variable), scales = "free") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "right")





# No covariates -----------------------------------------------------------

metrics = readRDS("data/hyper_avg/no_covariate.RDS")

metrics_long = melt(metrics, 
                    id.vars = c("method", "individual", 
                                "sim_idx", "covariate", "p", "weight1"))
metrics_long[method == "ss_ep2+varbvs", method := paste0(method, "_", weight1)]

plt = ggplot(metrics_long, aes(x = value, fill = method)) +
  geom_boxplot() +
  facet_wrap(vars(variable), scales = "free") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "right")
