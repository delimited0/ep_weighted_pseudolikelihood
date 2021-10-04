library(ggplot2)
library(data.table)

method_levels = 
  c("ss_ep2", 
    "ss_ep2+varbvs_0.75", "ss_ep2+varbvs_0.5", "ss_ep2+varbvs_0.25",
    "varbvs")

# Discrete covariate, independent -----------------------------------------

metrics = readRDS("data/hyper_avg/discrete_independent.RDS")

metrics_long = melt(metrics, 
                    id.vars = c("method", "individual", 
                                "sim_idx", "covariate", "p", "weight1"))
metrics_long[method == "ss_ep2+varbvs", method := paste0(method, "_", weight1)]
metrics_long[, method := factor(method, levels = method_levels)]

plt = ggplot(metrics_long, aes(x = value, fill = method)) +
  geom_boxplot() +
  facet_wrap(vars(variable), scales = "free") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "right")
plt

filename = paste0(Sys.Date(), "_discrete_independent_boxplot.pdf")
ggsave(filename, plt, "pdf", path = "output/hyperparam_avg", width = 6, height = 4)


# No covariates -----------------------------------------------------------

metrics = readRDS("data/hyper_avg/no_covariate.RDS")

metrics_long = melt(metrics, 
                    id.vars = c("method", "individual", 
                                "sim_idx", "covariate", "p", "weight1"))
metrics_long[method == "ss_ep2+varbvs", method := paste0(method, "_", weight1)]
metrics_long[, method := factor(method, levels = method_levels)]

plt = ggplot(metrics_long, aes(x = value, fill = method)) +
  geom_boxplot() +
  facet_wrap(vars(variable), scales = "free") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "right")
plt

filename = paste0(Sys.Date(), "_no_covariate_boxplot.pdf")
ggsave(filename, plt, "pdf", path = "output/hyperparam_avg", width = 6, height = 4)


# Discrete Dependent ------------------------------------------------------
paths = dir("data/hyper_avg/", full.names = TRUE)
metrics = 
  rbindlist(lapply(paths[stringr::str_detect(paths, "discrete_dependent_p=*")],
                   function(x) readRDS(x)))

metrics_long = melt(metrics, 
                    id.vars = c("method", "individual", 
                                "sim_idx", "covariate", "p", "weight1"))
metrics_long[method == "ss_ep2+varbvs", method := paste0(method, "_", weight1)]
metrics_long[, method := factor(method, levels = method_levels)]

plt = ggplot(metrics_long, aes(x = value, fill = method, linetype = as.factor(covariate))) +
  geom_boxplot() +
  facet_grid(rows = vars(p), cols = vars(variable)) + 
  coord_flip() +
  theme_bw() +
  theme(legend.position = "right")
plt

filename = paste0(Sys.Date(), "_discrete_dependent_boxplot.pdf")
ggsave(filename, plt, "pdf", path = "output/hyperparam_avg", width = 6, height = 4)



