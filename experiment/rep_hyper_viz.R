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

# Continuous covariate --------------------------------------------------------------

n = 180
p = 4

# ground truth graph (precision matrix)
Prec_cont = function(z) {
  
  STR = 1
  
  pr = matrix(0, p+1, p+1)
  diag(pr) = 2
  
  pr[2,3] = STR
  pr[1,2] = STR*((z>-1) && (z< -.33)) + (STR - STR*((z+.23)/.56)) * ((z>-0.23) && (z<0.33)) + (0)*((z>0.43) && (z<1))
  pr[1,3] = 0*((z>-1) && (z< -.33)) + (STR*((z+.23)/.56)) * ((z>-0.23) && (z<0.33)) + (STR)*((z>0.43) && (z<1))
  
  pr[2,1] = pr[1,2]
  pr[3,1] = pr[1,3]
  pr[3,2] = pr[2,3]
  
  return(pr)
}

Z = c(seq(-0.99, -0.331, (-.331+.99)/59), 
      seq(-0.229,0.329,(.329+.229)/59),
      seq(0.431,.99,(.99-.431)/59))
Z = matrix(Z, n, 1)

# read the full results
paths = dir("data/hyper_avg/continuous/", full.names = TRUE)
ep_paths = paths[stringr::str_detect(paths, "ep_")]
varbvs_paths = paths[stringr::str_detect(paths, "varbvs_")]

ep_results = lapply(ep_paths, readRDS)
varbvs_results = lapply(varbvs_paths, readRDS)

n_sim = length(ep_results)

metrics = rbindlist(
  mapply(function(ep_sim_results, varbvs_sim_results) {
    
    rbind(
      epwpl::incl_prob_123(ep_sim_results),
      epwpl::incl_prob_123(varbvs_sim_results),
      epwpl::incl_prob_123_combo(ep_sim_results, varbvs_sim_results, .25),
      epwpl::incl_prob_123_combo(ep_sim_results, varbvs_sim_results, .5),
      epwpl::incl_prob_123_combo(ep_sim_results, varbvs_sim_results, .75)
    )
  }, ep_results, varbvs_results, SIMPLIFY = FALSE)
)

precisions = rbindlist(
  lapply(Z, function(z) {
    prec_mat = Prec_cont(z)
    data.table(
      prec_12 = prec_mat[1, 2],
      prec_13 = prec_mat[1, 3]
    )
  })
)
precisions[, individual := 1:nrow(precisions)]

# metrics_long = melt(metrics, id.vars = c("individual", 
#                                          "incl_prob_12", "incl_prob_13",
#                                          "covariate", "p", "method"), 
#                     variable.name = "precision", value.name = "precision_value")

agg_metrics = 
  metrics[, .(mean_incl_prob_12 = mean(incl_prob_12),
              mean_incl_prob_13 = mean(incl_prob_13),
              q5_incl_prob_12 = quantile(incl_prob_12, .05),
              q95_incl_prob_12 = quantile(incl_prob_12, .95),
              q5_incl_prob_13 = quantile(incl_prob_13, .05),
              q95_incl_prob_13 = quantile(incl_prob_13, .95) 
  ), by = c("individual", "method")]

plt = ggplot(agg_metrics, aes(x = individual, y = mean_incl_prob_12)) +
  geom_point(shape = 20) +
  geom_ribbon(aes(ymin = q5_incl_prob_12, ymax = q95_incl_prob_12),
              alpha = .2) + 
  facet_wrap(vars(method)) +
  geom_path(data = precisions, aes(x=individual, y=prec_12)) +
  theme_bw()
filename = paste0(Sys.Date(), "_continuous_prec12.pdf")
ggsave(filename, plt, "pdf", path = "output/hyperparam_avg/", width = 6, height = 4)

plt = ggplot(agg_metrics, aes(x = individual, y = mean_incl_prob_13)) +
  geom_point(shape = 20) +
  geom_ribbon(aes(ymin = q5_incl_prob_13, ymax = q95_incl_prob_13),
              alpha = .2) + 
  facet_wrap(vars(method)) +
  geom_path(data = precisions, aes(x = individual, y = prec_13)) +
  theme_bw()
filename = paste0(Sys.Date(), "_continuous_prec13.pdf")
ggsave(filename, plt, "pdf", path = "output/hyperparam_avg/", width = 6, height = 4)

# why are the results so weird? Look at the raw data
raw_metrics = metrics[method %in% c("ss_ep2", "varbvs"), 
                      c("method", "individual", "incl_prob_12", "incl_prob_13")]
plt = ggplot(melt(raw_metrics, id.vars = c("method", "individual")),
             aes(x = individual, y = value)) +
  geom_point(pch = '.') +
  facet_grid(rows = vars(method), cols = vars(variable)) +
  theme_bw()
filename = paste0(Sys.Date(), "_continuous_diagnostic.pdf")
ggsave(filename, plt, "pdf", path = "output/hyperparam_avg/", width = 6, height = 4)

## prior inclusion diagnostic -------------------------------------------

precisions = rbindlist(
  lapply(Z, function(z) {
    prec_mat = Prec_cont(z)
    data.table(
      prec_12 = prec_mat[1, 2],
      prec_13 = prec_mat[1, 3]
    )
  })
)
precisions[, individual := 1:nrow(precisions)]

metrics = rbindlist(lapply(ep_results, function (ep_sim_results) {
  epwpl::alpha_123(ep_sim_results)
}))

metrics = rbindlist(lapply(1:length(varbvs_results), function (sim_idx) {
  vb_sim_results = varbvs_results[[sim_idx]]
  epwpl::alpha_123(vb_sim_results)[, sim_idx := sim_idx]
}))
metrics = melt(metrics[ , c("method", "individual",
                            "alpha_12", "alpha_13",
                            "logodds_12", "logodds_13", "sim_idx")],
               id.vars = c("method", "individual", "logodds_12", "logodds_13", "sim_idx"))
metrics[, c("variable", "element") := tstrsplit(variable, "_", fixed=TRUE)]

ggplot(data = precisions, aes(x=individual, y=prec_12)) +
  geom_path() +
# ggplot(metrics[element == 12 & plogis(logodds_12) == .5], 
#        aes(x = individual, y = value, color = logodds_12, fill = logodds_12)) +
  geom_point(data = metrics[element == 12 & plogis(logodds_12) == .4], 
             aes(x = individual, y = value, color = logodds_12, fill = logodds_12),
             alpha = .2)


metrics_12 = metrics[, c("alpha_12", "logodds_12", "individual")]

ggplot(metrics_12, aes(x = individual, y = alpha_12, color = logodds_12)) + 
  geom_hex()
  
  
