# comparing VB and EP on the same data

library(ggplot2)
library(data.table)

# discrete covariate, independent ----

# cov_indep = readRDS("data/discrete_independent/2021-06-23_covariate_independent.RDS")
files = dir("data/discrete_independent/", full.names = TRUE)
cov_indep = rbindlist(lapply(files, readRDS))

cov_indep_long = melt(cov_indep, 
                      id.vars = c("method", "individual", 
                                  "simulation", "covariate", "p", "w"))

cov_indep_long[method == "combo", method := paste(method, w, sep="_")]

plt = ggplot(cov_indep_long, aes(x = value, fill = method)) +
  geom_boxplot() +
  facet_wrap(vars(variable), scales = "free") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")

filename = paste0(Sys.Date(), "_discrete_independent_combo_boxplot.pdf")
ggsave(filename, plt, "pdf", path = "output", width = 6, height = 4)

# no covariate ----

files = dir("data/no_covariate/", full.names = TRUE)
no_cov = rbindlist(lapply(files, readRDS))

no_cov_long = melt(no_cov,
                   id.vars = c("method", "individual", 
                               "simulation", "covariate", "p", "w"))
no_cov_long[method == "combo", method := paste(method, w, sep="_")]

plt = ggplot(no_cov_long, aes(x = value, fill = method)) +
  geom_boxplot() +
  facet_wrap(vars(variable)) +
  coord_flip() +
  theme_bw() + 
  theme(legend.position = "bottom")
plt

filename = paste0(Sys.Date(), "_no_covariate_combo_boxplot.pdf")
ggsave(filename, plt, "pdf", path = "output", width = 6, height = 4)

  
# discrete covariate, dependent ----
files = dir("data/discrete_dependent", full.names = TRUE)

cov_dep = rbindlist(lapply(files, readRDS))

cov_dep_long = melt(cov_dep, 
                    id.vars = c("method", "individual", 
                                "simulation", "covariate", "p"))

plt = ggplot(cov_dep_long, aes(x = value, fill = method, color = covariate)) +
  geom_boxplot() + 
  facet_grid(cols = vars(p), rows = vars(variable), scales = "free_y") +
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "bottom")
plt

filename = paste0(Sys.Date(), "_discrete_dependent_boxplot.pdf")
ggsave(filename, plt, "pdf", path = "output", width = 6, height = 4)

# continuous covariate ----

n = 180
p = 4
n_sim = 50
MAXITER = 1
STR = 1
in_pr_13 = matrix(0, MAXITER, n)
in_pr_12 = in_pr_13
Var_cont = function(z) {
  pr = matrix(0, p+1, p+1)
  diag(pr) = 2
  #  pr[1,2]=STR*((z>0) && (z< .33)) + (STR - STR*((z-.33)/.33))*((z>0.33) && (z<0.66)) + (0)*((z>0.66) && (z<1))
  #  pr[1,3]=0*((z>0) && (z< .33)) + (STR*((z-.33)/.33))*((z>0.33) && (z<0.66)) + (STR)*((z>0.66) && (z<1))
  pr[2,3] = STR
  pr[1,2] = STR*((z>-1) && (z< -.33)) + (STR - STR*((z+.23)/.56)) * ((z>-0.23) && (z<0.33)) + (0)*((z>0.43) && (z<1))
  pr[1,3] = 0*((z>-1) && (z< -.33)) + (STR*((z+.23)/.56)) * ((z>-0.23) && (z<0.33)) + (STR)*((z>0.43) && (z<1))
  
  pr[2,1] = pr[1,2]
  pr[3,1] = pr[1,3]
  pr[3,2] = pr[2,3]
  
  
  Var = solve(pr)
  return(Var)
}

sensitivity_20 = matrix(0, MAXITER, 1)
specificity_20 = sensitivity_20
sensitivity_90 = sensitivity_20
specificity_90 = sensitivity_20
sensitivity_160 = sensitivity_20
specificity_160 = sensitivity_20

Z = c(seq(-0.99, -0.331, (-.331+.99)/59), 
      seq(-0.229,0.329,(.329+.229)/59),
      seq(0.431,.99,(.99-.431)/59))
# Z=seq(0.01,.99,.98/(n-1))
Z = matrix(Z, n, 1)
X = matrix(0, n, p+1)
for(i in 1:n) {
  X[i, ] = MASS::mvrnorm(1, rep(0, p+1), Var_cont(Z[i]))
}

precisions = rbindlist(
  lapply(Z, function(z) {
    cov_mat = Var_cont(z)
    prec_mat = solve(cov_mat)
    data.table(
      prec_12 = prec_mat[1, 2],
      prec_13 = prec_mat[1, 3]
    )
  })
)

files = dir("data/continuous/", full.names = TRUE)

cov_cont = rbindlist(lapply(files, readRDS))

# hack to add method labels
# cov_cont[, method := rep(rep(c("VB", "EP_fix", "EP_opt"), each = n), 50)]
# # hack to fix individual
# cov_cont[, individual := rep(rep(1:n, 3), 50)]
# # hack to fix simulation
# n_sim = 50
# cov_cont[, simulation := rep(1:n_sim, each = 3*n)]
cov_cont[, c("precision_12", "precision_13") := 
           list(rep(precisions$prec_12, 3),
                rep(precisions$prec_13, 3)), 
         by = "simulation"]


cont_cov_long = melt(cov_cont, id.vars = c("individual", "simulation",
                                           "incl_prob_12", "incl_prob_13",
                                           "covariate", "p", "method"), 
                     variable.name = "precision", value.name = "precision_value")
cont_cov_long = melt(cont_cov_long, 
                     id.vars = c("individual", "simulation", 
                                 "precision", "precision_value", "covariate", "p",
                                 "method"),
                     variable.name = "inclusion_prob", value.name = "inclusion_est")

summ_cont_cov = 
  cont_cov_long[, .(mean_incl_prob = mean(inclusion_est),
                    se_incl_prob = sd(inclusion_est) / sqrt(n_sim)),
                by = c("inclusion_prob", "covariate",
                       "individual", "method")]

plt = ggplot(summ_cont_cov, aes(x = individual, y = mean_incl_prob)) +
  geom_point() +
  geom_ribbon(aes(ymin = mean_incl_prob + 2*se_incl_prob,
                  ymax = mean_incl_prob - 2*se_incl_prob),
              alpha = .2) +
  facet_grid(rows = vars(inclusion_prob), cols = vars(method)) +
  
  theme_bw() +
  labs(y = "Inclusion Probability", x = "Subject Index")
  
filename = paste0(Sys.Date(), "_no_covariate_combo_boxplot.pdf")
ggsave(filename, plt, "pdf", path = "output", width = 6, height = 4)
