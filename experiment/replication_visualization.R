library(ggplot2)
library(data.table)

# Discrete covariate, independent -----------------------------------------

# legacy plot

cov_indep = readRDS("data/2021-05-29_covariate_independent.RDS")

n = 100
covariate = -.1*(1:n <= n/2)  + .1*(1:n > n/2)
cov_indep[, c("covariate", "individual") := 
            list(as.factor(covariate), 1:n), 
          by = "simulation"]

cov_indep_long = melt(cov_indep, id.vars = c("individual", "simulation", "covariate"))

summ_cov_indep = cov_indep_long[, .(mean_value = mean(value)),
                                by = c("individual", "variable", "covariate")]

plt = ggplot(summ_cov_indep, aes(x = mean_value, fill = covariate)) +
  geom_boxplot() +
  facet_wrap(vars(variable)) +
  coord_flip() +
  theme_bw()

filename = paste0(Sys.Date(), "_discrete_covariate_independent_boxplot.pdf")
ggsave(filename, plt, "pdf", width = 6, height = 4, path = "output")

# by simulation plot

cov_indep = readRDS("data/2021-05-31_covariate_independent.RDS")
cov_indep_long = melt(cov_indep, id.vars = c("individual", "simulation", "covariate"))

plt = ggplot(cov_indep_long, aes(x = value, fill = as.factor(covariate))) +
  geom_boxplot() + 
  facet_wrap(vars(variable)) +
  coord_flip() +
  theme_bw() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom") +
  labs(fill = "Covariate")

filename = paste0(Sys.Date(), "_discrete_covariate_independent_boxplot.pdf")
ggsave(filename, plt, "pdf", width = 6, height = 4, path = "output")

# No covariate ------------------------------------------------------------

# legacy plot

no_cov = readRDS("data/2021-05-29_no_covariate.RDS")

n = 100
covariate = -.1*(1:n <= n/2)  + .1*(1:n > n/2)
no_cov[, c("covariate", "individual") := 
         list(as.factor(covariate), 1:n), 
       by = "simulation"]

no_cov_long = melt(no_cov, id.vars = c("individual", "simulation", "covariate"))

summ_no_cov = no_cov_long[, .(mean_value = mean(value)),
                             by = c("individual", "variable", "covariate")]

plt = ggplot(summ_no_cov, aes(x = mean_value, fill = covariate)) +
  geom_boxplot() +
  facet_wrap(vars(variable)) +
  coord_flip() +
  theme_bw()

filename = paste0(Sys.Date(), "_no_covariate_boxplot.pdf")
ggsave(filename, plt, "pdf", width = 6, height = 4, path = "output")

# by simulation plot

no_cov = readRDS("data/2021-05-31_no_covariate.RDS")
no_cov_long = melt(no_cov, id.vars = c("individual", "simulation", "covariate"))

plt = ggplot(no_cov_long, aes(x = value, fill = as.factor(covariate))) +
  geom_boxplot() + 
  facet_wrap(vars(variable)) +
  coord_flip() +
  theme_bw() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom") +
  labs(fill = "Covariate")

filename = paste0(Sys.Date(), "_no_covariate_boxplot.pdf")
ggsave(filename, plt, "pdf", width = 6, height = 4, path = "output")

# Discrete covariate, dependent -------------------------------------------

# averaged by individual plot

all_simulations = dir("data", full.names = TRUE)
cov_dep_filenames = 
  all_simulations[stringr::str_detect(dir("data"), "_dependent_covariate")]

cov_dep_sims = lapply(cov_dep_filenames, readRDS) 
n = 100
covariate = c(rep(-.1, n/2), rep(.1, n/2))
cov_dep_sims = lapply(cov_dep_sims, function(sim) {
  sim[, c("individual", "covariate") := list(1:n, covariate), by = simulation]
})


cov_dep_dt = rbindlist(cov_dep_sims)

avg_cov_dep_dt = cov_dep_dt[, .(sensitivity = mean(sensitivity), 
                                specificity = mean(specificity),
                                covariate = as.factor(covariate)),
                            by = c("individual", "p")]

avg_cov_dep_long = melt(avg_cov_dep_dt, id.vars = c("individual", "p", "covariate"))

plt = ggplot(avg_cov_dep_long, aes(x = value, fill = covariate)) +
  geom_boxplot() +
  facet_grid(rows = vars(variable), cols = vars(p)) +
  coord_flip() +
  theme_bw()

filename = paste0(Sys.Date(), "_discrete_covariate_dependent_boxplot.pdf")
ggsave(filename, plt, "pdf", width = 8, height = 6, path = "output")

# by simulation plot

all_simulations = dir("data", full.names = TRUE)
cov_dep_filenames = 
  all_simulations[stringr::str_detect(all_simulations, 
                                      "2021-05-31_p=\\d\\d_dependent_covariate")]

cov_dep_sims = lapply(cov_dep_filenames, readRDS) 
cov_dep_dt = rbindlist(cov_dep_sims)
cov_dep_dt[, p := rep(c(10, 30, 50), each = 100)]  # forgot to add dimension

cov_dep_long = melt(cov_dep_dt, id.vars = c("individual", "p", "covariate", "simulation"))

plt = ggplot(cov_dep_long, aes(x = value, fill = as.factor(covariate))) +
  geom_boxplot() +
  facet_grid(rows = vars(variable), cols = vars(p)) +
  coord_flip() + 
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom") +
  labs(fill = "Covariate")

filename = paste0(Sys.Date(), "_discrete_covariate_dependent_boxplot.pdf")
ggsave(filename, plt, "pdf", width = 8, height = 6, path = "output")

# Continuous covariate ----------------------------------------------------

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

cont_cov = readRDS("data/2021-05-29_continuous_covariate.RDS")
cont_cov[, covariate := Z, by = "simulation"]
cont_cov[, c("precision_12", "precision_13") := 
           list(precisions$prec_12, precisions$prec_13), by = "simulation"]

cont_cov_long = melt(cont_cov, id.vars = c("individual", "simulation",
                                           "incl_prob_12", "incl_prob_13",
                                           "covariate"), 
                     variable.name = "precision", value.name = "precision_value")
cont_cov_long = melt(cont_cov_long, 
                     id.vars = c("individual", "simulation", 
                                 "precision", "precision_value", "covariate"),
                     variable.name = "inclusion_prob", value.name = "inclusion_est")

summ_cont_cov = 
  cont_cov_long[, .(mean_incl_prob = mean(inclusion_est),
                    se_incl_prob = sd(inclusion_est) / sqrt(n_sim)),
                by = c("inclusion_prob", "covariate",
                       "individual")]

plt = ggplot(summ_cont_cov, aes(x = individual, y = mean_incl_prob)) +
  geom_point() +
  geom_ribbon(aes(ymin = mean_incl_prob + 2*se_incl_prob,
                  ymax = mean_incl_prob - 2*se_incl_prob),
              alpha = .2) +
  facet_grid(rows = vars(inclusion_prob)) +
  theme_bw() +
  labs(y = "Inclusion Probability", x = "Subject Index")

filename = paste0(Sys.Date(), "_continuous_covariate_by_subject.pdf")
ggsave(filename, plt, "pdf", width = 6, height = 4, path = "output")

ggplot(cont_cov_long, aes(x = individual, y = precision_value)) +
  geom_point() + 
  facet_grid(rows = vars(precision)) +
  theme_bw()
