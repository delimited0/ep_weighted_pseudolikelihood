library(data.table)
source("experiment/simulation.R")
setDTthreads(1)
datadir = "data/pip_effect/"

# Discrete covariate, dependent -------------------------------------------

set.seed(1)
n = 100
n_sim = 50

pips = seq(.1, .9, .1)
sigma0 = 1
v_slab = 3 
tau = 1

for (p in c(10, 30, 50)) {
  
  Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
  Lam2 = c(rep(0, p-3), 3, 3, 3, 3) * 5
  
  Prec1 = Lam1 %*% t(Lam1) + diag(rep(10, p+1))
  Prec2 = Lam2 %*% t(Lam2) + diag(rep(10, p+1))
  
  Var1 = solve(Prec1)
  Var2 = solve(Prec2)
  
  Z = matrix(-.1*(1:n <= n/2)  + .1*(1:n > n/2), nrow = n, ncol = p, byrow = FALSE)
  
  weight_mat = weight_matrix(n, Z, tau)
  weight_mat_fit = weight_mat[c(1, n), ]
  
  true_graph_neg = Prec1 != 0
  diag(true_graph_neg) = 0
  true_graph_pos = Prec2 != 0
  diag(true_graph_pos) = 0
  
  # one data set per dimension
  X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
  X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
  data_mat = rbind(X1, X2)  

  for (p0 in pips) {
    
    print(paste0(
      "--- ", 
      "dim = ", p,
      ", p0 = ", p0 
    ))
     
    vsvb_result =
      wpl_vsvb_regression(data_mat, weight_mat,
                          sigma0, p0, v_slab, tune = TRUE, tune_p0 = FALSE)
    
    ep_result = 
      wpl_ep_regression(data_mat, weight_mat_fit, 
                        sigma0, p0, v_slab, opt = TRUE)
    
    metrics = rbind(
      score_model("VB", mean_symmetrize(vsvb_result$graphs[[1]]), true_graph_neg,
                  1, 1, -.1, p),
      score_model("VB", mean_symmetrize(vsvb_result$graphs[[n]]), true_graph_pos,
                  2, 1, .1, p),
      score_model("EP_opt", mean_symmetrize(ep_result$graphs[[1]]), true_graph_neg, 
                  1, 1, -.1, p),
      score_model("EP_opt", mean_symmetrize(ep_result$graphs[[2]]), true_graph_pos, 
                  2, 1, .1, p)
    )
    metrics[, p0 := p0
            ][, mlik := c(rep(vsvb_result$elbo, 2), 
                          rep(ep_result$llik, 2))]
    
    filename = paste0(datadir,
                      Sys.Date(),
                      "_dim=", p, 
                      "_p0=", p0,
                      "_dependent_covariate.RDS")
    saveRDS(metrics, file = filename)
  }
}


# Analyze results ---------------------------------------------------------

library(data.table)

files = dir("data/pip_effect/", full.names = TRUE)
pip_dt = rbindlist(lapply(files, readRDS))
# pip_dt_long = 
#   melt(pip_dt, 
#        id.vars = c("method", "individual", "covariate", "p", "p0", "mlik"),
#        measure.vars = list(metric = c("sensitivity", "specificity")))
pip_dt_long = melt(pip_dt,
                   id.vars = c("method", "individual", "covariate", "p", "p0"),
                   measure.vars = "mlik")

library(ggplot2)

# how does likelihood change with pip

ggplot(pip_dt, aes(x = p0, y = mlik)) +
  geom_point() +
  facet_grid(rows = vars(p), cols = vars(method), scales = "free")

ep_plot = ggplot(pip_dt[method == "EP_opt"], aes(x=p0, y=mlik)) +
  geom_point() + geom_path() +
  facet_grid(rows = vars(p), scales="free_y") #, nrow = 3)

filename = paste0(Sys.Date(), "_ep_mlik_vs_p0.pdf")
ggsave(filename, ep_plot, "pdf", width = 6, height=8, path = "output")

vb_plot = ggplot(pip_dt[method == "VB"], aes(x=p0, y=mlik)) +
  geom_point() + geom_path() +
  facet_grid(rows = vars(p), scales="free_y") #, nrow = 3)

filename = paste0(Sys.Date(), "_vb_elbo_vs_p0.pdf")
ggsave(filename, vb_plot, "pdf", width = 6, height=8, path = "output")

# matching sensitivity/specificity
pip_acc_dt = melt(pip_dt,
                  id.vars = c("method", "individual", "covariate", "p", "p0"),
                  measure.vars = list(sens = "sensitivity", 
                                      spec = "specificity"))
acc_plot = ggplot(pip_acc_dt, aes(x = p0, y = sens, 
                                  shape = as.factor(covariate))) +
  geom_point(color = "blue") +
  geom_point(aes(x = p0, y = spec), color = "red") +
  facet_grid(rows = vars(p), cols = vars(method)) +
  theme_bw() +
  theme(legend.position = "bottom")

filename = paste0(Sys.Date(), "_accuracy_vs_p0.pdf")
ggsave(filename, acc_plot, "pdf", width=6, height=8, path="output")
