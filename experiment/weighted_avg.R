
# discrete independent ----

files = dir("data/discrete_independent_graph/", full.names = TRUE)
graphs = lapply(files, readRDS)

length(graphs) = n_sim

weights = seq(0, )

for (sim_idx in 1:n_sim) {
  
  vsvb_result = graphs[[sim_idx]]$vsvb
  ep_fix_result = graphs[[sim_idx]]$ep_fix
  
  vsvb_neg = vsvb_result$graphs[[1]]
  vsvb_pos = vsvb_result$graphs[[100]]
  
  ep_fix_neg = ep_fix_result$graphs[[1]]
  ep_fix_pos = ep_fix_result$graphs[[2]]
  
}


 = graphs[[1]]
foo = lapply(graphs, function(g) g$vsvb$graphs)


