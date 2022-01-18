#' vectorized logsumexp function
#' element wise computation of logsumexp
log_sum_exp = function(u, v) {
  maxuv = pmax(u, v)
  maxuv + log(exp(u - maxuv) + exp(v - maxuv))
}

log1mexp = function(a) {
  if (a < log(.5))
    log1p(-exp(a))
  else 
    log(-expm1(a))
}
# 
# log_subtract = function(x, y) {
#   if (x <= y)
#     
# }

