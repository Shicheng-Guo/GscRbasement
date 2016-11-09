
# For two vectors: p & q are distributions so their elements should sum up to 1
p <- c(0.00029421, 0.42837957, 0.1371827, 0.00029419, 0.00029419,0.40526004, 0.02741252, 0.00029422, 0.00029417, 0.00029418)
q <- c(0.00476199, 0.004762, 0.004762, 0.00476202, 0.95714168,0.00476213, 0.00476212, 0.00476202, 0.00476202, 0.00476202)
m <- 0.5 * (p + q)
JS <- 0.5 * (sum(p * log(p / m)) + sum(q * log(q / m)))

# For more than 2 vectors: 
H <- function(v) {
      v <- v[v > 0]
      return(sum(-v * log(v)))
}
JSD <- function(w, m) {
  return(H(m %*% w) - apply(m, 2, H) %*% w)
}
JSD(w = c(1/3, 1/3, 1/3), m = cbind(p, q, m))

