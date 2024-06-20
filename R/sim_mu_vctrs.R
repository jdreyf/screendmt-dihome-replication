# Simulate vector of means
sim_mu_vctrs <- function(m, fac, pi0, pi1, pi2, mu){
  # true parameters
  m0 <- round(m*(pi0))
  m1 <- round(m*(pi1))
  # m2 <- round(m*pi2)
  m2 <- m-m0-m1
  
  # if unequal SNR, want same avg mu
  # (1+fac)x/2 = mu
  if (fac != 1) mu <- 2*mu/(1+fac)
  
  # randomly allocate m0, m1, & m2 rows
  mu1 <- mu2 <- numeric(m)
  
  # idx = index
  m1.idx <- sample(x=m, size=m1)
  m2.idx <- sample(x=m, size=m2)
  
  mu1[c(m1.idx, m2.idx)] <- sample(x=c(-mu, mu), size = m1+m2, replace = TRUE)
  
  m1.opp.idx <- if (length(m1.idx) == 1) m1.idx else sample(x=m1.idx, size=round(m1/2))
  if (length(m1.opp.idx) >= 1) mu2[m1.opp.idx] <- -1*mu1[m1.opp.idx]
  mu2[m2.idx] <- mu1[m2.idx]*fac # unequal signal
  
  stopifnot(length(mu1) == m, length(mu2) == m, sum(sign(mu1)*sign(mu2) == 1) == m2)
  list(mu1, mu2)
}