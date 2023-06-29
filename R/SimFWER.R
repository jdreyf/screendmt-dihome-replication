#  Simulation study to assess adjusted p-value and power
#
#  Args:
#   mu... positive value representing the SNR of the false hypotheses
#   pi0... proportion of H_i such that both H_i1 and H_i2 are true
#   pi1... proportion of H_i such that exactly one of H_i1 and H_i2 is true
#   pi2... proportion of H_i such that both H_i1 and H_i2 are false
#   alpha... level at which adjusted p-value is to be controlled
#   m...  number of hypotheses
#   B... the number of Monte Carlo runs
#   fac... the ratio between the SNR in the two columns of the p-values
#   sseed... starting seed for reproducibilty
#   nblocks... number of blocks of analytes that covary
#   rho... correlation of analytes within blocks
#   em.max.iter... max number of EM iterations in repFdr
#   meth.v... Vector of method names

SimFWER <- function(mu, pi0, pi1, pi2, alpha, m, B = 100, fac=1, sseed = 123, nblocks = 10, rho=0.5, em.max.iter=10^3,
                   meth.v=c("Pearson", "DMT", "ScreenDMT", "radjust_sym", "AF")){
  stopifnot(abs(pi0+pi1+pi2-1) < 10**-3, meth.v %in% c("Pearson", "DMT", "ScreenDMT", "radjust_sym", "AF"),
            m %% nblocks == 0)
  
  c(mu1, mu2) %<-% sim_mu_vctrs(m=m, fac=fac, pi0=pi0, pi1=pi1, pi2=pi2)
  
  # adafilter block diagonal Sigma w/ b blocks & rho = 0.5
  sigma.tmp <- matrix(rho, ncol=m/nblocks, nrow=m/nblocks) + diag(rep(1-rho, m/nblocks))
  Sigma <- kronecker(X=diag(x = nblocks), Y=sigma.tmp)
  
  # true and false hypotheses
  TH <- mu1*mu2 == 0
  FH <- !TH
  # stopifnot(sum(FH) == m2)
  
  # results matrix
  power <- FP <- matrix(NA, ncol = length(meth.v), nrow = B)
  colnames(FP) <- colnames(power) <- meth.v
  
  # results are calculated column-wise, so bh's need not be in same order
  holm.mat <- matrix(NA, nrow=m, ncol=length(meth.v))
  colnames(holm.mat) <- meth.v
  
  set.seed(sseed)
  for (i in 1:B){
    if (i%%100 == 0) print(i)
    if (rho == 0){
      x1 <- rnorm(m, mean = mu1, sd = 1)
      x2 <- rnorm(m, mean = mu2, sd = 1)
    } else {
      x1 = mvrnormArma(1, mu1, Sigma)[1,]
      x2 = mvrnormArma(1, mu2, Sigma)[1,]
    }
    stopifnot(length(x1) == m, length(x2) == m)
    
    p.left.mat <- cbind(pnorm(x1, lower.tail=TRUE), pnorm(x2, lower.tail=TRUE))
    p.right.mat <- cbind(pnorm(x1, lower.tail=FALSE), pnorm(x2, lower.tail=FALSE))
    
    p1 <- 2*pnorm(abs(x1), lower.tail=FALSE)
    p2 <- 2*pnorm(abs(x2), lower.tail=FALSE)
    
    minp <- pmin(p1, p2)
    maxp <- pmax(p1, p2)
    
    if ("Pearson" %in% meth.v){
      p.left.max <- apply(p.left.mat, 1, FUN=max)
      p.right.max <- apply(p.right.mat, 1, FUN=max)
      p.pearson <- 2*pmin(p.left.max, p.right.max)
      holm.mat[, "Pearson"] <- p.adjust(p.pearson, method="holm")
    }
    
    # lotman: assume E -> Y effect > 0
    if ("DMT" %in% meth.v){
      lm.p <- ifelse(x1*x2 > 0, yes = maxp/2, no=1)
      holm.mat[, "DMT"] <- p.adjust(lm.p, method="holm")
    }
    
    # hitman2
    if ("ScreenDMT" %in% meth.v){
      tab.hm2 <- screendmt(tab = data.frame(x1, p1, x2, p2), p.adj.rate = "FWER")
      holm.mat[, "ScreenDMT"] <- tab.hm2$FWER
    }
    
    # adafilter: alpha should be independent of adjusted p
    if ("AF" %in% meth.v){
      af.left <- adaFilter(p.left.mat, r = 2, type.I.err = "FWER", alpha = alpha)
      af.right <- adaFilter(p.right.mat, r = 2, type.I.err = "FWER", alpha = alpha)
      holm.mat[, "AF"] <- 2*pmin(af.left$adjusted.p, af.right$adjusted.p)
    }
    
    rej <- apply(holm.mat, 2, function(x) x <= alpha)
    FP[i, ] <- apply(rej, 2, function(x) ifelse(sum(x) > 0, yes = sum(x[TH]), no = 0))
    if (sum(FH) == 0){
      power[i, ] <- rep(0, ncol(holm.mat))
    } else {
      power[i, ] <- apply(rej, 2, function(x) sum(x[FH]) / sum(FH))
    }
  } #end for
  
  FWER <- colMeans(FP >= 1)
  power <- colMeans(power)
  rbind(FWER, power)
}