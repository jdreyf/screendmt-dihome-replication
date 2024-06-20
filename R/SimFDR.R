#  Simulation study to assess FDR and power
#
#  Args:
#   mu... positive value representing the SNR of the false hypotheses
#   pi0... proportion of null H_i such that both component nulls H_i1 and H_i2 are true
#   pi1... proportion of null H_i such that exactly one of H_i1 and H_i2 is true
#   pi2... proportion of null H_i such that both H_i1 and H_i2 are false
#   alpha... level at which FDR is to be controlled
#   m...  number of hypotheses
#   B... the number of simulations
#   fac... the ratio between the SNR in the two columns of the p-values
#   sseed... starting seed for reproducibilty
#   nblocks... number of blocks of analytes that covary
#   rho... correlation of analytes within blocks
#   em.max.iter... max number of EM iterations in repFdr
#   meth.v... Vector of method names

SimFDR <- function(mu, pi0, pi1, pi2, alpha, m, B = 100, fac=1, sseed = 123, nblocks = 10, rho=0.5, em.max.iter=10^3,
                   meth.v=c("JS", "Pearson", "DMT", "ScreenDMT", "radjust_sym", "AF")){
  stopifnot(abs(pi0+pi1+pi2-1) < 10**-3, meth.v %in% c("JS", "Pearson", "DMT", "ScreenDMT", "radjust_sym", "AF", "repfdr"),
            m %% nblocks == 0)
  
  # uses zeallot; mu's are means of stats
  c(mu1, mu2) %<-% sim_mu_vctrs(m=m, fac=fac, pi0=pi0, pi1=pi1, pi2=pi2, mu=mu)

  # adafilter block diagonal Sigma w/ b blocks & rho = 0.5
  sigma.tmp <- matrix(rho, ncol=m/nblocks, nrow=m/nblocks) + diag(rep(1-rho, m/nblocks))
  Sigma <- kronecker(X=diag(x = nblocks), Y=sigma.tmp)
  
  # true and false null hypotheses
  TH <- mu1*mu2 <= 0
  FH <- !TH
  
  # results matrix for power and false discovery proportion
  power <- FDP <- matrix(NA, ncol = length(meth.v), nrow = B)
  colnames(FDP) <- colnames(power) <- meth.v
  
  # matrix of BH FDRs
  # results are calculated column-wise, so bh's need not be in same order
  bh.mat <- matrix(NA, nrow=m, ncol=length(meth.v))
  colnames(bh.mat) <- meth.v
  
  set.seed(sseed)
  # Loop through simulations that each generate a dataset
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
    
    # for Pearson method
    p.left.mat <- cbind(pnorm(x1, lower.tail=TRUE), pnorm(x2, lower.tail=TRUE))
    p.right.mat <- cbind(pnorm(x1, lower.tail=FALSE), pnorm(x2, lower.tail=FALSE))
    
    p1 <- 2*pnorm(abs(x1), lower.tail=FALSE)
    p2 <- 2*pnorm(abs(x2), lower.tail=FALSE)
    
    minp <- pmin(p1, p2)
    maxp <- pmax(p1, p2)
    
    # js
    if ("JS" %in% meth.v) bh.mat[, "JS"] <- p.adjust(p=maxp, method = "BH")
    
    if ("Pearson" %in% meth.v){
      p.left.max <- apply(p.left.mat, 1, FUN=max)
      p.right.max <- apply(p.right.mat, 1, FUN=max)
      p.pearson <- 2*pmin(p.left.max, p.right.max)
      bh.mat[, "Pearson"] <- p.adjust(p.pearson, method="BH")
    }
    
    if ("DMT" %in% meth.v){
      lm.p <- ifelse(x1*x2 > 0, yes = maxp/2, no=1)
      bh.mat[, "DMT"] <- p.adjust(lm.p, method="BH")
    }
    
    if ("ScreenDMT" %in% meth.v){
      tab.hm2 <- screendmt(tab = data.frame(x1, p1, x2, p2), p.adj.rate = "FDR")
      bh.mat[, "ScreenDMT"] <- tab.hm2$FDR
    }
      
    # Heller 2018 radjust: paper suggests lambda=alpha, but this yields pi0 > 1 bug that resets lambda=0
    # ph <- MultiMed::medTest.SBMH(p1, p2, MCP.type="FDR", t1=alpha/2, t2=alpha/2, lambda=0)
    if ("radjust_sym" %in% meth.v){
      hell.bh <- rep(1, length(p1))
      pv1 <- ifelse(x1 > 0, p1/2, 1-p1/2)
      pv2 <- ifelse(x2 > 0, p2/2, 1-p2/2)
      hell <- suppressMessages(radjust_sym(pv1, pv2, input_type = "all_features", directional_rep_claim = TRUE, variant = "adaptive", 
                                           alpha=alpha)$results_tab)
      star.ind <- which(hell$Significant == "*")
      if (length(star.ind) > 0){
        hell.bh[ as.numeric(hell$name[star.ind]) ] <- 0
      }
      bh.mat[, "radjust_sym"] <- hell.bh
    }
    
    # repfdr: took minutes for 10 iterations of a single parameter combination!
    # this code addresses below error
    # Error in repfdr::repfdr(z2b$pdf.binned.z, z2b$binned.z.mat, non.null = "replication",  : 
    # Must have at least one study with estimated fraction of nulls below 1.
    # Error in if (p0theo >= 1) { : missing value where TRUE/FALSE needed
    if ("repfdr" %in% meth.v){
      bh.mat[, "repfdr"] <- safe_repfdr(x1, x2)
    }

    # adafilter: alpha should be independent of adjusted p
    # af <- adaFilter::adaFilter(cbind(p1, p2), r = 2, type.I.err = "FDR", alpha = alpha)$adjusted.p
    if ("AF" %in% meth.v){
      af.left <- adaFilter(p.left.mat, r = 2, type.I.err = "FDR", alpha = alpha)
      af.right <- adaFilter(p.right.mat, r = 2, type.I.err = "FDR", alpha = alpha)
      # there's no actual FDR for both tests, but we reject if either is below alpha/2 in SimFDR.R
      bh.mat[, "AF"] <- 2*pmin(af.left$adjusted.p, af.right$adjusted.p)
    }
    
    rej <- apply(bh.mat, 2, function(x) x <= alpha)
    FDP[i, ] <- apply(rej, 2, function(x) ifelse(sum(x) > 0, yes = sum(x[TH]) / sum(x), no = 0))
    if (sum(FH) == 0){
      power[i, ] <- rep(0, ncol(bh.mat))
    } else {
      power[i, ] <- apply(rej, 2, function(x) sum(x[FH]) / sum(FH))
    }
  } #end for
  
  FDR <- colMeans(FDP)
  power <- colMeans(power)
  rbind(FDR, power)
}