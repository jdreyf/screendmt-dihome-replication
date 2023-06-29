
safe_repfdr <- function(x1, x2){
  # fewer bins to avoid warning: glm.fit: fitted rates numerically 0 occurred
  # also number of bins here used for 3 statuses, instead of the usual 2
  nbins <- floor( sqrt(length(x1))/2 )
  # i'm seeing many warnings: In backend.ztobins(zmat, n.association.status, n.bins,  ... :
  # In column  2   x2 ,f(z) misfit =  5.6 . Rerun with increased df, but when I try df=20 from vignette get error in glm.fit, NA/NaN/Inf in 'x'
  z2b <- repfdr::ztobins(zmat = cbind(x1, x2), n.association.status = 3, n.bins = nbins)
  
  # from https://github.com/shay-y/repfdr/blob/master/vignettes/repfdr.Rmd
  # https://aosmith.rbind.io/2020/08/31/handling-errors/
  repfdr_safely <- purrr::safely(.f = repfdr::repfdr)
  res.repfdr = repfdr_safely(z2b$pdf.binned.z, z2b$binned.z.mat, non.null = "replication", 
                                               control = em.control(tol =  1e-4, verbose = FALSE))
  if (is.null(res.repfdr$error)){
    repfdr.q <- res.repfdr$result$mat[, "Fdr"]
  } else if (length(as.character(res.repfdr$error)) == 1 && 
             as.character(res.repfdr$error) == "Error in .f(...): Must have at least one study with estimated fraction of nulls below 1.\n"){
    repfdr.q <- rep(1, length(x1))
  } else {
    stop(as.character(res.repfdr$error))
  }
  repfdr.q
}
