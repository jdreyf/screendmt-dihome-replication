
sim_ggplot <- function(pwr_arr, meth.first=c("ScreenDMT", "DMT")){
  gg.df <- melt(pwr_arr) |> 
    dplyr::mutate(pi1 = as.numeric(gsub("pi1_", "", x=Var1))) |> 
    dplyr::mutate(M = as.character(gsub("m_", "M=", x=Var2))) |>
    dplyr::mutate(SNR = factor(gsub("fac_1", "Equal signal", x=gsub("fac_2", "Unequal signal (2x)", x=Var3)))) |>
    dplyr::mutate(Type = sub("power", "Power", as.character(Var4))) |>
    dplyr::mutate(Method = sub("AF", "AdaFilter", 
                               sub("radjust_sym", "radjust-sym", 
                                   sub("repfdr", "RepFdr", as.character(Var5))))) |>
    dplyr::mutate(Method = factor(Method, levels = c(meth.first, sort(setdiff(Method, meth.first))))) |> 
    dplyr::mutate(nb = factor(sub("nb_10", "Strong dependence", 
                           sub("nb_100", "Weak dependence", x=Var6)), levels=c("Weak dependence", "Strong dependence"))) |>
    dplyr::select(!Var1:Var6)
  
  type1 <- setdiff(gg.df$Type, "Power")
  
  ggp.lst <- list()
  for (Type.tmp in unique(gg.df$Type)[2:1]){
    for (M.tmp in unique(gg.df$M)){
      ggp <- gg.df |> filter(Type == Type.tmp & M == M.tmp) |> 
        ggplot() + geom_line(mapping = aes(pi1, value, color=Method, linetype=Method)) + facet_grid(cols=vars(SNR), rows = vars(nb)) + 
        labs(x="Percent of null analytes with signal", y=Type.tmp) + theme_bw() + 
        scale_x_continuous(trans=scales::pseudo_log_trans(base = 10, sigma=0.01), breaks = pi1.v, labels = 100*pi1.v) + 
        geom_hline(yintercept=alpha, linetype="dotted")
      
      if (Type.tmp == "FDR"){
        ggp <- ggp + labs(title = "False discovery rate")
      } else if (Type.tmp == "FWER"){
        ggp <- ggp + labs(title = "Familywise error rate")
      } else {
        ggp <- ggp + labs(title = "Probability of detecting true associations")
      }
      # want to have corresponding colors in both PDFs
      if (type1 == "FWER") ggp <- ggp + scale_color_manual(values = scales::hue_pal()(6)[1:4])
      
      ggp.lst[[ Type.tmp ]] <- ggp
    } # end for M.tmp
  }
  ggp.lst
}