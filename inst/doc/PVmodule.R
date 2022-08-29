## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  
#  ## Load the pv module data set
#  data(PVmodule)
#  
#  ## Run netSEMp1 model
#  ans1 <- netSEMp1(PVmodule)
#  #plot network model with latent variable
#  plot(ans1, style = F, cutoff = c(0.3, 0.6, 0.9),
#       latent=c('EVA_hyd', 'Hac'))
#  
#  ## Run netSEMp2 model
#  ans2 <- netSEMp2(PVmodule)
#  ## Plot the network model for principle 2
#  plot(ans2, style = F, cutoff = c(0.3, 0.6, 0.9),
#       latent=c('EVA_hyd', 'Hac'))

## ---- out.width="675px", echo=FALSE, fig.cap="PVmodule netSEMp1 model"--------
knitr::include_graphics("PV1.png")

## ---- out.width="675px", echo=FALSE, fig.cap="PVmodule netSEMp2 model"--------
knitr::include_graphics("PV2.png")

