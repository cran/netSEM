#' Backsheet PET/PET/EVA Degradation
#' 
#' A backsheet is a mulit-layer polymer laminate used as the back surface of a crystalline silicon photovoltaic (PV) module that is designed to protect the inner components of the PV module.
#' Typical backsheets consist of three layers of polymers to control diffusion of water and oxygen, as well as to protect human beings from electrical shock.
#' It is critical that your PV module solar panel has a backsheet that is of high quality and can withstand various environmental elements for its 25 year expected product lifetime.
#' Poly(ethylene-terephthalate) (PET) is an important material, and it mainly used as the core (middle) layer of a backsheets to provide mechanical stability and electrical isolation. It is sometimes also used as the outer, or air side, layer of a backsheet. 
#' However, Backsheets using PET air side layers are highly susceptible to moisture and ultraviolet (UV) irradiance.
#' Degradation of PV module backsheets can cause severe safety issues and economic losses.
#' Indoor accelerated exposures are used to study backsheets degradation using higher irradiance and shorter exposure times so as to predict the performance of backsheets exposed to the real world conditions. 
#' 
#' PV module backsheet degradation can be characterized by changes in the  yellowness index (YI), which represents the color change of polymer and is associated with chemical change due to irradiance, high temperature and other process.
#' The YI value of a PV backsheet also relates to the module efficiency.
#' In addition, the Fourier-transform infrared spectroscopy (FTIR) is an effective method to study the chemical change of polymers. 
#'This non-destructive FTIR measurement enables researchers to obtain qualitative information about polymer functional groups and the relative contribution of each group in the sample to degradation.
#'
#' In this example, a dataset containing the YI of PET based backsheets (PET/PET/EVA) exposed to 1,500 hours of Damp Heat with no irradiance and a relative humidity (RH) of 85\% and a chamber temperature of 85 degree celsius. 
#'The PET samples were removed from the exposure chamber every 500 hours and the YI measurement was conducted at every time step.
#'The stress variable is exposure time with unit of hours.
#'Mechanistic variables from FTIR are included to track chemical changes in the materials related to polymer chain scission due to hydrolysis, crystallinity change and oxidation.
#' 
#' 
#' 
#' @docType data
#' @usage data(backsheet)
#' @name backsheet
#' @author Yu Wang, Wei-heng Huang, Laura S. Bruckman, Roger H. French
#' 
#' @format A data frame with 110 rows and 5 variables:
#' \describe{
#'   \item{YI}{Yellowness index of PET outer layer measured at every exposure step}
#'   \item{Hours}{Exposure time of PET in Damp Heat condition}
#'   \item{oxidation}{Formation of conjugated structures in the polymer during oxidation (FTIR peak at 1552 cm-1)}
#'   \item{hydrolysis}{Formation of acid and alcohol structures during hydrolysis (FTIR peak 2900 - 3300 cm-1)}
#'   \item{crystallization}{Crystallinity calculated as the ratio of FTIR peak at 1340 cm-1 to 1370 cm-1}
#' }
#' @source Solar Durability and Lifetime Extension (SDLE) Research Center, Case Western
#' Reserve University
"backsheet"