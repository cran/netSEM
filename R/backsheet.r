#' Backsheet PET/PET/EVA Degradation
#' 
#' A backsheet is a polymeric cover of photovoltaic (PV) module and is designed to protect the inner components of module.
#' Typical backsheets consist of three layers of polymers to prohibit diffusion of water and oxygen, as well as to protect human beings from electrical shock.
#' It is critical that your solar panel has a backsheet that is of high quality and can withstand various environmental elements for 25 long years.
#' Polyethylene terephthalate (PET) is an important material, and it mainly used as core and outer layer of backsheets and provide mechanical stability and electrical isolation.
#' However, PET based Backsheets are highly susceptible to moisture and ultraviolet (UV) irradiance.
#' Degradation of backsheets will causes severe economic loss and safty issue.
#' Indoor accelerated exposures are used to study backsheets degradation within short time and predict the performance of backsheets exposed to the real world.
#' 
#' The backsheet degradation can be characterized with yellowness index (YI), which represents the color change of polymer and is associated with chemical change due to irradiance, high temperature and other process.
#' The YI value of a PV backsheet also relates to the module efficiency.
#'In addition, the Fourier-transform infrared spectroscopy (FTIR) is a effective method to study the chemical change of polymers.
#'This non-destructive measurement enable people to obtain qualitative information of polymer functional groups and the relative amount of each group in the sample.
#'
#' In this example, a dataset containing the YI of PET based backsheets (PET/PET/EVA) exposed to 1,500 hours of Damp Heat with no irradiance and an relative humidity of 85\% and a chamber temperature of 85 degree celsius. 
#'The PET samples were pull out from exposure chamber every 500 hours and the YI measurement was conducted at every step.
#'The stress variable is exposure time with unit of hour.
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
#'   \item{oxidation}{Formed conjugation structure during oxidation (FTIR peak at 1552 cm-1)}
#'   \item{hydrolysis}{Fromed acid and alcohol structure during hydrolysis (FTIR peak 2900 - 3300 cm-1)}
#'   \item{crystallization}{Crystallinity calculated as the ratio of FTIR peak at 1340 cm-1 to 1370 cm-1}
#' }
#' @source Solar Durability and Lifetime Extension (SDLE) Research Center, Case Western
#' Reserve University
"backsheet"