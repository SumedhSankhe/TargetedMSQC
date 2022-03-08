## ----setup, include = FALSE---------------------------------------------------
# knitr::opts_chunk$set(
#   collapse = TRUE,
#   comment = "#>"
# )

## ----install, eval = FALSE, echo=TRUE,message=FALSE---------------------------
#  library(devtools)
#  devtools::install_github("shadieshghi/TargetedMSQC",build_vignettes = TRUE)

## ----library, echo = TRUE,message=FALSE---------------------------------------
# library(TargetedMSQC)

## ----CleanUpChromatograms, eval = FALSE---------------------------------------
#  # Set the path to the subdirectories of chromatogram and peak boundary files
#  extdata.path <- system.file("extdata",package = "TargetedMSQC")
#  project.folder.name <- "CSF_Panel"
#  project.path <- file.path(extdata.path,project.folder.name)
#  chromatogram.path <- file.path(project.path,"Chromatograms")
#  peak.boundary.path <- file.path(project.path,"Peak_boundary")
#
#  # CleanUpChromatograms reformats the Skyline exported files into a data frame
#  data.CSF <- CleanUpChromatograms(chromatogram.path = chromatogram.path,
#                               peak.boundary.path = peak.boundary.path,
#                               endogenous.label = "light",
#                               standard.label = "heavy")
#

