library(data.table)
library(microbenchmark)

source("./R/metrics.R")
source("./R/dataprep.R")
source("./R/functions.R")
extdata.path <- system.file("extdata",package = "TargetedMSQC")
project.folder.name <- "CSF_Panel"
project.path <- file.path(extdata.path,project.folder.name)
chromatogram.path <- file.path(project.path,"Chromatograms")
peak.boundary.path <- file.path(project.path,"Peak_boundary")
d <- readRDS('test.rds')
load("~/Documents/GitHub/TargetedMSQC/data/data.CSF.rda")
res1 <- microbenchmark(
  ExtractFeatures(
    data = copy(d$data),
    export.features = FALSE,
    intensity.threshold = 1000
  ),
  TargetedMSQC::ExtractFeatures(data = data.CSF$data,
                                export.features = FALSE,
                                intensity.threshold = 1000),
  times = 10
)
saveRDS(res1,'result_benchmark.rds')

