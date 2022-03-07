#' This file contains the functions that clean up the input data and calculate
#' the features associated with peaks and peak groups.

#' Convert the exported skyline chromatogram and peak boundary reports to a
#' dataframe format compatible with TargetedMSQC functions.
#'
#' The function takes the directories that contain chromatogram and peak boundary
#' files of Skyline documents and merges them. The rows with NA values for peak
#' boundaries, rows corresponding to peptides with fewer than 3 transitions and
#' transitions missing a standard isotope pair are separated into a dataframe
#' that is returned by the function in output$removed and excluded from subsequence
#' QC steps. The remainder of the rows are returned by the function in output$data
#' for subsequence analysis. Transition peaks of each peptide in each run are
#' formed into peak groups of custom class peakObj and stored in the ChromGroup
#' column. The ApplyPeakBoundary function is applied to chromatograms to restrict
#' them to the provided peak boundaries and the results are stored in the PeakGroup
#' column. Additionally, rows are added for the "sum" transition, which are is
#' the sum of peptides transitions with the same isotope label in each sample.
#'
#' @param chromatogram.path Path to the directory containing the .tsv files of
#' the peak chromatograms. For each Skyline document, this file is exported from
#' Skyline through File > Export > Chromatograms. Here, check runs of interest
#' and include Precursors, Products, Base Peaks and TICs. Each chromatogram .tsv
#' file corresponds to a single Skyline document, which may contain any number
#' of runs. Multiple chromatogram files, corresponding to multiple Skyline documents
#' can be copied into the chromatogram.path directory. For each chromatogram file
#' in this folder, there should be a peak boundary file with an identical name
#' in peak.boundary.path directory.
#' @param peak.boundary.path  Path to the directory containing the .csv files of
#' the peak boundaries. For each Skyline document, this file is exported from
#' Skyline through File > Export > Report. Here, select Peak Boundaries. Each
#' peak boundary .csv file corresponds to a single Skyline document, which may
#' contain any number of runs. Multiple peak boundary files, corresponding to
#' multiple Skyline documents can be copied into the peak.boundary.path directory.
#' For each peak boundary file in this folder, there should be a peak chromatogram
#' file with an identical name in chromatogram.path directory.
#' @param labkey.url.base URL of the labkey server, if data is being imported
#' from labkey (labkey = TRUE). This feature has not been implemented yet.
#' @param labkey.url.path Path to the directory containing the Skyline documents
#' on the labkey server, if data is being imported from labkey (labkey = TRUE).
#' This feature has not been implemented yet.
#' @param labkey Logical parameter to indicate whether data will be imported from
#' labkey. This feature has not been implemented yet. Set labkey = FALSE.
#' @param endogenous.label Label of the endogenous analyte: default is "light"
#' @param standard.label Label of the spiked-in isotopically labeled analyte: default is "heavy"
#' @param iRT.list List of iRT standards used in the experiment. These peptides
#' will be removed from the training set.
#'
#' @return A list with the following objects:
#'               data: A dataframe that contains all rows of the input data with
#'               assigned peak boundaries. This dataframe also includes columns
#'               for peakObj objects created for each peak group.
#'               removed: A dataframe that contains the rows with missing peak
#'               boundary values, too few transitions and missing isotope pairs.
#'               These rows are removed from downstream feature extraction and
#'               QC analysis.
#'
#' @export
#'
#' @import data.table
#'
#' @examples
#'
#' extdata.path <- system.file("extdata",package = "TargetedMSQC")
#' project.folder.name <- "CSF_Panel"
#' project.path <- file.path(extdata.path,project.folder.name)
#' chromatogram.path <- file.path(project.path,"Chromatograms")
#' peak.boundary.path <- file.path(project.path,"Peak_boundary")
#' data <- CleanUpChromatograms(chromatogram.path = chromatogram.path,
#'                             peak.boundary.path = peak.boundary.path,
#'                             endogenous.label = "light",
#'                             standard.label = "heavy",
#'                             iRT.list = iRTList())

CleanUpChromatograms <- function(chromatogram.path = NULL,
                                 peak.boundary.path = NULL,
                                 labkey.url.base = NULL,
                                 labkey.url.path = NULL,
                                 labkey = FALSE,
                                 endogenous.label = "light",
                                 standard.label = "heavy" ,
                                 iRT.list = iRTList(), ... ) {

  if (labkey) {
   errorReporting("Labkey worklfow not implemented")
  }else{
    # chromgram.path and peak.boundary.path can only be NULL if the code is being run on labkey.
    if (is.null(chromatogram.path)){
      errorReporting("Input chromatogram file should not be empty")
    }
    if (is.null(peak.boundary.path)){
      errorReporting("Input peak boundary file should not be empty")
    }

    # load data from the provided chromatogram and peak boundary paths
    input.data <- loadData(chrom.Path = chromatogram.path,
                           peak.Path = peak.boundary.path)
    # combine the chromatogram and peak boundary data in the right format
    chrom.data <- combineChromPeak(chrom = input.data$chrom,
                                   peak = input.data$peak.boundary,
                                   lvls = c(endogenous.label, standard.label))
    rm(input.data)
    gc()
  }

  # Rows where peak boundaries are NA are separated into a data frame and returned as output.
  removed <- chrom.data[Remove == T | PeptideModifiedSequence %in% iRT.list,.
                        (File, FileName, PeptideModifiedSequence, PrecursorCharge,
                          FragmentIon, ProductCharge, IsotopeLabelType, MinStartTime,
                          MaxEndTime)]
  chrom.data <- chrom.data[Remove != T & !PeptideModifiedSequence %in% iRT.list,
                           .(FileName, PeptideModifiedSequence, PrecursorCharge,
                            FragmentIon, ProductCharge, IsotopeLabelType, TotalArea,
                            Times, Intensities, MinStartTime, MaxEndTime, File)]

  if (nrow(chrom.data) == 0) {
    errorReporting("Input data is not appropriately formatted for TargetedMSQC. Please check the RemoveIsotopePair and RemovePeakBoundary columns in output$removed.")
  } else {
    # substitute NA values in the TotalArea column with 0.
    if (nrow(chrom.data[is.na(TotalArea)]) > 0) {
      warnings("NA values are detected in the peak area column of the input and replaced with 0")
      chrom.data[is.na(TotalArea), TotalArea := 0]
    }

    # build the chromatogram group for each (File,FileName,peptide,precursorcharge) trio
    chrom.data[, ChromGroup := buildPeakGroup( c(.BY, .SD)),
               .(File, FileName, PeptideModifiedSequence, PrecursorCharge)]

    # build the peak group for each (File,FileName,peptide,precursorcharge) trio:
    #the only difference with the chromatogram group is that the peak boundaries are applied.
    packageStartupMessage("Applying Peak Boundaries ...", appendLF = F)
    chrom.data[, PeakGroup := ApplyPeakBoundary(chromGroup = ChromGroup[[1]],
                                                minStartTime = MinStartTime[1],
                                                maxEndTime = MaxEndTime[1]),
               .(File, FileName, PeptideModifiedSequence, PrecursorCharge)]
    packageStartupMessage(" Done")

    # remove those rows where the peak is not a proper peak object. This can
    # happen if the peak boundaries are too narrow and the peak has only 3 or
    # fewer points across it.
    removed.na.peaks <- chrom.data[is.na(PeakGroup),
                                   .(File, FileName, PeptideModifiedSequence,
                                     PrecursorCharge, FragmentIon, ProductCharge,
                                     IsotopeLabelType, MinStartTime, MaxEndTime)]

    if (nrow(removed.na.peaks) > 0) {
      removed.na.peaks[,':='(RemoveIsotopePair = FALSE,
                             RemovePeakBoundary = TRUE,
                             Remove = TRUE)]
      removed <- rbind(removed,removed.na.peaks)
    }

    chrom.data <- chrom.data[!is.na(PeakGroup)]
    message("Adding Peak Area to data")
    # add the peak area column to the data frame
    chrom.data[, id := .I]
    chrom.data[, PeakArea := peakArea(c(.SD, .BY)), by = id]
    chrom.data[, id := NULL]
    message("Calculating sum of AUC's for individual transitions")
    # calculate the sum of AUCs of individual transitions and assign to SumArea
    chrom.data[, SumArea := sum(PeakArea),
               .(PeptideModifiedSequence, PrecursorCharge, IsotopeLabelType, FileName,
                 File)]

    # calculate sum of transitions with the FragmentIon name of "sum". Add the
    # new "sum" transition rows to the dataframe
    message("Calculating sum of transitions with the FragmentIon")
    tmp <- chrom.data[, .(PeakGroup = CalculateTransitionSum(peak = PeakGroup[[1]]),
                          MinStartTime = MinStartTime[1],
                          MaxEndTime = MaxEndTime[1]),
               .(File, FileName, PeptideModifiedSequence, PrecursorCharge)]
    tmp[, ':='(FragmentIon = 'sum', ProductCharge = 0L, Times = NA, Intensities = NA,
               ChromGroup = NA)]
    tmp <- rbind(
      copy(tmp)[,IsotopeLabelType := endogenous.label],
      copy(tmp)[,IsotopeLabelType := standard.label]
    )

    tmp[, SumArea := peakArea(PeakGroup, req = paste0('sum.0.', IsotopeLabelType)),
        by = .(File, FileName, PeptideModifiedSequence, PrecursorCharge, IsotopeLabelType)]
    tmp[, PeakArea := SumArea]

    chrom.data[, TotalArea := NULL]
    chrom.data <- rbind(chrom.data, tmp)
  }
  # return output
  return(list(data = chrom.data, removed = removed))
}


#' Calculate an ensemble of QC features/metrics for a dataset that has been
#' cleaned up by CleanUpChromatograms.
#'
#' The function takes the output of CleanUpChromatograms as input and for each
#' transition pair of each peptide, calculates a number of QC metrics that are
#' used as features for training a predictive peak QC model.
#'
#' @param data A dataframe that contains identifier columns and peak groups of
#' class peakObj for each transition peak. This dataframe is the output of
#' CleanUpChromatograms (output$data).
#' @param blanks A dataframe with two columns of File and FileName. The FileName
#' column contains the names of blank runs for each Skyline File in the File
#' column. The values in File and FileName should be consistent with corresponding
#' columns in input data. Example: If Skyline document "SkylineFile1" includes
#' three blank runs of "Blank1", "Blank2" and "Blank3",
#' blanks <- data.frame(File = "SkylineFile1", FileName = c("Blank1","Blank2","Blank3")).
#' The function uses the blank runs to estimate the peak intensity at LLOQ for
#' each pepide and transition. If blanks = NA, the function uses the intensity.threshold
#' by default as approximation of intensity at LLOQ for all peaks. This threshold
#' can be adjusted. (Not Fully implemented in V2)
#' @param intensity.threshold The threshold that the function uses as an
#' approximation of the intensity at LLOQ for all peptides. If blank samples are
#' not provided to estimate the intensity at LLOQ, this value is used as an
#' approximation to determine transition peaks that are below limit of quantitation.
#' @param endogenous.label Label of the endogenous analyte: default is "light"
#' @param standard.label Label of the spiked-in isotopically labeled analyte:
#' default is "heavy"
#' @param export.features Logical parameter to indicate whether calculated features
#' should be exported as a .csv file.
#' @param feature.path Path to the directory, where the features should be saved
#' if export.features = TRUE.
#'
#' @return A list with the following objects:
#'               features: A dataframe with columns that contain the calculated
#'               QC features for each transition pair in the input data.
#'
#' @import data.table
#'
#' @examples
#'
#' data.features <- ExtractFeatures(data = data.CSF$data,
#'                                 export.features = FALSE,
#'                                 intensity.threshold = 1000)
#'
ExtractFeatures <- function(data, blanks = NA, intensity.threshold = 1000,
                            endogenous.label = "light",
                            standard.label = "heavy", export.features = FALSE,
                            feature.path = "", ...) {

  if(!is.na(blanks)){
    stop('What is blank?')
  } else {
    data[, aboveThreshold := PeakArea > intensity.threshold]
  }
  # this data is already captured in peak group object
  data[,':='(Times = NULL, Intensities = NULL)]
  labs <- list('endogenous' = endogenous.label, 'standard' = standard.label)
  # create a id column
  # rowwise operations performed in TMSQC V1.0, idea continued with the newer
  # implementation untill more clarity achieved
  data[, id := .I]

  message(Sys.time(),' : Extracting Jaggedness Features')
  jag.cols <- c('PeakGroupJaggedness.m','TransitionJaggedness','IsotopeJaggedness')
  data[, (jag.cols) := CalculatePeakJaggedness(
    PeakGroup[[1]],
    FragmentIon=FragmentIon,
    ProductCharge=ProductCharge,
    IsotopeLabelType=as.character(IsotopeLabelType)), by = id]
  message(Sys.time(),' : Done')

  message(Sys.time(),' : Extracting Symmetry Features')
  sym.cols <- c('PeakGroupSymmetry.m', 'TransitionSymmetry', 'IsotopeSymmetry')
  data[, (sym.cols) := CalculatePeakSymmetry(
    PeakGroup[[1]],
    FragmentIon=FragmentIon,
    ProductCharge=ProductCharge,
    IsotopeLabelType=as.character(IsotopeLabelType)), by = id]
  message(Sys.time(),' : Done')

  message(Sys.time(),' : Extracting Similarity Features')
  sim.cols <- c('PeakGroupSimilarity.m','PairSimilarity','IsotopeSimilarity')
  data[,(sim.cols) := CalculatePeakShapeSimilarity(
    PeakGroup[[1]],
    FragmentIon=FragmentIon,
    ProductCharge=ProductCharge,
    IsotopeLabelType=IsotopeLabelType), by = id]
  message(Sys.time(), ' : Done')

  message(Sys.time(),' : Extracting Elution Shift Features')
  els.cols <- c('PeakGroupShift','PairShift','IsotopeShift','TransitionShift')
  data[,(els.cols) := CalculatePeakElutionShift(
    PeakGroup[[1]],FragmentIon=FragmentIon,
    ProductCharge=ProductCharge,
    IsotopeLabelType=IsotopeLabelType), by = id]
  message(Sys.time(), ' : Done')

  message(Sys.time(),' Extracting FWHM Feature')
  fwhm.cols <- c('PeakGroupFWHM.m','PeakGroupFWHM2base.m',
                 'TransitionFWHM2base','TransitionFWHM','IsotopeFWHM2base',
                 'IsotopeFWHM')
  data[, (fwhm.cols) := CalculateFWHM(peak = PeakGroup[[1]],
                                      FragmentIon=FragmentIon,
                                      ProductCharge=ProductCharge,
                                      IsotopeLabelType=as.character(IsotopeLabelType)),
       by = id]
  message(Sys.time(), ' : Done')

  message(Sys.time(), ' : Extracting Modality Features')
  mod.cols <- c('PeakGroupModality.m','TransitionModality','IsotopeModality')
  data[, (mod.cols) := CalculateModality(
    PeakGroup[[1]],
    FragmentIon=FragmentIon,
    ProductCharge=ProductCharge,
    IsotopeLabelType=as.character(IsotopeLabelType)), by = id]
  message(Sys.time(),' : Done')

  message(Sys.time(),' : Calculating Peak Max Intensity for each transition')
  data[, TransitionMaxIntensity := CalculatePeakMaxIntensity(
    PeakGroup[[1]],
    FragmentIon=FragmentIon,
    ProductCharge=ProductCharge,
    IsotopeLabelType=as.character(IsotopeLabelType)),
    by =id]
  message(Sys.time(),' : Done')

  message(Sys.time(),'  :Calculating Max Boundary Intensity for each transition')
  data[,TransitionMaxBoundaryIntensity := CalculateMaxBoundaryIntensity(
    PeakGroup[[1]],
    FragmentIon=FragmentIon,
    ProductCharge=ProductCharge,
    IsotopeLabelType=as.character(IsotopeLabelType)),
    by =id]
  message(Sys.time(), ' : Done')

  message(Sys.time(),' : Extracting other derived features')
  data[, TransitionMaxBoundaryIntensityNormalized := TransitionMaxBoundaryIntensity/data$TransitionMaxIntensity]
  data[!is.finite(TransitionMaxBoundaryIntensityNormalized), TransitionMaxBoundaryIntensityNormalized := 0]

  merge.cols <- c('File', 'FileName', 'PeptideModifiedSequence', 'PrecursorCharge',
                  'ProductCharge', 'FragmentIon')

  #calculate area and peakcenter features, remove unwanted columns
  data[, ':='(PeakCenter = (MaxEndTime+MinStartTime)/2,
              Area2SumRatio = PeakArea/SumArea, id = NULL, PeakGroup = NULL,
              ChromGroup = NULL, MinStartTime = NULL, MaxEndTime = NULL)]
  data[, Area2SumRatio := impute_nonfinite(Area2SumRatio, '0')]
  # calculate 'mean' features
  data <- merge(data,
                data[aboveThreshold == T,.(MeanArea2SumRatio = mean(Area2SumRatio, na.rm = T),
                                           MeanTransitionFWHM = mean(TransitionFWHM, na.rm = T),
                                           MeanPeakCenter = mean(PeakCenter, na.rm = T)),
                     .(PeptideModifiedSequence, PrecursorCharge, File, FragmentIon, ProductCharge, IsotopeLabelType)],
                by = c('PeptideModifiedSequence', 'PrecursorCharge', 'File',
                       'FragmentIon', 'ProductCharge', 'IsotopeLabelType'),
                all.x = T)

  data[,':='(
    MeanIsotopeRatioConsistency = abs(Area2SumRatio - MeanArea2SumRatio)/MeanArea2SumRatio,
    MeanIsotopeFWHMConsistency = abs(TransitionFWHM - MeanTransitionFWHM)/MeanTransitionFWHM,
    MeanIsotopeRTConsistency = abs(PeakCenter - MeanPeakCenter)/MeanPeakCenter
  )]
  # impute NA, NaN, Null, 0 values
  data[,':='(
    MeanIsotopeRatioConsistency = impute_nonfinite(MeanIsotopeRatioConsistency, 'max'),
    MeanIsotopeFWHMConsistency = impute_nonfinite(MeanIsotopeFWHMConsistency,'max'),
    MeanIsotopeRTConsistency = impute_nonfinite(MeanIsotopeRTConsistency, 'max')
  )]

  data[,Area2SumRatioCV := sd(Area2SumRatio, na.rm = TRUE)/mean(Area2SumRatio, na.rm = TRUE),
       .(File, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge,
         IsotopeLabelType)]

  tmp <- dcast(data = data[,c(merge.cols, 'IsotopeLabelType', 'Area2SumRatio'), with = F],
               formula = ...~IsotopeLabelType, value.var = 'Area2SumRatio')
  tmp[, PairRatioConsistency_endogenous := abs(get(endogenous.label)-get(standard.label))/get(standard.label)]
  tmp[, PairRatioConsistency_endogenous := impute_nonfinite(PairRatioConsistency_endogenous, 'max')]

  tmp <- merge(tmp,
               tmp[FragmentIon != 'sum', .(PeakGroupRatioCorr_endogenous = cor(
                 get(endogenous.label), get(standard.label),method = 'pearson')),
                 .(PeptideModifiedSequence, PrecursorCharge, FileName,File)],
               by = c('PeptideModifiedSequence', 'PrecursorCharge', 'FileName',
                      'File'), all.x = T)

  tmp[is.na(PeakGroupRatioCorr_endogenous), PeakGroupRatioCorr_endogenous := 0]
  tmp[,(c(endogenous.label, standard.label)) := NULL]

  tmp <- merge(tmp,
               dcast(data = data[,c(merge.cols, 'IsotopeLabelType',
                                    'TransitionFWHM'), with = F],
                     formula = ...~IsotopeLabelType,value.var = 'TransitionFWHM'),
               by = merge.cols, all.x = T)

  tmp[,PairFWHMConsistency_endogenous := abs(get(endogenous.label)-get(standard.label))/get(standard.label)]
  tmp[,(c(endogenous.label, standard.label)) := NULL]
  tmp[, PairFWHMConsistency_endogenous:= impute_nonfinite(PairFWHMConsistency_endogenous,'max')]

  tmp <- merge(tmp,
               dcast(data = data[,c(merge.cols,'IsotopeLabelType', 'PeakArea'), with = F],
                     formula = ...~IsotopeLabelType, value.var = 'PeakArea'),
               by = merge.cols, all.x = T)
  tmp[, Endogenous2StandardRatio := (get(endogenous.label))/(get(standard.label))]
  tmp[, c(standard.label, endogenous.label) := NULL]
  message(Sys.time(),' : Done')

  cast.cols <- names(data)[!names(data) %in% c('File', 'FileName',
                                               'PeptideModifiedSequence',
                                               'FragmentIon', 'PrecursorCharge',
                                               'ProductCharge')]

  #TODO check if formula  = ...~IsotopeLabelType works in the cast below
  features <- dcast(
    data = data,
    File+FileName+PeptideModifiedSequence+FragmentIon+PrecursorCharge+ProductCharge~IsotopeLabelType,
    value.var = cast.cols[-c(1,2,4)]
  )
  change_name <- names(features)
  change_name <- gsub(endogenous.label,'endogenous',change_name)
  change_name <- gsub(standard.label,'standard',change_name)
  names(features) <- change_name

  features <- merge(features, tmp, merge.cols)
  features <- features[, req.features(), with = F]

  if(export.features){
    message(Sys.time(),' : Exporting features')
    if(!dir.exists(feature.path)){
      dir.create(feature.path)
    }
    fwrite(features, file.path(feature.path,"features.csv"))
  }
  message(Sys.time(),' : Feature Extraction Complete')
  invisible(features)
}



#' Extract the peak from a chromatogram by restricting the chromatogram to peak
#' boundaries.
#'
#' The function takes a peak group object and the peak boundary as
#' input. The input peak group is usually a chromatogram with a time range that
#' is wider than the peak of interest. This function extracts peak of interest
#' from the chromatogram by limiting the time and signal intensities of the peak
#' to a range that is within the provided peak boundaries.
#'
#' @param peak A peak group object
#' @param boundary a numeric vector of size two with min start time and max end
#' time of the peak (peak boundaries)
#'
#' @return A peak group object with time and intensity vectors that
#' are within the boundaries of the peak
#'
#' @export
#' @import pracma
#' @examples
#'
#' chrom <- data.CSF$data$ChromGroup[[1]]
#' PlotChromPeak(chrom)
#' peak <- ApplyPeakBoundary(chrom,c(data.CSF$data$MinStartTime[[1]],data.CSF$data$MaxEndTime[[1]]))
#' PlotChromPeak(peak)

ApplyPeakBoundary <- function(chromGroup, minStartTime, maxEndTime){
  tryCatch({
    idx <-  which(chromGroup$time > minStartTime & chromGroup$time < maxEndTime)
    if (length(idx) == 3) {
      if (idx[1] > 1){
        idx <- c(idx[1] - 1,idx)
      } else{
        idx <- c(idx,tail(idx,1) + 1)
      }
    }

    time <-  chromGroup$time[idx]
    peak.sig <- chromGroup$sig[idx,]
    area <- sapply(peak.sig, function(x) pracma::trapz(time,x))

    obj <- list(time = time, sig = peak.sig, area = area)
    if(check_peak_group(t = time, s = peak.sig, a = area)){
      NA
    } else{
      list(list(obj))
    }
  }, error = function(e){
    traceback()
    errorReporting(e$message)
  })
}
