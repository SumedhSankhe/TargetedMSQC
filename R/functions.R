iRTList <- function(){
  invisible(c("LGGNETQVR","AGGSSEPVTGLADK","VEATFGVDESANK","YILAGVESNK",
              "TPVISGGPYYER","TPVITGAPYYER","GDLDAASYYAPVR","DAVTPADFSEWS0K",
              "TGFIIDPGGVIR","GTFIIDPAAIVR","FLLQFGAQGSPLFK","LGGNEQVTR",
              "GAGSSEPVTGLDAK","VEATFGVDESNAK","YILAGVENSK","TPVISGGPYEYR",
              "TPVITGAPYEYR","DGLDAASYYAPVR","ADVTPADFSEWSK","GTFIIDPGGVIR",
              "GTFIIDPAAVIR","LFLQFGAQGSPFLK"))
}

req.features <- function(){
  invisible(
    c('File','FileName','PeptideModifiedSequence','PrecursorCharge',
      'FragmentIon','ProductCharge','Area2SumRatioCV_endogenous',
      'Area2SumRatioCV_standard','PeakGroupRatioCorr_endogenous',
      'PairFWHMConsistency_endogenous','PairRatioConsistency_endogenous',
      'PeakGroupJaggedness.m_endogenous','PeakGroupSymmetry.m_endogenous',
      'PeakGroupSimilarity.m_endogenous','PeakGroupShift_endogenous',
      'PeakGroupFWHM.m_endogenous','PeakGroupFWHM2base.m_endogenous',
      'PeakGroupModality.m_endogenous','TransitionJaggedness_endogenous',
      'TransitionJaggedness_standard','TransitionSymmetry_endogenous',
      'TransitionSymmetry_standard','PairSimilarity_endogenous','PairShift_endogenous',
      'TransitionFWHM2base_endogenous','TransitionFWHM2base_standard',
      'TransitionFWHM_endogenous','TransitionFWHM_standard','TransitionModality_endogenous',
      'TransitionModality_standard','IsotopeJaggedness_endogenous',
      'IsotopeJaggedness_standard','IsotopeSymmetry_endogenous',
      'IsotopeSymmetry_standard','IsotopeSimilarity_endogenous',
      'IsotopeSimilarity_standard','IsotopeShift_endogenous',
      'IsotopeShift_standard','IsotopeFWHM2base_endogenous','IsotopeFWHM2base_standard',
      'IsotopeFWHM_endogenous','IsotopeFWHM_standard','IsotopeModality_endogenous',
      'IsotopeModality_standard','MeanIsotopeRatioConsistency_endogenous',
      'MeanIsotopeRatioConsistency_standard','MeanIsotopeFWHMConsistency_endogenous',
      'MeanIsotopeFWHMConsistency_standard','MeanIsotopeRTConsistency_endogenous',
      'MeanIsotopeRTConsistency_standard','TransitionMaxIntensity_endogenous',
      'TransitionMaxIntensity_standard','TransitionMaxBoundaryIntensity_endogenous',
      'TransitionMaxBoundaryIntensity_standard','TransitionMaxBoundaryIntensityNormalized_endogenous',
      'TransitionMaxBoundaryIntensityNormalized_standard','TransitionShift_endogenous',
      'TransitionShift_standard'))
}

errorReporting <- function(e){
  f <- deparse(sys.calls()[[sys.nframe()-1]])
  stop(simpleError(e,call = f))
}

checkFileNames <- function(fnames){
  f <- gsub(".csv|.tsv","",fnames)
  f <- basename(f)
  if(length(f)/2 != length(unique(f))){
    errorReporting("For each chromatograph file in chromatogram.path, there should be a peak boundary file in peak.boundary.path with an identical file name and vice versa")
  }else{
    message("CheckFilenames: Success!")
  }
}

readData <- function(files, sep){
  message("Reading Data")
  data.table::rbindlist(lapply(files, function(f){
    d <- data.table::fread(file = f, sep = sep)
    d$File <- basename(f)
    names(d) <- gsub(" ","", names(d))
    d
  }))
}

loadData <- function(chrom.Path, peak.Path){

  chrom.Files <- list.files(path = chrom.Path, pattern = ".tsv", full.names = T)
  peak.Files <- list.files(path = peak.Path, pattern = ".csv", full.names = T)

  checkFileNames(list("chrom.Files" = chrom.Files, "peak.Files" = peak.Files))

  list("chrom" = readData(files = chrom.Files, sep = "\t"),
       "peak.boundary" = readData(files = peak.Files, sep = ","))
}

formatChromData <- function(dt,lvls){
  packageStartupMessage("Formatting Chromatograms ...", appendLF = F)
  dt[, IsotopeLabelType := tolower(IsotopeLabelType)]
  dt[, IsotopeLabelType := factor(IsotopeLabelType, levels = lvls)]
  dt[, grp := .GRP, .(File, FileName, PeptideModifiedSequence, PrecursorCharge,
                      FragmentIon, ProductCharge)]

  dc <- data.table::dcast(dt, grp~IsotopeLabelType, value.var = 'TotalArea')
  dc[is.na(get(lvls[1])) | is.na(get(lvls[2])) ,#| is.na(IsotopLabelType),
     RemoveIsotopePair := T]
  dc[is.na(RemoveIsotopePair), RemoveIsotopePair := F]
  dt <- dt[dc, on = 'grp']
  dt[, ':='(grp = NULL, File = gsub(".tsv","",File))]
  packageStartupMessage(" Done")
  dt
}

formatPeakData <- function(dt){
  # mark the peaks with missing peak boundaries or the ones with multiple peak
  # boundaries for removal. Note that sometimes skyline exports multiple lines
  # for a peak, where one line is NA. In these cases, I keep the non-NA line
  # and remove the NA line as redundant.
  packageStartupMessage("Formatting Peak Boundaries ...", appendLF = F)
  dt <- unique(dt)
  dt[, n := .N, .(FileName, File, PeptideModifiedSequence)]
  dt[(is.na(MinStartTime) | is.na(MaxEndTime)) & n > 1, Redundant := TRUE]
  dt <- dt[is.na(Redundant)]
  dt[, n := .N, .(FileName, File, PeptideModifiedSequence)]
  dt[is.na(MinStartTime) | is.na(MaxEndTime) | n > 1, RemovePeakBoundary := TRUE]
  dt[is.na(RemovePeakBoundary), RemovePeakBoundary := F]
  dt[, ':='(n = NULL, Redundant = NULL, File = gsub(".csv","",File))]
  packageStartupMessage(" Done")
  dt
}

combineChromPeak <- function(chrom, peak, lvls){
  cData <- formatChromData(dt = chrom, lvls = lvls)
  pData <- formatPeakData(dt = peak)

  packageStartupMessage("Combining Chromatograms and Peak Boundaries ...", appendLF = F)
  cData <- cData[pData, on = c('File','FileName','PeptideModifiedSequence','PrecursorCharge')]
  cData[, Remove := RemoveIsotopePair | RemovePeakBoundary]
  packageStartupMessage(" Done")
  cData
}

pad_zero <- function(x, t){
  d <- length(t) - length(x)
  if(d>0){
    x <- c(x, rep(0, length(t) - length(x)))
  } else if (d < 0){
    x <- x[1:length(t)]
  } else{
    'Welp'
  }
  x
}

check_peak_group <- function(t, s, a){
  length(t) * dim(s)[1] * dim(s)[2] * length(a) == 0
}

buildPeakGroup <- function(dt){
  message(Sys.time(),' : Building peak group for : ',
          paste(dt$File,dt$FileName,dt$PeptideModifiedSequence,dt$PrecursorCharge,
                collapse = "_"))

  time <-  lapply(strsplit(dt$Times, ','), as.numeric)
  time <- colMeans(do.call('rbind', time))

  itn <- lapply(strsplit(dt$Intensities,','), as.numeric)

  peak.sig <-lapply(seq_along(itn), function(x){
    pad_zero(itn[[x]], time)
  })

  area <- lapply(seq_along(itn), function(x) {
    pracma::trapz(time, peak.sig[[x]])
  })


  all.isotopes <- unique(as.character(dt$IsotopeLabelType))
  columns <- paste(dt$FragmentIon, dt$ProductCharge, dt$IsotopeLabelType,
                   sep = '.')

  single.peaks <- columns[which(
    mapply(function(x,columns){
      sum(grepl(substr(x,1,regexpr("\\.[^\\.]*$", x)),columns))
    }, columns, list(columns)) != 2
  )]

  if(length(single.peaks) > 0){
    message('Missing peaks identified')
    missing.isotopes <- mapply(function(x,all.isotopes){
      all.isotopes[which(all.isotopes != substr(x,regexpr("\\.[^\\.]*$", x) + 1,nchar(x)))]
    } ,single.peaks,list(all.isotopes))

    # names of the missing peaks
    missing.peaks <- sapply(single.peaks,function(x) substr(x,1,regexpr("\\.[^\\.]*$", x)))
    # colnames of the missing peaks that are to be imputated by zero
    imputed.peaks <- paste(missing.peaks,missing.isotopes,sep = "")

    message('Setting missing peak with 0 peak.sig and 0 area')
    pos <- length(dt$peak.sig)
    lapply(seq_along(imputed.peaks), function(x){
      area[[x+pos]] <- 0
      peak.sig[[x+pos]] <- rep(0, length(dt$peak.sig[[1]]))
    })
    columns <- c(columns, imputed.peaks)
  }

  peak.sig <- do.call('cbind', peak.sig)
  colnames(peak.sig) <- columns

  names(area) <- columns
  area <- do.call('c', area)

  obj <- list(
    list(
      list(time = time,#numeric
           sig =  as.data.frame(peak.sig),#data.frame
           area = area#numeric
      )))

  if(check_peak_group(t = time, s = peak.sig, a = area)){
    NA
  } else{
    obj
  }

}

applyPeakBoundary <- function(chromGroup, minStartTime, maxEndTime){
  pg <- lapply(seq_along(chromGroup), function(x){
    m <- minStartTime[x]
    M <- maxEndTime[x]
    idx <-  which(chromGroup[[x]]$time > m & chromGroup[[x]]$time < M)
    time <-  chromGroup[[x]]$time[idx]

    if (length(idx) == 3) {
      if (idx[1] > 1) idx <- c(idx[1] - 1,idx)
      else idx <- c(idx,tail(idx,1) + 1)
    }

    peak.sig <- chromGroup[[x]]$sig[idx,]
    area <- sapply(peak.sig, function(x) pracma::trapz(time,x))

    obj <- list(time = time, sig = peak.sig, area = area)
    if(check_peak_group(t = time, s = peak.sig, a = area)){
      NA
    } else{
      obj
    }
  })
  list(pg)

}

peakArea <- function(x, req=NA){
  req <- unique(req)
  if(is.na(req)){
    req <- paste(x$FragmentIon, x$ProductCharge, x$IsotopeLabelType, sep = '.')
    x$PeakGroup[[1]][['area']][[req]]
  }else{
    do.call('c',lapply(seq_along(x), function(y){
      x[[y]][['area']][[req]]
    }))
  }
}

plot_chrom <- function(peak, font.size = 14, transition.list = NA, label.list = NA,
                       split.label = TRUE){

  pdt <- data.table( peak$sig, 'time' = peak$time)
  pdt <- melt(pdt, id.vars = 'time')
  pdt[, transition := gsub('(.*)\\.\\w+', '\\1', variable)]
  pdt[, label := gsub(".*\\.","", variable)]

  if(!is.na(transition.list)){
    pdt <- pdt[transition %in% transition.list]
  }

  if(!is.na(label.list)){
    pdt <- pdt[label %in% label.list]
  }

  p <- ggplot(data = pdt,
              aes(x = time, y = value, color = transition, linetype = label))+
    geom_line()+
    theme_bw()+
    labs(x = "Retention Time (mins)", y = "Intensities", color = 'Transition',
         linetype = '')+
    scale_linetype_manual(values = c("solid", "twodash")) +
    theme(text = element_text(size = font.size))

  if(split.label){
    p <- p+facet_wrap(~toupper(label), scales = 'free')
  }
  p
}

ExtractFeatures <- function(data, blanks = NA, intensity.threshold = 1000,
                            endogenous.label = "light", standard.label = "heavy",
                            export.features = FALSE, feature.path = "", ...) {

  if(!is.na(blanks)){
    stop('What is blank?')
  } else {
    data[, aboveThreshold := PeakArea > intensity.threshold]
  }

  labs <- list('endogenous' = endogenous.label, 'standard' = standard.label)

  data[, id := .I]
  packageStartupMessage(Sys.time()," Extracting Jaggedness Features...", appendLF = F)
  jag.cols <- c('PeakGroupJaggedness.m','TransitionJaggedness','IsotopeJaggedness')
  data[, (jag.cols) := CalculatePeakJaggedness(PeakGroup[[1]],
                                                 FragmentIon=FragmentIon,
                                                 ProductCharge=ProductCharge,
                                                 IsotopeLabelType=as.character(IsotopeLabelType)), by = id]
  packageStartupMessage('Done')

  packageStartupMessage(Sys.time()," Extracting Symmetry Features...", appendLF = F)
  sym.cols <- c('PeakGroupSymmetry.m', 'TransitionSymmetry', 'IsotopeSymmetry')
  data[, (sym.cols) := CalculatePeakSymmetry(PeakGroup[[1]],
                                              FragmentIon=FragmentIon,
                                              ProductCharge=ProductCharge,
                                              IsotopeLabelType=as.character(IsotopeLabelType)), by = id]
  packageStartupMessage('Done')


  packageStartupMessage(Sys.time()," Extracting Similarity Features...", appendLF = F)
  sim.cols <- c('PeakGroupSimilarity.m','PairSimilarity','IsotopeSimilarity')
  data[,(sim.cols) := CalculatePeakShapeSimilarity(PeakGroup[[1]],
                                           FragmentIon=FragmentIon,
                                           ProductCharge=ProductCharge,
                                           IsotopeLabelType=IsotopeLabelType), by = id]
  packageStartupMessage('Done')

  packageStartupMessage(Sys.time()," Extracting Elution Shift Features...", appendLF = F)
  els.cols <- c('PeakGroupShift','PairShift','IsotopeShift','TransitionShift')
  data[,(els.cols) := CalculatePeakElutionShift(PeakGroup[[1]],FragmentIon=FragmentIon,
                                               ProductCharge=ProductCharge,
                                               IsotopeLabelType=IsotopeLabelType), by = id]
  packageStartupMessage('Done')

  packageStartupMessage(Sys.time()," Extracting FWHM Feature...", appendLF = F)
  fwhm.cols <- c('PeakGroupFWHM.m','PeakGroupFWHM2base.m',
               'TransitionFWHM2base','TransitionFWHM','IsotopeFWHM2base',
               'IsotopeFWHM')
  data[, (fwhm.cols) := CalculateFWHM(peak = PeakGroup[[1]],
                                    FragmentIon=FragmentIon,
                                    ProductCharge=ProductCharge,
                                    IsotopeLabelType=as.character(IsotopeLabelType)), by = id]
  packageStartupMessage('Done')

  packageStartupMessage(Sys.time()," Extracting Modality Features...", appendLF = F)
  mod.cols <- c('PeakGroupModality.m','TransitionModality','IsotopeModality')
  data[, (mod.cols) := CalculateModality(PeakGroup[[1]],
                                        FragmentIon=FragmentIon,
                                        ProductCharge=ProductCharge,
                                        IsotopeLabelType=as.character(IsotopeLabelType)), by = id]
  packageStartupMessage('Done')

  packageStartupMessage(Sys.time()," Calculating Peak Max Intensity for each transition...",
                        appendLF = F)
  data[, TransitionMaxIntensity := CalculatePeakMaxIntensity(PeakGroup[[1]],
                                                             FragmentIon=FragmentIon,
                                                             ProductCharge=ProductCharge,
                                                             IsotopeLabelType=as.character(IsotopeLabelType)),
       by =id]
  packageStartupMessage('Done')

  packageStartupMessage(Sys.time()," Calculating Max Boundary Intensity for each transition...",
                        appendLF = F)
  data[,TransitionMaxBoundaryIntensity := CalculateMaxBoundaryIntensity(PeakGroup[[1]],
                                                                    FragmentIon=FragmentIon,
                                                                    ProductCharge=ProductCharge,
                                                                    IsotopeLabelType=as.character(IsotopeLabelType)),
       by =id]
  packageStartupMessage('Done')

  packageStartupMessage(Sys.time()," Extracting other derived features...",
                        appendLF = F)
  data[, TransitionMaxBoundaryIntensityNormalized := TransitionMaxBoundaryIntensity/data$TransitionMaxIntensity]
  data[!is.finite(TransitionMaxBoundaryIntensityNormalized), TransitionMaxBoundaryIntensityNormalized := 0]

  merge.cols <- c('File', 'FileName', 'PeptideModifiedSequence', 'PrecursorCharge',
                'ProductCharge', 'FragmentIon')

  data[, ':='(PeakCenter = (MaxEndTime+MinStartTime)/2,
              Area2SumRatio = PeakArea/SumArea, id = NULL, PeakGroup = NULL,
              ChromGroup = NULL, Times = NULL, Intensities = NULL,
              MinStartTime = NULL, MaxEndTime = NULL)]
  data[!is.finite(Area2SumRatio), Area2SumRatio := 0]

  data[aboveThreshold == T,':='(MeanArea2SumRatio = mean(Area2SumRatio),
                                MeanTransitionFWHM = mean(TransitionFWHM),
                                MeanPeakCenter = mean(PeakCenter)),
       .(PeptideModifiedSequence, PrecursorCharge, File, FragmentIon, ProductCharge, IsotopeLabelType)]
  data[,':='(
    MeanIsotopeRatioConsistency = abs(Area2SumRatio - MeanArea2SumRatio)/MeanArea2SumRatio,
    MeanIsotopeFWHMConsistency = abs(TransitionFWHM - MeanTransitionFWHM)/MeanTransitionFWHM,
    MeanIsotopeRTConsistency = abs(PeakCenter - MeanPeakCenter)/MeanPeakCenter
  )]

  data$MeanIsotopeRatioConsistency[!is.finite(data$MeanIsotopeRatioConsistency)] <- max(data$MeanIsotopeRatioConsistency[is.finite(data$MeanIsotopeRatioConsistency)],na.rm = TRUE)
  data$MeanIsotopeFWHMConsistency[!is.finite(data$MeanIsotopeFWHMConsistency)] <- max(data$MeanIsotopeFWHMConsistency[is.finite(data$MeanIsotopeFWHMConsistency)],na.rm = TRUE)
  data$MeanIsotopeRTConsistency[!is.finite(data$MeanIsotopeRTConsistency)] <- max(data$MeanIsotopeRTConsistency[is.finite(data$MeanIsotopeRTConsistency)],na.rm = TRUE)

  data[,Area2SumRatioCV := sd(Area2SumRatio, na.rm = TRUE)/mean(Area2SumRatio, na.rm = TRUE),
       .(File, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge,
         IsotopeLabelType)]
  cast.cols <- names(data)[!names(data) %in% c('File', 'FileName',
                                                'PeptideModifiedSequence',
                                                'FragmentIon', 'PrecursorCharge',
                                                'ProductCharge')]
  features <- dcast(
    data = data,
    File+FileName+PeptideModifiedSequence+FragmentIon+PrecursorCharge+ProductCharge~IsotopeLabelType,
    value.var = cast.cols[-c(1,2,4)]
  )
  change_name <- names(features)
  change_name <- gsub(endogenous.label,'endogenous',change_name)
  change_name <- gsub(standard.label,'standard',change_name)
  names(features) <- change_name

  features <- features[, names(features)[names(features) %in% req.features()], with = F]


  tmp <- dcast(data = data[,c(merge.cols, 'IsotopeLabelType', 'Area2SumRatio'), with = F],
                    formula = ...~IsotopeLabelType, value.var = 'Area2SumRatio')
  tmp[, PairRatioConsistency_endogenous := abs(get(endogenous.label)-get(standard.label))/get(standard.label)]
  maxval <- max(tmp$PairRatioConsistency_endogenous[is.finite(tmp$PairRatioConsistency_endogenous)],
                na.rm = T)
  tmp[!is.finite(PairRatioConsistency_endogenous), PairRatioConsistency_endogenous := maxval]

  tmp[FragmentIon != 'sum', PeakGroupRatioCorr_endogenous := cor(get(endogenous.label), get(standard.label),
                                       method = 'pearson'),
           .(PeptideModifiedSequence, PrecursorCharge, FileName,File)]
  tmp[is.na(PeakGroupRatioCorr_endogenous), PeakGroupRatioCorr_endogenous := 0]
  tmp[,(c(endogenous.label, standard.label)) := NULL]

  tmp <- merge(tmp,
                    dcast(data = data[,c(merge.cols, 'IsotopeLabelType',
                                         'TransitionFWHM'), with = F],
                          formula = ...~IsotopeLabelType,value.var = 'TransitionFWHM'),
                    by = merge.cols)

  tmp[,PairFWHMConsistency_endogenous := abs(get(endogenous.label)-get(standard.label))/get(standard.label)]
  tmp[,(c(endogenous.label, standard.label)) := NULL]
  maxval <- max(tmp$PairFWHMConsistency_endogenous[is.finite(tmp$PairFWHMConsistency_endogenous)],
                na.rm = T)
  tmp[!is.finite(PairFWHMConsistency_endogenous), PairFWHMConsistency_endogenous := maxval]

  tmp <- merge(tmp,
               dcast(data = data[,c(merge.cols,'IsotopeLabelType', 'PeakArea'), with = F],
                     formula = ...~IsotopeLabelType, value.var = 'PeakArea'),
               by = merge.cols)
  tmp[, Endogenous2StandardRatio := (get(endogenous.label))/(get(standard.label))]
  tmp[, c(standard.label, endogenous.label) := NULL]
  packageStartupMessage('Done')

  features <- merge(features, tmp, merge.cols)
  features <- features[, req.features(), with = F]

  if(export.features){
    message(Sys.time()," Exporting features")
    if(!dir.exists(feature.path)){
      dir.create(feature.path)
    }
    fwrite(feature.path, "features.csv")
  }
  message(Sys.time()," Done")
  invisible(features)
}

plot_qc_summary <- function(data, runs = 'all', features = NULL, labels = NULL,
                            font.size = 14){

  if(runs=='all'){
    runs <- unique(data$FileName)
  }
  id.vars <-  c('File', 'FileName', 'PeptideModifiedSequence',
                'FragmentIon', 'PrecursorCharge',
                'ProductCharge')
  if(!is.null(features)){
    data <- data[, c(id.vars,features), with = F]
  }

  vplots <- lapply(runs, function(x){

    tmp <- melt(data = data[FileName == x], id.vars = id.vars)

    p <- ggplot(data = tmp, aes(x = variable, y = value))+
      geom_violin(fill = "#F2F2F2",scale = "width") +
      coord_flip() +
      labs(x = 'Features', y = 'Score', title = x)+
      theme_bw()+
      theme(text = element_text(size = font.size))

    if(!is.null(labels)){
      p <- p+
        geom_jitter(aes_string( color = labels),
                    position = position_jitter(0.25), size = 0.5)
    }
    p
  })

  tmp <- melt(data = data, id.vars = id.vars)
  p <- ggplot(data = tmp, aes(x = variable, y = value))+
    geom_violin(fill = "#F2F2F2",scale = "width") +
    coord_flip() +
    labs(x = 'Features', y = 'Score', title = 'All Runs')+
    theme_bw()+
    theme(text = element_text(size = font.size))

  if(!is.null(labels)){
    p <- p+
      geom_jitter(aes_string( color = labels),
                  position = position_jitter(0.25), size = 0.5)
  }
  vplots <- append(vplots,list(p))
  vplots
}
