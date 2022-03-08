## This file contains the functions used to build and evaluate the training model.


#' Create a training set template that can be shared with an analyst for manual
#' annotation.
#'
#' The function takes the directory that contains chromatogram of Skyline
#' documents and creates a template for manual annotation of the peaks. The
#' template is saved as a .csv file in the directory specified by template.path.
#' This templae should be shared with an analyst for manual annotation and the
#' "Status" column should be filled with "flag" or "ok" for low and high quality
#' peaks, respectively. Any comments about the manual annotation can be entered
#' in the "Notes" column.
#'
#' @param chromatogram.path Path to the directory containing the .tsv files of
#' the peak chromatograms. For each Skyline document, this file is exported from
#' Skyline through File > Export > Chromatograms. Here, check runs of interest
#' and include Precursors, Products, Base Peaks and TICs. Each chromatogram .tsv
#' file corresponds to a single Skyline document, which may contain any number
#' of runs. Multiple chromatogram files, corresponding to multiple Skyline
#' documents can be copied into the chromatogram.path directory. For each
#' chromatogram file in this folder, there should be a peak boundary file with
#' an identical name in peak.boundary.path directory.
#' @param template.path Path to the directory, where the template file will be
#' saved.
#' @param training.filename.list List of the runs that are going to be used for
#' training. if set to defaul of "all", all the runs in the chromatogram file
#' will be used for training.
#' @param endogenous.label Label of the endogenous analyte: default is "light"
#' @param standard.label Label of the spiked-in isotopically labeled analyte:
#' default is "heavy"
#' @param iRT.list List of iRT standards used in the experiment. These peptides
#' will be removed from the training set.
#'
#' @return Saves the template .csv file in template.path
#'
#' @export
#'
#' @examples
#'
#' extdata.path <- system.file("extdata",package = "TargetedMSQC")
#' project.folder.name <- "CSF_Panel"
#' project.path <- file.path(extdata.path,project.folder.name)
#' chromatogram.path <- file.path(project.path,"Chromatograms")
#' template.path <- file.path(project.path,"Templates")
#' MakeTemplate(chromatogram.path = chromatogram.path,
#'               template.path = template.path,
#'               endogenous.label = "light",standard.label = "heavy")

MakeTemplate <- function(chromatogram.path, template.path = NULL,
                         training.filename.list = "all",
                         endogenous.label = "light",
                         standard.label = "heavy" ,
                         iRT.list = iRTList(), ...) {
  # read all the tsv files in the chromatogram directory
  chromatogram.files <- sort(list.files(chromatogram.path,pattern = "\\.tsv", full.names = T))

  lapply(chromatogram.files, function(file){
    dt <- readData(file, sep = '\t')
    # remove the rows with an isotope label that is not included in the pair
    dt[, IsotopeLabelType := factor(tolower(IsotopeLabelType),
                                    levels = c(endogenous.label, standard.label))]
    # remove the rows that correspond to iRTs
    dt <- dt[!is.na(IsotopeLabelType) | !PeptideModifiedSequence %in% iRT.list]
    # assign the skyline file name as a column, empty status and notes columns
    dt[,':='(File = file, Status = NA_character_, Notes = NA_character_)]
    # QC.data holds File,FileName, Peptide, Transition and Isotopelabel Info
    dt <- dt[,.(File, FileName, PeptideModifiedSequence, PrecursorCharge, FragmentIon,
                ProductCharge, Status, Notes)]
    dt <- unique(dt)
    # filter to keep the rows in the training.filename.list
    if(!identical(training.filename.list, 'all')){
      dt <- dt[, (training.filename.list) := NULL]
    }
    # if the template path does not exist, create it
    if(!is.null(template.path)){
      if(!dir.exists(template.path)){
        dir.create(template.path)
      }
    }else{
      errorReporting("Please provide a path to save the template, use Arg: 'template.path'")
    }

    template.file <- file.path(template.path,
                               paste(gsub("\\..*", "", basename(file)),
                                     'training_template.csv', sep = '_'))
    fwrite(dt, file = template.file)
    message(Sys.time(), 'Template file created : \n',template.file)
  })
}


#' Create the dataset used as input to the TrainQCModel function.
#'
#' The function takes the feature dataframe (containing the calculated ensemble
#' of QC features) and the annotated training dataframe and merges them to create
#' the input of the TrainQCModel function.
#'
#' @param feature.data A dataframe that contains peak identifiers (File,FileName,
#' PeptideModifiedSequence,FragmentIon,IsotopeLabelType,PrecursorCharge and ProductCharge)
#' as well as  QC metrics calcualted for each transition pair. This dataframe is
#' the output of ExtractFeatures function (output$features).
#' @param feature.path Alternative to providing the feature.data, a path to the
#' directory that contains the features.csv file exported by ExtractFeatures
#' function can be provided. This file can be generated by by setting
#' export.features = TRUE and specifying the feature.path in ExtractFeatures
#' function. If both feature.data and feature.path are provided, feature.path is
#' ignored.
#' @param training.path Path to the directory containing the annotated template
#' file. The template file is generated using the MakeTemplate function should be
#' manually annotated by an expert analyst and saved as a separate .csv file. The
#' path to this annotated file should be provided. Please note that the directory
#' should contain only annotated .csv files that are meant to be in the study.
#' @param training.data Alternative to providing a path, the training data can
#' be provided explicitly.
#'
#'
#' @return A list with the following objects:
#'               data.merged: A dataframe that is product of merging and cleaning
#'               up feature.data and training.data.
#'               feature.data: The input feature.data
#'               training.data: The annotated input training.data
#'               data.training.feature: A dataframe that is product of merging
#'               and cleaning up feature.data and training.data. This dataframe
#'               contain only peaks that have been manually annotated and are
#'               included in training.data This data can be used by the
#'               TrainQCModel for training a predictive peak QC model.
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
#' training.path <- file.path(project.path,"Training")
#' data.set <- MakeDataSet(feature.data = data.features.CSF$features,training.path = training.path)
#' feature.path <- file.path(project.path,"Features")
#' data.set <- MakeDataSet(feature.path = feature.path,training.path = training.path)

MakeDataSet <- function(feature.path = NULL, training.path = NULL,
                        feature.data = NULL, training.data = NULL){

  non.numeric.cols <- c("IsotopeLabelType","File","FileName",
                        "PeptideModifiedSequence", "FragmentIon","PrecursorCharge",
                        "ProductCharge")

  if(is.null(feature.data) && !is.null(feature.path)){
    feature.files <- sort(list.files(feature.path, pattern = '.csv'))
    feature.data <- rbindlist(lapply(feature.files, fread))

    col.class <- unlist(feature.data[, lapply(.SD, class)])
    numeric.cols <- names(col.class)[!names(col.class) %in% non.numeric.cols]
    change.class <- names(which(col.class[numeric.cols] != 'numeric'))

    if(length(change.class) != 0){
      feature.data[, (change.class) := lapply(.SD, as.numeric), .SDcols = change.class]
    }
  }

  if(!is.null(training.path) && is.null(training.data)){
    training.files <- sort(list.files(training.path, pattern = '.csv'))
    training.data <- rbindlist(lapply(training.files, fread))
  }

  feature.data <- feature.data[FragmentIon != 'sum']
  feature.data[, ':='(
    FragmentIon = factor(FragmentIon),
    File = factor(File, levels = sort(unique(File)))
  )]

  training.data[,':='(Status = tolower(Status), Notes = tolower(Notes))]
  training.data <- training.data[Status != '#N/A' | !is.na(Status)]
  training.data[Notes == '#N/A' | is.na(Notes), Notes := 'ok']

  train.data.tmp <- training.data[PeptideModifiedSequence %in% feature.data$PeptideModifiedSequence]
  # Factors by default are avoided in data.table this logic can be revisited later on
  # train.data.tmp[, ':='(
  #   PeptideModifiedSequence = droplevels(PeptideModifiedSequence)
  # )]

  data.merged <- merge(feature.data, train.data.tmp, by = non.numeric.cols[-1],
                       all.x = T)

  list('data.merged' = data.merged, 'feature.data' = feature.data,
       'training.data' = training.data, 'data.training.feature' = data.merged[!is.na(Status)])
}


#' Train a binary classification model to flag peaks with poor chromatography or
#' interference.
#'
#' The function acts as a wrapper to several functions from the caret package to
#' train and optimize a binary predictive peak QC model for the provided training
#' data. Twenty percent of the training dataset is randomly selected as validation
#' set and left out from the training process to estimate the performance of the
#' models on unseen data. The features are mean centered and scaled by diving by
#' the standard deviation before being used for training. Repeated 10-fold cross
#' validation (3 repeats) is applied to the remainder of the training set to
#' minimize over-fitting. The model offering the highest accuracy is used and
#' returned by the function.
#'
#' @param data.merged A dataframe that contains peak identifiers (File,FileName,
#' PeptideModifiedSequence,FragmentIon,IsotopeLabelType,PrecursorCharge and ProductCharge),
#' the calculated QC metrics as well as the Status assigned by the expert analyst
#' to each transition pair. data.merged is usually output of MakeDataSet function
#' (output$data.merged or output$data.training.feature).
#' @param response.var This variable indicates the name of the column that stores
#' the "ok" and "flag" labels for the transition pairs in the training data.
#' @param description.columns If the input dataframe contains columns corresponding
#' to description variables (such as Notes), it should be indicated here.
#' Description and identifier columns will be removed from the data before training
#' the model.
#' @param method The machine learning algorithm for training the classifier. The
#' algorithm can be chosen from the list of available packages in caret
#' \url{https://topepo.github.io/caret/available-models.html}. The following have
#' been tested: RRF, regLogistic, svmLinear3, svmPoly, kknn. Before using TrainQCModel
#' with any of these packages, you will need to first install the machine learning
#' package using the install.packages command.
#' @param metric The performance metric used for selecting the model. The metric
#' can be Accuracy or ROC (AUC).
#' @param tuneGrid Use this parameter of you want to specify  tuneGrid for the
#' caret train method. Otherwise, set tuneGrid to NULL. See the caret package help
#' for more details: \url{https://topepo.github.io/caret/model-training-and-tuning.html}.
#' @param random.seed To fix the random seed for splitting the dataset into training
#' and validation and the data splitting for cross validation, provide a vector
#' of length 2 e.g. random.seed = c(1000,2000). This is particularly useful if
#' you want to compare multiple models with the same data split.
#' @param export.model A Logical parameter to indicate whether the model should
#' be saved. If export.model = TRUE the model will be saved in model.path.
#' @param model.path Path to the directory where the model will be saved if export.model = TRUE.
#'
#' @return A list with the following objects:
#'               model: Trained model to flag peaks with poor chromatography or
#'               interference.
#'               performance.testing: Confusion matrix of applying the model on
#'               the unseen validation data (20% of the input data). This parameter
#'               can help evaluate the performance of the model on unseen data
#'               and identify potential overfitting issues.
#'               model.file.path: If export.model = TRUE and the model is saved,
#'               the path and file name for the model is stored in this field.
#'
#' @export
#'
#' @import caret
#' @import RRF
#' @import doParallel
#' @examples
#'
#' rrf.grid <-  expand.grid(mtry = c(2,10),
#'                          coefReg = c(0.5,1),
#'                          coefImp = c(0))
#'
#' model.rrf <- TrainQCModel(data.set.CSF$data.training.feature,
#'                           response.var = c("Status"),
#'                           description.columns = c("Notes"),
#'                           method = "RRF",
#'                           metric = "Accuracy",
#'                           tuneGrid = rrf.grid,
#'                           random.seed = c(100,200))
#'

TrainQCModel <- function(data.merged, response.var = c("Status"),
                         description.columns = c("Notes"), method = "RRF",
                         metric = c("Accuracy","ROC"),tuneGrid = NULL,
                         random.seed = NULL, export.model = FALSE,
                         model.path = NA, ...) {

  identifier.columns = c("File","FileName","PeptideModifiedSequence","FragmentIon",
                         "PrecursorCharge","ProductCharge")

  # The identifier/description/response columns are removed.
  rm.cols <- c(identifier.columns, description.columns, response.var)
  feature.only <-  data.merged[, !rm.cols, with = F]

  # response vector
  resp_vector <- data.merged[,get(response.var)]

  # a 10-fold repeated cross validation (3 repeats) is used for model optimization.
  # if metric = ROC, the trainControl
  if (metric == "ROC") {
    train_control <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 3,
                                  verboseIter = FALSE,sampling = "up",grid,
                                  classProbs = TRUE,summaryFunction = twoClassSummary)
  } else {
    train_control <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 3,
                                  verboseIter = FALSE,sampling = "up",grid)
  }

  # split the data into training and testing sets. Training set is used to optimize
  # model parameters and train the model. Testing set is used as unseen data to
  # estimate model performance and evaluate overfitting
  if (!is.null(random.seed)) set.seed(random.seed[1])
  trainIndex <- caret::createDataPartition(resp_vector,p = .8, list = F, times = 1)

#  datasetTrain <- data.merged.feature.only.transformed[trainIndex,]
  datasetTrain <- feature.only[trainIndex,]

  # if a seed is provided for controlling the randomness of the algorithm, apply here
  if (!is.null(random.seed)) set.seed(random.seed[2])

  #parallel system check
  dots <- list(...)
  if(dots$parallel & !is.na(dots$workers)){
    message(Sys.time(),' : Setting up ',dots$workers,' workers for parallel model training')
    cl <- parallel::makePSOCKcluster(dots$workers)
    doParallel::registerDoParallel(cl)
  }

  # Train the model using the training dataset
  if (!is.null(tuneGrid)) {
    model <- caret::train(feature.only[trainIndex,], resp_vector[trainIndex],
                   method = method, preProcess = c("center","scale"),
                   trControl = train_control, importance = TRUE, metric = metric,
                   tuneGrid = tuneGrid)

  } else {
      model <- caret::train(feature.only[trainIndex,], resp_vector[trainIndex],
                     method = method, preProcess = c("center","scale"),
                     trControl = train_control, importance = TRUE, metric = metric)

  }
  message(Sys.time(), ' : Model Training Complete')
  if(dots$parallel){
    message(Sys.time(),' : Stopping Cluster')
    parallel::stopCluster(cl)
  }

  # predicting the response on inputs
  response.prediction <- predict(model, newdata = feature.only)

  # model performance based on confusion matrix for unseen data (testing set)
  performance.testing <- confusionMatrix(
    as.factor(response.prediction[-trainIndex]),
    as.factor(resp_vector[-trainIndex]))

  QC.model = list(model = model, performance.testing = performance.testing,
                  model.file.path = model.path)

  # if export.features is true save the features in the csv file specified by feature.path
  if (export.model == TRUE) {
    # template file name
    tmp_file <- paste0("model_",method,format(Sys.time(), "_%Y%m%d_%H%M"),".rda")
    model.file <- file.path(model.path, tmp_file)
    # if the template path does not exist, create it
    if (!dir.exists(model.path)) {
      dir.create(model.path)
    }
    QC.model$model.file.path <- model.file
    save(QC.model, file = model.file)
  }
  QC.model
}

#' Apply a binary predictive peak QC model to flag peaks that suffer from poor
#' chromatography or interference.
#'
#' The function takes the feature dataframe (output of ExtractFeatures output$features)
#' and the trained predictive QC model (output of TrainQCModel) and applies the
#' model to the input feature data using the function provided in the caret package.
#'
#' @param data.feature A dataframe that contains peak identifiers (File,FileName,
#' PeptideModifiedSequence,FragmentIon,IsotopeLabelType,PrecursorCharge and ProductCharge)
#' as well as  QC metrics calcualted for each transition pair. This dataframe is
#' the output of ExtractFeatures function (output$features).
#' @param model  The predictive model of peak QC. The model is the output of TrainQCModel.
#' @param response.var If the input dataframe contains columns corresponding to
#' response variables, it should be indicated here. it should be indicated here.
#' Response and description columns as well as identifier columns will be removed
#' from the data before applying the model.
#' @param description.columns If the input dataframe contains columns corresponding
#' to description variables (Such as Notes), it should be indicated here. Response
#' and description columns as well identifier columns will be removed from the
#' data before applying the model.
#' @param flag.prob.threshold A numeric value between 0 and 1 which determines
#' the cut-off threshold for assigning classes to each peak based on corresponding
#' class probabilities. By default, the caret package uses a probability threshold
#' of 0.5. This parameter can be used to override the default probability threshold.
#' @param standard.intensity.threshold This parameter can be used to set an
#' intensity threshold to identify and flag transitions where the spiked-in standard
#' is too low. If the numerical value of desired intensity threshold is provided,
#' it is used to flag any transition whose standard signal intensity is below this
#' threshold. For such transitions, this will override the model output.
#' @param type If type = "prob", the function will return class probabilities for
#' the binary classification. This feature can be used only if the model supports
#' classification probabilities e.g. logistic regression and random forest.
#'
#'
#' @return A dataframe of the predicted response (final class and/or class probabilities)
#' appended to the input data.feature.
#'
#' @export
#'
#' @import caret
#'
#' @examples
#'
#' response.data <- ApplyQCModel(data.set.CSF$feature.data,
#'                               model.rrf.CSF,
#'                               response.var = c("Status"),
#'                               description.columns = c("Notes"),
#'                               flag.prob.threshold = 0.5,
#'                               type = "prob")


ApplyQCModel <- function(data.feature, model, response.var = c("Status"),
                         description.columns = c("Notes"), flag.prob.threshold = 0.5,
                         standard.intensity.threshold = NULL, type = NULL, ...) {

  if (!is.null(standard.intensity.threshold) & !is.numeric(standard.intensity.threshold)){
    errorReporting("standard.intensity.threshold must be numeric")
  }

  # default identifier.columns, which will be removed from data.merged
  identifier.columns = c("File","FileName","PeptideModifiedSequence","FragmentIon",
                         "IsotopeLabelType","PrecursorCharge","ProductCharge")

  rm.cols <- c(identifier.columns, description.columns, response.var)
  feature.only <-  data.feature[, !rm.cols, with = F]
  # predicting the response on inputs
  response.prob <- as.data.table(predict(model$model, newdata = feature.only,
                                                    type = "prob"))

  response.prob[, Status.prediction := ifelse(flag <= flag.prob.threshold,
                                                         'ok', 'flag')]
  response.prob[, Status.prediction := factor(Status.prediction,
                                              levels = c('flag', 'ok'))]

  # merge QC predictions with original data
  data.feature.pred <- cbind(data.feature,response.prob)

  # if (!is.null(type) && type == "prob") {
  #   response.prediction.prob <- data.frame(flag.prob.prediction = response.prediction.prob[,1],
  #                                          ok.prob.prediction = response.prediction.prob[,2])
  #   # merge QC predictions with original data
  #   data.feature.pred <- cbind(data.feature.pred,response.prediction.prob)
  # }

  # if a threshold for standard signal intensity is provided, flag any transition
  # whose max intensity is below the threshold:
  if (!is.null(standard.intensity.threshold)) {
    data.feature.pred[TransitionMaxIntensity_standard < standard.intensity.threshold,
                      Status.prediction := 'flag']
  }
  data.feature.pred
}
