require(methods)
require(fcluster)
require(class)
require(parallel)
require(doParallel)
require(Matrix)
require(readxl) ##Requires version on github, not CRAN
require(tools)
require(readr)

## Function definitions -------------------------------------------------------
pixelInfo <- setClass("pixelInfo", slots = c(name="character", className="character", patNum="numeric",
                                             scanNum="numeric", x="numeric", y="numeric", z="numeric", tic="numeric", ctic="numeric"))

makePixelInfo <- function(name, className, patNum, scanNum, x, y, z, tic, ctic) {
  pixelInfo <- setClass("pixelInfo", slots = c(name="character", className="character", patNum="numeric",
                                               scanNum="numeric", x="numeric", y="numeric", z="numeric", tic="numeric", ctic="numeric"))
  return(pixelInfo(name=name, className=className, patNum = patNum, scanNum=scanNum, x=x, y=y, z=z, tic=tic, ctic=ctic))
}

read_files <- function(fileList) {
  # Reads in CSV files into a list of tables
  #
  # Args:
  #   fileList: List of CSV filenames
  #
  # Returns:
  #   List of data frames

  # tableList <- vector("list", length(fileList))
  # for (i in 1:length(fileList)) {
  #   rawTable <- read.csv(file.path(fileList[i]), colClasses = "numeric", skip=4)
  #   rawTable <- Filter(function(x)!all(is.na(x)), rawTable)
  #   tableList[[i]] <- rawTable
  # }

  # tableList <- mclapply(fileList, function(x) read_file(x), mc.cores = 4)

  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  tableList <- foreach(i = 1:length(fileList), .export=c("read_file_data","findFirstNumeric"), .packages=c("doParallel","tools", "readxl", "readr")) %dopar% {
    read_file_data(fileList[[i]])
  }
  stopCluster(cl)
  registerDoSEQ()
  return(tableList)
}

read_headers <- function(fileList, className) {
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  tableList <- foreach(i = 1:length(fileList), .export=c("read_file_header", "findFirstNumeric", "makePixelInfo"), .packages=c("doParallel","tools", "readxl", "readr")) %dopar% {
    header <- read_file_header(fileList[[i]])
    numPixels <- ceiling(ncol(header)/2)
    baseName <- tail(strsplit(fileList[i], .Platform$file.sep)[[1]], n = 1)
    ptName <- file_path_sans_ext(baseName)

    foreach(k = 1:numPixels) %do% {
      scan <- NA_integer_
      x <- NA_integer_
      y <- NA_integer_
      z <- NA_integer_
      tic <- NA_integer_
      ctic <- NA_integer_


      for(j in 1:nrow(header)){
        if(header[j,k*2-1] == "Scan") {
          scan <- as.numeric(header[j,k*2])
        } else if(header[j,k*2-1] == "X") {
          x <- as.numeric(header[j,k*2])
        } else if(header[j,k*2-1] == "Y") {
          y <- as.numeric(header[j,k*2])
        } else if(header[j,k*2-1] == "Z") {
          z <- as.numeric(header[j,k*2])
        } else if(header[j,k*2-1] == "TIC") {
          tic <- as.numeric(header[j,k*2])
        }
      }
      makePixelInfo(name=ptName, className=className, scan=scan, x=x,
                    y=y, z=z, tic=tic, patNum=NA_integer_, ctic=ctic)
    }
  }
  stopCluster(cl)
  registerDoSEQ()
  return(tableList)
}

read_file_data <- function(filePath) {
  headerLines <- findFirstNumeric(filePath) - 2
  if(file_ext(filePath) == 'csv'){
    rawTable <- read_csv(filePath, col_types = cols(.default=col_double()), col_names=FALSE, skip=headerLines+1)
  }
  else if(file_ext(filePath) == 'xlsx' || file_ext(filePath) == 'xls'){
    rawTable <- read_excel(filePath, col_types = "numeric", skip=headerLines)
  }
  rawTable <- Filter(function(x)!all(is.na(x)), rawTable)
  return(rawTable)
}

read_file_header <- function(filePath) {
  headerLines <- findFirstNumeric(filePath) - 2
  if(file_ext(filePath) == 'csv'){
      if(!headerLines){
        rawTable <- read_csv(filePath, col_types = cols(.default=col_character()), col_names=FALSE)
      } else {
        rawTable <- read_csv(filePath, col_types = cols(.default=col_character()), col_names=FALSE, n_max = headerLines+1)
      }
    rawTable <- Filter(function(x)!all(x==""), rawTable)
  }
  else if(file_ext(filePath) == 'xlsx' || file_ext(filePath) == 'xls'){
      if(!headerLines){
        rawTable <- read_excel(filePath, col_types = "text", col_names=FALSE)
      } else {
        rawTable <- read_excel(filePath, col_types = "text", col_names=FALSE, n_max = headerLines)
      }
  }
  return(rawTable)
}

findFirstNumeric <- function(filePath) {
  if(file_ext(filePath) == 'csv'){
    rawTable <- read_csv(filePath, col_types = cols_only(X1=col_character()), col_names=FALSE)
  }
  else if(file_ext(filePath) == 'xlsx' || file_ext(filePath) == 'xls'){
    rawTable <- read_excel(filePath, col_types = "text", col_names=FALSE)
  }
  numericLocations <- suppressWarnings(sapply(rawTable[,1], function(x) !is.na(as.numeric(x))))
  return(which.max(numericLocations))
}

extract_peaks <- function(classData) {
  # Extract all peaks for data from class (i.e. cancer, normal)
  #
  # Args:
  #   classData: List of data frames
  #
  # Returns:
  #   Vector of m/z
  numFiles <- length(classData)
  fullMZ <- NULL

  # for (k in 1:numFiles) {
  #   fileData <- classData[[k]]
  #   numCols <- ncol(fileData)
  #   mz <- unlist(fileData[c(T,F)])
  #   mz <- mz[!is.na(mz)]
  #   fullMZ <- c(fullMZ, mz)
  # }

  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  fullMZ <- foreach(i = 1:numFiles, .combine=c) %dopar% {
    fileData <- classData[[i]]
    numCols <- ncol(fileData)
    mz <- unlist(fileData[c(T,F)])
    mz[!is.na(mz)]
  }
  stopCluster(cl)
  registerDoSEQ()
  return(fullMZ)
}

get_data_matrix_binning <- function(classData, filteredMZ, type = "sum", decimal = 2) {
  # Match most frequent peaks to individual pixels to form data matrices
  #
  # Args:
  #   classData: List of data frames
  #   filteredMZ: vector of most frequent peaks
  #
  # Returns:
  #   [1] matrix of intensities corresponding to most frequent peaks
  numColsAll <- sapply(classData, ncol)
  data <- matrix(NA, sum(numColsAll) / 2, length(filteredMZ))
  ii <- 0
  for (i in 1:length(classData)) {
    fileData <- as.data.frame(classData[[i]])
    for (k in seq(1, numColsAll[i], 2)) {
      pixelMZ <- fileData[, k]
      pixelMZ <- round(as.numeric(as.character(pixelMZ)), decimal)
      intens <- as.numeric(as.character(fileData[, k+1]))
      df <- data.frame(pixelMZ, intens)
      df <- aggregate(df['intens'], by = df['pixelMZ'], FUN = type)
      pixelMZ <- unname(unlist(df['pixelMZ']))
      intens <- unname(unlist(df['intens']))
      filterCounts <- match(pixelMZ, filteredMZ)
      intens <- intens[!is.na(filterCounts)]
      filterCounts <- filterCounts[!is.na(filterCounts)]
      intensOut <- rep(NA, length(filteredMZ))
      intensOut[filterCounts] <- intens
      ii <- ii + 1
      data[ii, ] <- intensOut
    }
  }
  data[is.na(data)] <- 0
  return(data)
}

get_data_matrix_binning_var <- function(classData, bins, filteredMZ, type = "sum") {
  # Match most frequent peaks to individual pixels to form data matrices
  #
  # Args:
  #   classData: List of data frames
  #   bins: vector of bins
  #
  # Returns:
  #   [1] matrix of intensities corresponding to bins
  numColsAll <- sapply(classData, ncol)
  data <- matrix(NA, sum(numColsAll) / 2, length(filteredMZ))
  ii <- 0
  for (i in 1:length(classData)) {
    fileData <- as.data.frame(classData[[i]])
    for (k in seq(1, numColsAll[i], 2)) {
      pixelMZ <- fileData[, k]
      pixelMZ <- bins[findInterval(round(as.numeric(as.character(pixelMZ)), 4), bins)]
      intens <- as.numeric(as.character(fileData[, k+1]))
      df <- data.frame(pixelMZ, intens)
      df <- aggregate(df['intens'], by = df['pixelMZ'], FUN = type)
      pixelMZ <- unname(unlist(df['pixelMZ']))
      intens <- unname(unlist(df['intens']))
      filterCounts <- match(pixelMZ, filteredMZ)
      intens <- intens[!is.na(filterCounts)]
      filterCounts <- filterCounts[!is.na(filterCounts)]
      intensOut <- rep(NA, length(filteredMZ))
      intensOut[filterCounts] <- intens
      ii <- ii + 1
      data[ii, ] <- intensOut
    }
  }
  data[is.na(data)] <- 0
  return(data)
}

get_data_matrix_clustering <- function(classData) {
  # Combine data matrixes for each file into single data matrix.
  #
  # Args:
  #   fileList: List of filenames
  #   classData: List of data frames
  #
  #
  # Returns:
  #   [1] matrix of intensities of m/z's
  #   [2] vector of patient labels
  # pt <- NULL
  data <- NULL

  numRowsAll <- sapply(classData,function(x) nrow(x$all))

  for (i in 1:length(classData)) {
    # baseName <- tail(strsplit(fileList[i], .Platform$file.sep)[[1]], n = 1)
    # ptName <- substring(baseName, 1, nchar(baseName) - 4)
    if(is.null(data)){
      data <- classData[[i]]$all
    }
    else{
      data <- rBind(data, classData[[i]]$all)
    }
    # pt <- c(pt, rep(ptName, numRowsAll[i]))
  }
  return(data)
}

get_cluster_matrix <- function(data, mz) {
  #Converts list containing data into list containing matrix corresponding to mz

  # matrix <- vector("list",length(data))
  # for(k in 1:length(data)){
  #   cat(k,fill=T)
  #   matrix[[k]]=formDataMatrix(data[[k]],mz)
  # }
  # return(matrix);

  # matrix <- mclapply(data, function(x, y) formDataMatrix(x, y), y=mz, mc.cores = 4)

  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  matrix <- foreach(i = 1:length(data), .export = "formDataMatrix", .packages=c("Matrix","class")) %dopar% {
    formDataMatrix(data[[i]], mz)
  }
  stopCluster(cl)
  registerDoSEQ()
  return(matrix)
}

formDataMatrix = function(d,
                          cen,
                          h = .05,
                          ncol = 2,
                          round.decimals = NULL) {
  # Function imported from Rob's ovarian sc.R. Do not call directly.
  # Use wrapper method get_cluster_matrix
  #
  # Want to avoid modifying this to ensure max compatibility with Rob.
  # May rewrite at later date. - jqlin 10/13/16
  #
  # @param d a data frame of data from one file
  # @param cen vector; m/z's produced by clustering
  # @param h height at which clustering was cut
  # @param ncol width in columns of each independent pizel's m/z and intensity list
  # @param round.decimals integer; the digits to round to
  # @return List of data matrix with intensities matched to cluster m/z's for each pizel
  #         and "pat" which numbers pixels in a file. "pat" is never used; unknown purpose.
  numscans = (ncol(d) + 1) / ncol

  all = Matrix(0, numscans, length(cen))
  pat = NULL
  num = 0
  ii = 0
  nc = ncol(d)
  for (k in seq(1, nc, ncol)) {
    num = num + 1
    x = d[, k]
    x = x[!is.na(x)]
    if (!is.null(round.decimals))
      x = round(x, round.decimals)
    set.seed(10)
    a = knn1(matrix(cen, ncol = 1), matrix(x, ncol = 1), 1:length(cen))
    # x=round(as.numeric(as.character(x)),2)
    xhat = cen[a]
    dd = abs(xhat - x) > h
    y = d[, k + 1]
    y = y[!is.na(y)]
    y[dd] = 0  # remove obs more than h away from a centroid
    yout = rep(0, length(cen))
    yout[a] = y
    ii = ii + 1
    all[ii, ] = yout
    pat = c(pat, num)
  }
  return(list(all = all, pat = pat))
}

assignPatNum <- function(pixelInfos){
  pixelInfos <- do.call(c, pixelInfos)
  values <- foreach(i = 1:length(pixelInfos), .combine=c) %do% {
    lapply(pixelInfos[[i]], function(k) makePixelInfo(name=k@name, className = k@className, patNum = i, scanNum = k@scanNum, x = k@x, y = k@y, z = k@z, tic=k@tic, ctic=k@ctic))
    }
  return(values)
}

calculateCTIC <- function(pixelInfos, data){
  for(i in 1:length(pixelInfos)){
    k <- pixelInfos[[i]]
    calcCTIC = sum(data[,i*2], na.rm=TRUE)
    pixelInfos[[i]] <- makePixelInfo(name=k@name, className = k@className, patNum = k@patNum, scanNum = k@scanNum, x = k@x, y = k@y, z = k@z, tic=k@tic, ctic=calcCTIC)
  }
  return(pixelInfos)
}
