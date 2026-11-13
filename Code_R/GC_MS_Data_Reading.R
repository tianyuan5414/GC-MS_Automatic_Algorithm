#This script aims at reading original files from GC-MS, and 
  #isolate targeted signals

#Load necessary packages
{
  #Data reading
  library(here)
  
  #Pakages for GCMS data reading and analysis
  library(xcms)
  library(RColorBrewer)
  
  #Peak searching
  library(seewave)
  
  #for parallel processing
  library(foreach)
  library(doParallel)
  
  #List manipulating
  library(tidyverse)
  
  #Vector comparison
  library(dplyr)
  
  #Derivative spectra
  library(prospectr)
  
  #Judgying on invalid values
  library(gtools)
  
  #Plot
  library(ggplot2)
  library(gridExtra)
  
  #Date melt
  library(reshape2)
  
  #Spectra comparison
  library(OrgMassSpecR)
  
  #Peak searching
  library(photobiology)
}

#Declare necessary functions
{
  #Function to replace all negative values to 0, while leaving 0 and positve
    #values as they are
  nonNegative <- function(num) {
    if(num < 0){
      return(0)
    } else {
      return(1)
    }
  }
  
  #Function to find the tipping point of a spectra, used to locate the start and 
    #end retention time for a derivative spectra
  tipPoint <- function(retentionTime, deriSpectrum, oriSpectrum,
                       threshould) {
    #Vector storing all start (1) and end (2) points
    tempVec1 <- c()
    tempVec2 <- c()
    tempVec3 <- c()
    tempVec4 <- c()
    
    #Loop to find all start and end points
    loopI <- 1
    while (loopI <= (length(deriSpectrum) - 1)) {
      
      if ((deriSpectrum[loopI] == 0) && (deriSpectrum[loopI + 1] > 0)) {
        tempVec1 <- c(tempVec1, retentionTime[loopI])
        
        tempVec3 <- c(tempVec3, loopI)
      }
      
      loopI <- loopI + 1
    }
    
    loopI <- 1
    while (loopI <= (length(deriSpectrum) - 1)) {
      
      if ((deriSpectrum[loopI + 1] == 0) && (deriSpectrum[loopI] < 0)) {
        tempVec2 <- c(tempVec2, retentionTime[loopI])
        
        tempVec4 <- c(tempVec4, loopI)
      }
      
      loopI <- loopI + 1
    }
    
    #Create a matrix to store all start and end retention time
    tempMat <- matrix(nrow = length(tempVec1),
                      ncol = 4,
                      dimnames = list(c(),
                                      c('start', 'end',
                                        'start_ind', 'end_ind')))
    
    #Fill in the matrix 
    tempMat[, 1] <- tempVec1
    tempMat[, 2] <- tempVec2
    tempMat[, 3] <- tempVec3
    tempMat[, 4] <- tempVec4
    
    #Check all peaks to exclude noises that with smaller signals than threshold
    #Vector for recording qualified peak position
    tempVec5 <- c()
    
    loopI <- 1
    while (loopI <= nrow(tempMat)) {
      
      if (max(oriSpectrum[tempMat[loopI,3]:tempMat[loopI,4]]) >= 
          (threshould * max(oriSpectrum))) {
        
        tempVec5 <- c(tempVec5, loopI)
        
      }
      
      loopI <- loopI + 1
    }
    
    #Pick out the qualified peak position
    tempMat <- tempMat[tempVec5, -3:-4]
    
    tempMat2 <- tempMat
    
    #Loop to extend the peak range so that it suits more to the original spectrum
    for (loopI in 1:nrow(tempMat)) {
      
      tempMat2[loopI,] <- expandStaEnd(retentionTime,
                                       oriSpectrum,
                                       tempMat[loopI, 1],
                                       tempMat[loopI, 2],
                                       tolMax = 3)
      
    }
    
    #Merge the overlapped peak start and end
    tempMat2 <- as.matrix(merge_intervals(as.data.frame(tempMat2)))
    
    #Output results
    return(tempMat2)
  }
  
  #Expand the current start and end point of the found peak, so that it preserve most of the 
    #signals in the original spectra
  expandStaEnd <- function(retentionTime,
                           oriSpectrum, 
                           start, 
                           end,
                           tolMax = 3) {
    tempLeft <- start
    tempLeftOri <- start
    tempRight <- end
    tempRightOri <- end
    
    tempLeftSig <- oriSpectrum[which(retentionTime == tempLeft)]
    tempRightSig <- oriSpectrum[which(retentionTime == tempRight)]
    
    #Extended left boundary
    tempIndcator <- 1
    loopI <- 1
    
    while (tempIndcator == 1 && loopI <= tolMax) {
      loopJ <- 1
      while (loopJ <= 10*tolMax) {
        
        if ((which(as.numeric(retentionTime) == 
                   tempLeftOri
        ) - loopJ) < 1) {
          
          tempLeft <- retentionTime[which(as.numeric(retentionTime) == 
                                            tempLeftOri
          ) - loopJ + 1]
          
          tempIndcator <- 0
          
          
          break
          
        }else {
          
          #Search for the local minimum for left section of the current peak
          if (tempLeftSig > oriSpectrum[which(as.numeric(retentionTime) == 
                                              tempLeftOri
          ) - loopJ]) {
            
            tempLeftSig <- oriSpectrum[which(as.numeric(retentionTime) == 
                                               tempLeftOri
            ) - loopJ]
            
            tempLeft <- retentionTime[which(as.numeric(retentionTime) == 
                                              tempLeftOri
            ) - loopJ]
            
            #Reassign the terminate indicator to 1, to continue the loop
            tempIndcator <- 1
            
            loopJ <- loopJ + 1
            
          }else {
            
            loopI <- loopI + 1
            loopJ <- loopJ + 1
            
          }
        }
      }
    }
    
    #Extended right boundary
    tempIndcator <- 1
    loopI <- 1
    
    while (tempIndcator == 1 && loopI <= tolMax) {
      loopJ <- 1
      while (loopJ <= 10*tolMax) {
        
        if ((which(as.numeric(retentionTime) == 
                   tempRightOri
        ) + loopJ) > length(retentionTime)) {
          
          tempRight <- retentionTime[which(as.numeric(retentionTime) == 
                                             tempRightOri
          ) + loopJ - 1]
          
          tempIndcator <- 0
          
          break
          
        }else {
          
          #Search for the local minimum for right section of the current peak
          if (tempRightSig > oriSpectrum[which(as.numeric(retentionTime) == 
                                               tempRightOri
          ) + loopJ]) {
            
            tempRightSig <- oriSpectrum[which(as.numeric(retentionTime) == 
                                                tempRightOri
            ) + loopJ]
            
            tempRight <- retentionTime[which(as.numeric(retentionTime) == 
                                               tempRightOri
            ) + loopJ]
            
            #Reassign the terminate indicator to 1, to continue the loop
            tempIndcator <- 1
            
            loopJ <- loopJ + 1
            
          }else {
            
            loopI <- loopI + 1
            loopJ <- loopJ + 1
            
          }
        }
      }
    }
    
    
    return(c(tempLeft, tempRight))
  }
  
  #Merge the multiple ranges into inoverlapped sections
  merge_intervals <- function(df) {
    
    # Sort by start value
    df <- df[order(df$start), ]
    
    merged <- list()
    
    current_start <- df$start[1]
    current_end   <- df$end[1]
    
    for (i in 2:nrow(df)) {
      
      s <- df$start[i]
      e <- df$end[i]
      
      # Overlap detected
      if (s <= current_end) {
        
        current_end <- max(current_end, e)
        
      } else {
        
        merged[[length(merged) + 1]] <-
          data.frame(start = current_start,
                     end   = current_end)
        
        current_start <- s
        current_end   <- e
      }
    }
    
    # Add final interval
    merged[[length(merged) + 1]] <-
      data.frame(start = current_start,
                 end   = current_end)
    
    do.call(rbind, merged)
  }
  
  #Combine a list of matrix into one with combined columns
  spectra_list_to_matrix <- function(spec_list, bin_size = 0.1) {
    
    # all mz values
    all_mz <- unlist(lapply(spec_list, function(x) x[, "mz"]))
    
    min_mz <- floor(min(all_mz) / bin_size) * bin_size
    max_mz <- ceiling(max(all_mz) / bin_size) * bin_size
    
    bins <- seq(min_mz, max_mz, by = bin_size)
    
    result_mat <- matrix(
      0,
      nrow = length(spec_list),
      ncol = length(bins),
      dimnames = list(names(spec_list), sprintf("%.4f", bins))
    )
    
    for (i in seq_along(spec_list)) {
      
      spec <- as.data.frame(spec_list[[i]])
      
      # fixed-width binning
      spec$bin_idx <- floor((spec$mz - min_mz) / bin_size) + 1
      
      spec_bin <- aggregate(
        intensity ~ bin_idx,
        data = spec,
        FUN = sum
      )
      
      result_mat[i, spec_bin$bin_idx] <- spec_bin$intensity
    }
    
    result_mat
  }
}

#Read .cdf files
{
  # Get the full path to the CDF files
  cdfs <- dir(file.path('Data_R','Emre_26'))
  
  #Crate a list to store all cdf files in xcmsraw format
  gcmsRawXcms <- list()
  
  #Setup environment for paralell processing
  # For doParallel using all available cores
  registerDoParallel(cores = 14)
  
  #Loop to read all cdfs as xcmsraw format
  gcmsRawXcms <- foreach(i = 1:length(cdfs)) %dopar% {
    xcms::xcmsRaw(file.path('Data_R', 'Emre_26',cdfs[i]), profstep = 0.1)
  }
  
  # For doParallel
  stopImplicitCluster() 
  
  if (!file.exists(file.path('Output_R', 'Emre_26.RData'))) {
    #Save the spectra data read from CDF
    save(gcmsRawXcms, file = file.path('Output_R', 'Emre_26.RData'))
  }
}

#Extract signals of given Mass-to-charge ratio
{
  # #Create a matrix to store the GC spectra under given m/z
  # eicmz <- c(161, 196)
  # eicRange <- c(0.3, 0.7)
  timeRange <- c(4, 28) *60
  
  # #Prepare a matrix for extracting the EIC
  # eicMat <- matrix(nrow = length(eicmz), ncol = 2,
  #                  dimnames = list(c(), c('mzmin', 'mzmax')))
  # 
  # #Calculate the m/z range for each interested value
  # eicMat[, 1] <- eicmz - eicRange[1]
  # eicMat[, 2] <- eicmz + eicRange[2]
  
  #Create a list to store all eic
  gcEIC <- list()
  
  #Loop to Get the EIC
  loopI <- 1
  while(loopI <= length(gcmsRawXcms)) {
    
    gcEIC[[loopI]] <- getEIC(gcmsRawXcms[[loopI]], 
                             # eicMat, 
                             step = 1, 
                             rtrange = timeRange)@eic[["xcmsRaw"]]
    
    loopI <- loopI + 1
  }
  
  #Define the s/n to clean the noise in the EIC
  snRatio <- 30
  
  #Copy a list of gcEIC to store the cleanned EIC
  gcEICCleaned <- gcEIC
  
  #Clean up noises in EIC
  loopI <- 1
  while (loopI <= length(gcEIC)) {
    
    #Loop to go through every interested m/z
    loopJ <- 1
    while (loopJ <= length(gcEIC[[loopI]])) {
      
      #Calculate the background noises under given s/n
      tempBackNoise <- quantile(gcEIC[[loopI]][[loopJ]][, 2], probs = 1- 1/snRatio)
      
      #exclude all noises according to given s/n, referring to the corresponding
        #quantile as baseline value
      gcEICCleaned[[loopI]][[loopJ]] <- sapply(gcEIC[[loopI]][[loopJ]][,2] - tempBackNoise,
                                     FUN = nonNegative) * gcEIC[[loopI]][[loopJ]][,2]
      
      loopJ <- loopJ + 1
    }
    loopI <- loopI + 1
  }
  
  #Calculate the derivative spectra basing on the cleaned EIC spectra
  #Create a list for derivative spectra
  gcEICDeri <- gcEICCleaned
  gcEICRetention <- gcEICCleaned
  
  loopI <- 1
  while (loopI <= length(gcEICCleaned)) {
    
    #Loop to go through all EICs of every interested m/z
    loopJ <- 1
      while (loopJ <= length(gcEICCleaned[[loopI]])) {
        
        gcEICDeri[[loopI]][[loopJ]] <- diff(gcEICCleaned[[loopI]][[loopJ]], 
                                            differences = 1)
        
        #Find the start and end retention time of all peaks for targeted m/z
        gcEICRetention[[loopI]][[loopJ]] <-
          tipPoint(gcEIC[[loopI]][[loopJ]][-1, 1],
                   gcEICDeri[[loopI]][[loopJ]],
                   gcEIC[[loopI]][[loopJ]][-1, 2],
                   0.01)
        
        loopJ <- loopJ + 1
      }
    
    loopI <- loopI + 1
  }
  
  #Create a list to record peak statistics
  gcmsPeakStatis <- matrix(nrow = 1, ncol = 9,
                           dimnames = list(c(), 
                                           c('No',
                                             'mz',
                                             'mzmin',
                                             'mzmax',
                                             'rt',
                                             'rtmin',
                                             'rtmax',
                                             'into',
                                             'maxo')))
  
  #Loop to extract the peak information of all interested m/z
  loopI <- 1
  while (loopI <= length(gcEICRetention)) {
    
    #Create a matrix for storing all peak statistics for current sample
    tempMat <- matrix(nrow = 1, ncol = 4,
                       dimnames = list(c(),
                                       c('mzmin', 'mzmax',
                                         'rtmin', 'rtmax')))
    
    loopJ <- 1
    # while (loopJ <= nrow(eicMat)) {
    while (loopJ <= 1) {
      
      if (!invalid(gcEICRetention[[loopI]][[loopJ]])) {
        
        if (is.vector(gcEICRetention[[loopI]][[loopJ]])) {
          #Create a matrix for recording the m/z range and retention time range
          tempMat2 <- matrix(nrow = 1, ncol = 4,
                            dimnames = list(c(),
                                            c('mzmin', 'mzmax',
                                              'rtmin', 'rtmax')))
          
          #Fill in the m/z and retention time ranges
          # tempMat2[1, 1:2] <- eicMat[loopJ,]
          tempMat2[1, 3:4] <- gcEICRetention[[loopI]][[loopJ]]
          
        }else {
          
          tempMat2 <- matrix(nrow = nrow(gcEICRetention[[loopI]][[loopJ]]), ncol = 4,
                            dimnames = list(c(),
                                            c('mzmin', 'mzmax',
                                              'rtmin', 'rtmax')))
          
          # tempMat2[, 1:2] <-
          #   matrix(rep(eicMat[loopJ,], nrow(gcEICRetention[[loopI]][[loopJ]])),
          #          nrow = nrow(gcEICRetention[[loopI]][[loopJ]]),
          #          ncol = length(eicMat[loopJ,]),
          #          byrow = TRUE)
          
          tempMat2[, 3:4] <- gcEICRetention[[loopI]][[loopJ]]
          
        }
      }
      
      #Combine the current peak information matrix with the previous one
      tempMat <- rbind(tempMat, tempMat2)
      
      loopJ <- loopJ + 1
    }
    
    #Delete the first row of the information matrix
    tempMat <- tempMat[-1, ]
    
    #Get the peak statistics for current sample
    gcmsPeakStatis <- rbind(gcmsPeakStatis,
                            cbind(rep(loopI, nrow(tempMat)),
                                  getPeaks(gcmsRawXcms[[loopI]], tempMat)))
    
    loopI <- loopI + 1
  }
  
  #Remove the additional row of the last matrix
  gcmsPeakStatis <- gcmsPeakStatis[-1, ]
  
  if (!file.exists(file.path('Output_R', 'Emre_26_GC_Peak.csv'))) {
    
    #Output the peak statistics table
    write.csv(gcmsPeakStatis, file = file.path('Output_R', 'Emre_26_GC_Peak.csv'))
    
  }
}

#Extract all mass spectra at all peaks identified
{
  #Create a list to store all mass spectra
  gcmsMS <- vector(mode = 'list', length = length(gcmsRawXcms))
  
  for (loopI in 1:length(gcmsMS)) {
    
    #Get the number of peaks identified
    tempNumPeaks <- nrow(gcmsPeakStatis[gcmsPeakStatis[, 1] == loopI, ])
    gcmsMS[[loopI]] <- vector(mode = 'list', length = tempNumPeaks)
    
    #Go through all individual peaks in each sample
    for (loopJ in 1:tempNumPeaks) {
      
      gcmsMS[[loopI]][[loopJ]] <- getSpec(gcmsRawXcms[[loopI]], 
                                          rtrange = gcmsPeakStatis[gcmsPeakStatis[, 1] == loopI, ][loopJ, 5])
      
    }
  }
  
  #Unlist the second structure in the last list
  gcmsMSUnSec <- vector(mode = 'list', length = nrow(gcmsPeakStatis))
  tempInd <- 0
  for (loopI in 1:length(gcmsMS)) {
    
    gcmsMSUnSec[(tempInd + 1):(tempInd + length(gcmsMS[[loopI]]))] <-
      gcmsMS[[loopI]]
    
    tempInd <- tempInd + length(gcmsMS[[loopI]])
    
  }
  
  #Combine all MS matrix together, following the same order as the GC peak matrix
  gcmsMSMat <- spectra_list_to_matrix(gcmsMSUnSec)
  
  gcmsFullInfo <- as.data.frame(cbind(gcmsPeakStatis, gcmsMSMat))
  
  if (!file.exists(file.path('Output_R', 'Emre_26_GC_Peak_GC+MS.csv'))) {
    
    #Read batch information
    gcmsBatch <- read.csv(file.path('Data_R', 'Emre_26_info.csv'))
    
    # Match No
    idx <- match(gcmsFullInfo$No, gcmsBatch$No)
    
    gcmsFullInfo <- cbind(No = gcmsFullInfo[, 1], 
                          file_id = gcmsBatch$file_id[idx],
                          batch_id = gcmsBatch$batch_id[idx],
                          gcmsFullInfo[, -1])
    
    # Save result
    write.csv(
      gcmsFullInfo[, -4:-6],
      file.path('Output_R',
                "Emre_26_GC_Peak_GC+MS.csv"),
      row.names = FALSE
    )
  }
}

#Data output
{
  #Export the necessary results
  #Peak ratio
  write.csv(gcmsPeakRatio, file = paste(here(), '/gcmsPeakRatio.csv',
                                        sep = '',
                                        collapse = NULL))
}

#Code for presentation on the results
{
  #Showing the functionality of m/z range selection
  {
    #Prescribe the m/z ratio and retention time for example
    eicmz <- c(161)
    timeRange <- c(4, 28) *60
    
    mzRangeSel <- list()
    
    mzRangeSelData <- matrix(nrow = 15, ncol = 2)
    
    mzRangeSelData[1:10, 1] <- seq(0.1, 1, by = 0.1)
    
    mzRangeSelData[1:10, 2] <- 1 - seq(0.1, 1, by = 0.1)
    
    mzRangeSelData[11:15, 1] <- seq(0.1, 0.5, by = 0.1)
    
    mzRangeSelData[11:15, 2] <- 0.5 - seq(0.1, 0.5, by = 0.1)
    
    #Loop to go through all combinations of m/z boundaries
    loopI <- 1
    while (loopI <= nrow(mzRangeSelData)) {
      
      #Create a matrix to store the GC spectra under given m/z
      eicRange <- mzRangeSelData[loopI,]
      
      #Prepare a matrix for extracting the EIC
      eicMat <- matrix(nrow = length(eicmz), ncol = 2,
                       dimnames = list(c(), c('mzmin', 'mzmax')))
      
      #Calculate the m/z range for each interested value
      eicMat[, 1] <- eicmz - eicRange[1]
      eicMat[, 2] <- eicmz + eicRange[2]
      
      #Create a list to store all eic
      testa <- list()
      
      #Loop to Get the EIC
      loopJ <- 1
      while(loopJ <= length(gcmsRawXcms)) {
        
        testa[[loopJ]] <- getEIC(gcmsRawXcms[[loopJ]], eicMat, step = 1, 
                                 rtrange = timeRange)@eic[["xcmsRaw"]]
        
        loopJ <- loopJ + 1
      }
      
      mzRangeSel[[loopI]] <- testa[[2]][[1]]
      
      loopI <- loopI + 1
    }
    
    #Visualise the example EIC under different settings
    #Rearrange all data into one matrix
    #Set up an offset for differnt lines
    offset <- 10000
    
    #The matrix for storing all data for visualisation
    mzRangeSelMat <- matrix(ncol = length(mzRangeSel) + 1,
                            nrow = nrow(mzRangeSel[[1]]))
    
    mzRangeSelMat[, 1] <- mzRangeSel[[1]][, 1]
    

    
    for (loopI in 1:length(mzRangeSel)) {
      mzRangeSelMat[,loopI + 1] <- mzRangeSel[[loopI]][, 2] +
        offset * (loopI - 1)
    }
    
    #Prepare a vector for the group information
    mzRangeSelGroup <- paste(-mzRangeSelData[, 1], '~', mzRangeSelData[, 2], 
                             sep = '')
    
    #Add the group name to the matrix, and change the matrix's column names,
      #and transfer it into a data frame
    mzRangeSelMat <- t(rbind(c(NA, mzRangeSelGroup), mzRangeSelMat))
    
    colnames(mzRangeSelMat) <- mzRangeSelMat[1, ]
    
    mzRangeSelMat <- mzRangeSelMat[-1, ]
    
    mzRangeSelMat <- as.data.frame(mzRangeSelMat)
    
    colnames(mzRangeSelMat)[1] <- 'Range'
    
    
    #Melt the data matrix
    mzRangeSelMatMelt <- 
      melt(mzRangeSelMat, id.vars = "Range", variable.name = "Time", 
                      value.name = "EIC")
    
    mzRangeSelMatMelt[, 1] <- factor(mzRangeSelMatMelt$Range, 
                                     levels = rev(mzRangeSelGroup))
    
    mzRangeSelMatMelt[, 2:3] <- apply(mzRangeSelMatMelt[, 2:3],
                                      c(1,2),
                                      as.numeric)
    
    #Plot
    mzRangeSelPlot <-
      ggplot(mzRangeSelMatMelt, aes(x = Time, y = EIC, color = Range)) +
      geom_line(size = 0.4) +
      labs(title = "Different range of m/z around 161",
           x = "Time (s)",
           y = "",
           color = "m/z range") +
      theme_classic() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
      )
    
    
  }
  
  #Showing the influence of S/N ratio on filtered spectra
  {
    #Define the list of S/N
    snRatioSelData <- seq(100, 600, by = 100)
    
    #Prepare a list for storing all filtered example spectra
    snRatioSel <- list()
    
    #Loop to extract all data under different S/N
    for (loopI in 1:length(snRatioSelData)) {
      
      #Define the s/n to clean the noise in the EIC
      snRatio <- snRatioSelData[loopI]
      
      #Copy a list of gcEIC to store the cleanned EIC
      testb <- gcEIC
      
      #Clean up noises in EIC
      for (loopJ in 1:length(gcEIC)) {
        
        for (loopK in 1:length(gcEIC[[loopJ]])) {
          
          #Calculate the background noises under given s/n
          tempBackNoise <- quantile(gcEIC[[loopJ]][[loopK]][, 2], probs = 1- 1/snRatio)
          
          #exclude all noises according to given s/n, referring to the corresponding
          #quantile as baseline value
          testb[[loopJ]][[loopK]] <- sapply(gcEIC[[loopJ]][[loopK]][,2] - tempBackNoise,
                                                   FUN = nonNegative) * gcEIC[[loopJ]][[loopK]][,2]
        }
      }
      
      snRatioSel[[loopI]] <- testb[[2]][[1]]
    }
    
    #Loop to plot all test results into seperate plots
    for (loopI in 1:length(snRatioSel)) {
      
      #Plot
      snRatioRawPlot <-
        ggplot() +
        geom_line(aes(x = as.numeric(testa[[2]][[1]][, 1]), y = as.numeric(testa[[2]][[1]][, 2])),
                  size = 1) +
        geom_hline(yintercept = quantile(testa[[2]][[1]][, 2], probs = 1- 1/snRatioSelData[[loopI]]),
                   size  = 0.5, color = 'red')+
        labs(title = paste("Original EIC + signal to noise level at ", snRatioSelData[loopI],
                           sep = '', collapse = NULL),
             x = "Time (s)",
             y = "") +
        theme_classic() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
        )
      
      #Plot
      snRatioPlot <-
        ggplot() +
        geom_line(aes(x = as.numeric(testa[[2]][[1]][, 1]), 
                      y = as.numeric(snRatioSel[[loopI]])),
                  size = 1) +
        labs(title = paste("Signal to noise level = ", snRatioSelData[loopI],
                           sep = '', collapse = NULL),
             x = "Time (s)",
             y = "") +
        theme_classic() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
        )
      
      #Create a pdf for the current plot
      pdf(paste("SN", snRatioSelData[loopI], '.pdf',
                sep = '', collapse = NULL), width = 16, height = 7)
      
      grid.arrange(snRatioRawPlot, snRatioPlot)
      
      #Close the PDF device
      dev.off()
    }
  }
  
  #Showing processes of find peaks through derivative spectra
  {
    #Compare the cleaned EIC and its derivative spectra
    {
      #Plot
      peakSelRawPlot <-
        ggplot() +
        geom_line(aes(x = as.numeric(testa[[2]][[1]][, 1]), y = as.numeric(gcEICCleaned[[2]][[1]])),
                  size = 1) +
        labs(title = "Cleaned EIC",
             x = "Time (s)",
             y = "") +
        theme_classic() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
        )
      
      #Plot
      peakSelDeriPlot <-
        ggplot() +
        geom_line(aes(x = as.numeric(testa[[2]][[1]][-1, 1]), y = as.numeric(gcEICDeri[[2]][[1]])),
                  size = 1) +
        labs(title = '1st derivative of EIC',
             x = "Time (s)",
             y = "") +
        theme_classic() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
        )
      
      #Create a pdf for the current plot
      pdf('Cleaned_Derivative.pdf', width = 16, height = 7)
      
      grid.arrange(peakSelRawPlot, peakSelDeriPlot)
      
      #Close the PDF device
      dev.off()
    }
    
    #Compare the position of identified peaks
    {
      #Melt the original peak statistic matrix for plotting
      peakSelSta <- gcmsPeakStatis[which(gcmsPeakStatis[, 3] == 160.7),c(1,5,6,7)]
      
      peakSelSta <- as.data.frame(cbind(peakSelSta, seq(1, nrow(peakSelSta), by = 1)))
      
      #Plot
      peakSelStaPlot <-
        ggplot(data = peakSelSta,
               aes(group = V5, y = No)) +
        geom_ribbon(aes(xmin = rtmin, xmax = rtmax), 
                    fill = "black", alpha = 0.5, orientation = "y") +
        scale_x_continuous(limits = timeRange) +
        labs(title = "Peak position for pCA on gas chromatography",
             x = "Time (s)",
             y = "") +
        theme_classic() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
        )
      
      #Density plot for central retention time
      peakSelStaDenPlot <-
        ggplot(data = peakSelSta,
               aes(x = rt)) +
        geom_density()+
        scale_x_continuous(limits = timeRange) +
        labs(title = "Distributional density - retention time of local signal maxima",
             x = "Time (s)",
             y = "") +
        theme_classic() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
        )
      
      #Create a pdf for the current plot
      pdf('Peak_Statistics.pdf', width = 16, height = 7)
      
      grid.arrange(peakSelStaPlot, peakSelStaDenPlot)
      
      #Close the PDF device
      dev.off()
      
    }
  }
  
  #Showing processes of searching for mass spectra for given retention time
  {
    #Get the example mass spectra
    testc <- gcCandi[[2]][[1]]
    
    #Plot the mass spectra (multiple plot if possible)
    #Loop to go through all plots
    massSpecSelPlot <- list()
    
    for (loopI in 1:length(testc)) {
      
      massSpecSelPlot[[loopI]] <-
        ggplot(data = testc[[loopI]],
               aes(x=mz, ymax=intensity, ymin=0)) +
        geom_linerange() +
        scale_x_continuous(limits = c(40, 210),
                           expand=c(0,0),
                           breaks = seq(40, 210, by = 10)) +
        scale_y_continuous(limits = c(0, 26000),
                           expand=c(0,0)) +
        labs(title = paste("Mass spectra of peak ",
                           loopI, sep = '',
                           collapse = NULL),
             x = "m/z",
             y = "Intensity") +
        theme_classic() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
        )
      
      #Plot the mass spectra together with the gc spectra
      #Create a pdf for the current plot
      pdf(paste('Mass_Spectra_', loopI,
                '.pdf', sep = '',
                collapse = NULL), width = 16, height = 7)
      
      grid.arrange(peakSelRawPlot, massSpecSelPlot[[loopI]])
      
      #Close the PDF device
      dev.off()
    }
    
    #Compare the example spectra with the reference spectra
    #Loop to go through all plots
    massSpecSelComPlot <- list()
    
    for (loopI in 1:length(testc)) {
      
      massSpecSelComPlot[[loopI]] <-
        SpectrumSimilarity(as.data.frame(testc[[loopI]]), 
                           msReferenceSpectra[[1]],
                           xlim = c(50, 210),
                           output.list = TRUE)
      #Plot the mass spectra together with the gc spectra
      #Create a pdf for the current plot
      pdf(paste('Mass_Spectra_Compare_', loopI,
                '.pdf', sep = '',
                collapse = NULL), width = 16, height = 7)
      
      grid.arrange(peakSelRawPlot, massSpecSelComPlot[[loopI]]$plot)
      
      #Close the PDF device
      dev.off()
    }
  }
}