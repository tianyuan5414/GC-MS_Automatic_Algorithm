#This script aims at matching mass spectra extracted from GC-MS files with known
  #spectral library

#Load necessary packages
{
  library(MsBackendMsp)
  library(Spectra)
  library(Matrix)
}

#Declare necessary functions
{
  #Normalize a spectrum
  normalize_spectrum <- function(intensity) {
    intensity / max(intensity)
  }
  
  #Cosine similarity function
  cosine_similarity <- function(
    mz_query,
    int_query,
    mz_ref,
    int_ref,
    mz_tol = 0.5
  ){
    
    int_ref <- int_ref / max(int_ref)
    
    q <- numeric()
    r <- numeric()
    
    for(i in seq_along(mz_query)) {
      
      hit <- which(abs(mz_ref - mz_query[i]) <= mz_tol)
      
      if(length(hit) > 0){
        
        best <- hit[
          which.min(abs(mz_ref[hit] - mz_query[i]))
        ]
        
        q <- c(q, int_query[i])
        r <- c(r, int_ref[best])
      }
    }
    
    if(length(q) < 3)
      return(0)
    
    sum(q * r) /
      sqrt(sum(q^2) * sum(r^2))
  }
  
  search_mona <- function(
    query_mz,
    query_int,
    mona_index,
    top_n = 20,
    n_filter_peaks = 5,
    min_peak_matches = 3,
    mz_tol = 1
  ){
    
    ############################################
    # Check input
    ############################################
    
    if(inherits(mona_index, "Spectra")) {
      
      stop(
        "Input is a Spectra object. ",
        "Run build_mona_index() first."
      )
      
    }
    
    ############################################
    # Query preprocessing
    ############################################
    
    query_mz <- round(query_mz)
    
    keep <- !is.na(query_mz) &
      !is.na(query_int)
    
    query_mz <- query_mz[keep]
    query_int <- query_int[keep]
    
    ############################################
    # OPTION 2
    # Candidate filtering
    ############################################
    
    peak_order <- order(
      query_int,
      decreasing = TRUE
    )
    
    filter_peaks <- query_mz[
      peak_order[
        seq_len(
          min(
            n_filter_peaks,
            length(query_mz)
          )
        )
      ]
    ]
    
    peak_index <- mona_index$peak_index
    
    candidate_logical <- vapply(
      
      peak_index,
      
      FUN.VALUE = logical(1),
      
      function(lib_mz){
        
        if(is.null(lib_mz))
          return(FALSE)
        
        if(length(lib_mz) == 0)
          return(FALSE)
        
        n_matches <- sum(
          
          vapply(
            
            filter_peaks,
            
            FUN.VALUE = logical(1),
            
            function(p){
              
              any(
                abs(lib_mz - p) <= mz_tol
              )
              
            }
            
          )
          
        )
        
        n_matches >= min(
          min_peak_matches,
          length(filter_peaks)
        )
        
      }
      
    )
    
    candidate_idx <- which(
      candidate_logical
    )
    
    if(length(candidate_idx) == 0) {
      
      message(
        "No candidate spectra found."
      )
      
      return(NULL)
    }
    
    ############################################
    # OPTION 3
    # Sparse matrix search
    ############################################
    
    mz_max <- mona_index$mz_max
    
    query_vec <- numeric(mz_max)
    
    keep <- query_mz >= 1 &
      query_mz <= mz_max
    
    mz <- query_mz[keep]
    int <- query_int[keep]
    
    if(length(mz) == 0) {
      
      stop(
        "No m/z values within mz_max."
      )
      
    }
    
    if(mona_index$weighting) {
      
      int <- sqrt(int) * mz^2
      
    }
    
    query_vec[mz] <- int
    
    query_vec <- query_vec /
      sqrt(sum(query_vec^2))
    
    scores <- as.numeric(
      
      mona_index$matrix[
        candidate_idx,
        ,
        drop = FALSE
      ] %*% query_vec
      
    )
    
    ############################################
    # Build result table
    ############################################
    
    results <- data.frame(
      SpectrumIndex = candidate_idx,
      Similarity = scores
    )
    
    results <- results[
      order(
        results$Similarity,
        decreasing = TRUE
      ),
    ]
    
    results <- head(
      results,
      top_n
    )
    
    ############################################
    # Name column detection
    ############################################
    
    meta <- mona_index$metadata
    
    name_col <- grep(
      "name",
      colnames(meta),
      ignore.case = TRUE,
      value = TRUE
    )[1]
    
    if(!is.na(name_col)) {
      
      results$Compound <-
        meta[
          results$SpectrumIndex,
          name_col
        ][[1]]
      
    }
    
    ############################################
    # Formula column detection
    ############################################
    
    formula_col <- grep(
      "formula",
      colnames(meta),
      ignore.case = TRUE,
      value = TRUE
    )[1]
    
    if(!is.na(formula_col)) {
      
      results$Formula <-
        meta[
          results$SpectrumIndex,
          formula_col
        ][[1]]
      
    }
    
    ############################################
    # Relative probability
    ############################################
    
    results$Probability <- round(
      100 *
        results$Similarity /
        sum(results$Similarity),
      2
    )
    
    rownames(results) <- NULL
    
    results
  }
  
  build_mona_index <- function(
    spectra_db,
    mz_max = 1200,
    weighting = TRUE
  ){
    
    pd <- peaksData(spectra_db)
    
    n_spec <- length(pd)
    
    lib_mat <- Matrix(
      0,
      nrow = n_spec,
      ncol = mz_max,
      sparse = TRUE
    )
    
    peak_index <- vector("list", n_spec)
    
    for(i in seq_len(n_spec)) {
      
      spec <- pd[[i]]
      
      mz <- round(spec[,1])
      int <- spec[,2]
      
      keep <- mz >= 1 & mz <= mz_max
      
      mz <- mz[keep]
      int <- int[keep]
      
      peak_index[[i]] <- mz
      
      if(length(mz) == 0)
        next
      
      if(weighting)
        int <- sqrt(int) * mz^2
      
      int <- int / sqrt(sum(int^2))
      
      lib_mat[i,mz] <- int
    }
    
    list(
      matrix = lib_mat,
      metadata = spectraData(spectra_db),
      peak_index = peak_index,
      mz_max = mz_max,
      weighting = weighting
    )
  }
}

#Data input and preparation
{
  #Raw mass spectra
  rawMSData <- read.csv(file = file.path('Output_R', 'Emre_26_GC_Peak_GC+MS.csv'))
  
  #Adjust the column names for raw  data to reflect the true m/z values
  colnames(rawMSData)[-1:-9] <- seq(50, 429, by = 0.1)
  
  #Loop to split the raw ms into list
  splitMSData <- vector(mode = 'list', length = nrow(rawMSData))
  for (loopI in 1:nrow(rawMSData)) {
    
    splitMSData[[loopI]] <- list(as.numeric(colnames(rawMSData)[-1:-9][which(rawMSData[loopI, -1:-9] > 0)]),
                                 unlist(rawMSData[loopI, -1:-9][which(rawMSData[loopI, -1:-9] > 0)]))
    
  }
  
  #Normalise the raw mass spectra
  norMSData <- t(apply(rawMSData[-1:-9], 1, normalize_spectrum))
  
  #Load MONA library
  sp <- Spectra(file.path('Data_R', 'MoNA-export-GC-MS_Spectra.msp'),
                source = MsBackendMsp())
  
  mona_index <- build_mona_index(sp)
}

#Conduct comparison
{
  #Matrix for matching results
  matchMSRes <- matrix(nrow = nrow(rawMSData),
                       ncol = 13,
                       dimnames = list(c(),
                                       c(colnames(rawMSData)[1:9],
                                         c('SpectrumIndex', 'Similarity',
                                           'Compound', 'Formula'))))
  
  matchMSRes[,1:9] <- as.matrix(rawMSData[,1:9])
  
  #Loop to fill the matrix
  for (loopI in 1:nrow(rawMSData)) {
    matchMSRes[loopI, 10:13] <-
      unlist(search_mona(
        splitMSData[[loopI]][[1]],
        normalize_spectrum(splitMSData[[loopI]][[2]]),
        mona_index,
        mz_tol = 0.5,
        top_n = 1
      ))[1:4]
  }
  
  #Export the matching results
  write.csv(matchMSRes, file = file.path('Output_R', 'matchMSRes_Emre_2026.csv'))
}