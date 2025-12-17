########################################################################
# CLEAN FUNCTIONS FOR TOPOLOGICAL SUMMARY COMPARISON
# Comparing: VAB, PL, PI, HNAV, HWNAV, OW-HNPV

########################################################################

library(igraph)
library(TDA)
library(tidyverse)
library(depthTools)
library(randomForest)
library(caret)
library(fda.usc)
library(pROC)
library(lubridate)
library(ggplot2)

options(dplyr.summarise.inform = FALSE)

########################################################################
# 1. GRAPH FEATURES
########################################################################

featureGraph <- function(networkDF, graphDir) {
  
  periodList <- unique(networkDF$time)
  
  featurePerDay <- function(inputDay) {
    periodDF <- networkDF %>% filter(time == inputDay)
    G <- graph.data.frame(periodDF[, c("from", "to")], directed = FALSE)
    
    tibble(
      day = inputDay,
      vertexNum = vcount(G),
      edgeNum = ecount(G),
      clusterCoef = transitivity(G)
    )
  }
  
  res <- map_dfr(periodList, featurePerDay)
  saveRDS(res, file.path(graphDir, "graphFeature_allTokens.rds"))
  cat("Graph features saved!\n")
  return(res)
}

########################################################################
# 2. PERSISTENCE DIAGRAM COMPUTATION
########################################################################

computePD <- function(networkDF, topRank, pdDir, filtration = "sublevel") {
  
  selectTop <- function(dat, topRank) {
    frqFrom <- table(dat$from) %>% sort(decreasing = TRUE) %>% names() %>% .[1:topRank]
    dat <- subset(dat, from %in% frqFrom)
    frqTo <- table(dat$to) %>% sort(decreasing = TRUE) %>% names() %>% setdiff(frqFrom) %>% .[1:topRank]
    topList <- union(frqFrom, frqTo)
    dat <- subset(dat, (to %in% topList) & (from %in% topList))
    return(dat)
  }
  
  periodList <- unique(networkDF$time)
  nDays <- length(periodList)
  
  for (i in 1:nDays) {
    
    periodDF <- networkDF %>% filter(time == periodList[i])
    periodDF <- selectTop(periodDF, topRank)
    periodDF <- periodDF %>% 
      group_by(from, to) %>% 
      summarise(value = mean(value), .groups = 'drop')
    
    if (nrow(periodDF) == 0) {
      cat(paste("Warning: No data for day", periodList[i], "- skipping\n"))
      next
    }
    
    G <- graph.data.frame(periodDF[, c("from", "to")], directed = FALSE)
    G <- igraph::simplify(G)
    
    # Compute simplices
    zeroSimplx <- as.list(V(G))
    oneSimplx <- data.frame(t(as_edgelist(G, names = FALSE)))
    twoSimplx <- data.frame(matrix(triangles(G), nrow = 3))
    cmplx <- c(zeroSimplx, oneSimplx, twoSimplx)
    
    # Node values
    vertValues <- matrix(0, nrow = length(V(G)), ncol = 2)
    rownames(vertValues) <- unique(c(periodDF$from, periodDF$to))
    colnames(vertValues) <- c('sent', 'received')
    
    valuesFrom <- tapply(periodDF$value, periodDF$from, mean)
    vertValues[names(valuesFrom), 1] <- valuesFrom
    
    valuesTo <- tapply(periodDF$value, periodDF$to, mean)
    vertValues[names(valuesTo), 2] <- valuesTo
    
    vertValues <- rowSums(vertValues)
    if (max(vertValues) > 0) {
      vertValues <- vertValues / max(vertValues)
    }
    
    # Filtration
    Flt <- funFiltration(FUNvalues = vertValues, 
                         cmplx = cmplx, 
                         sublevel = (filtration == 'sublevel'))
    
    PD <- filtrationDiag(filtration = Flt, maxdimension = 1, 
                         library = 'Dionysus', location = TRUE)$diagram
    
    # Replace infinity
    PD[PD[, 3] == Inf, 3] <- 1.01
    
    saveRDS(PD, file = file.path(pdDir, paste0('PD_', filtration, '_', periodList[i], '.rds')))
    
    if (i %% 30 == 0) cat(paste(i, 'out of', nDays, 'PDs processed\n'))
    
    # Memory cleanup
    rm(periodDF, G, zeroSimplx, oneSimplx, twoSimplx, cmplx, vertValues, Flt, PD)
  }
  
  gc(verbose = FALSE)
  cat("All persistence diagrams computed!\n")
}

########################################################################
# 3. VAB (Vector of Averaged Bettis)
########################################################################

computeVAB <- function(networkDF, pdDir, vabDir, num_intervals, filtration = "sublevel") {
  
  scale_seq <- seq(0, 1.01, length.out = num_intervals + 1)
  
  extractBetti <- function(D) {
    if (nrow(D) == 0) return(rep(0, num_intervals))
    
    x <- D[, 1]  # birth
    y <- D[, 2]  # death
    delta <- diff(scale_seq)
    n <- length(delta)
    betti <- numeric(n)
    
    for (k in 1:n) {
      b <- pmin(scale_seq[k + 1], y) - pmax(scale_seq[k], x)
      betti[k] <- sum(pmax(0, b)) / delta[k]
    }
    return(betti)
  }
  
  periodList <- unique(networkDF$time)
  nDays <- length(periodList)
  
  B0 <- B1 <- matrix(0, nrow = nDays, ncol = num_intervals)
  
  for (i in 1:nDays) {
    pd_file <- file.path(pdDir, paste0('PD_', filtration, '_', periodList[i], '.rds'))
    if (!file.exists(pd_file)) next
    
    PD <- readRDS(pd_file)
    B0[i, ] <- extractBetti(PD[PD[, 1] == 0, 2:3, drop = FALSE])
    B1[i, ] <- extractBetti(PD[PD[, 1] == 1, 2:3, drop = FALSE])
    
    if (i %% 30 == 0) cat(paste(i, 'out of', nDays, 'VAB processed\n'))
  }
  
  colnames(B0) <- paste0('int', 1:num_intervals)
  colnames(B1) <- paste0('int', 1:num_intervals)
  
  # NOTE: Column names are "vab0" and "vab1" (no trailing underscore)
  saveRDS(data.frame(Time = paste0('vab0.', periodList), B0), 
          file.path(vabDir, 'allTokens_vab0.rds'))
  saveRDS(data.frame(Time = paste0('vab1.', periodList), B1), 
          file.path(vabDir, 'allTokens_vab1.rds'))
  
  cat("VAB computation complete!\n")
}

########################################################################
# 4. PERSISTENCE LANDSCAPES (PL)
########################################################################

computePL <- function(networkDF, pdDir, plDir, num_intervals, filtration = "sublevel") {
  
  scale_seq <- seq(0, 1.01, length.out = num_intervals)
  
  periodList <- unique(networkDF$time)
  nDays <- length(periodList)
  
  PL0 <- PL1 <- matrix(0, nrow = nDays, ncol = num_intervals)
  
  for (i in 1:nDays) {
    pd_file <- file.path(pdDir, paste0('PD_', filtration, '_', periodList[i], '.rds'))
    if (!file.exists(pd_file)) next
    
    PD <- readRDS(pd_file)
    PL0[i, ] <- landscape(PD, dimension = 0, KK = 1, tseq = scale_seq)
    PL1[i, ] <- landscape(PD, dimension = 1, KK = 1, tseq = scale_seq)
    
    if (i %% 30 == 0) cat(paste(i, 'out of', nDays, 'PL processed\n'))
  }
  
  colnames(PL0) <- paste0('int', 1:num_intervals)
  colnames(PL1) <- paste0('int', 1:num_intervals)
  
  saveRDS(data.frame(Time = paste0('pl0.', periodList), PL0), 
          file.path(plDir, 'allTokens_pl0.rds'))
  saveRDS(data.frame(Time = paste0('pl1.', periodList), PL1), 
          file.path(plDir, 'allTokens_pl1.rds'))
  
  cat("PL computation complete!\n")
}

########################################################################
# 5. PERSISTENCE IMAGES (PI)
########################################################################

computePI <- function(networkDF, pdDir, piDir, res = 6, filtration = "sublevel") {
  
  PI_func <- function(D, res, sig, maxB, maxP) {
    if (nrow(D) == 0) return(rep(0, res^2))
    
    dy <- maxP / res
    y_lower <- seq(0, maxP - dy, length.out = res)
    y_upper <- y_lower + dy
    
    dx <- maxB / res
    x_lower <- seq(0, maxB - dx, length.out = res)
    x_upper <- x_lower + dx
    
    PSurfaceHk <- function(point) {
      x <- point[1]
      y <- point[2]
      out1 <- pnorm(x_upper, mean = x, sd = sig) - pnorm(x_lower, mean = x, sd = sig)
      out2 <- pnorm(y_upper, mean = y, sd = sig) - pnorm(y_lower, mean = y, sd = sig)
      wgt <- y / maxP * (y < maxP) + 1 * (y >= maxP)
      return(out1 %o% out2 * wgt)
    }
    
    Psurf_mat <- apply(D, 1, PSurfaceHk)
    return(rowSums(Psurf_mat))
  }
  
  periodList <- unique(networkDF$time)
  nDays <- length(periodList)
  
  # First pass: find max values
  maxB0 <- maxB1 <- maxP0 <- maxP1 <- numeric(nDays)
  PD_list <- list()
  
  for (i in 1:nDays) {
    pd_file <- file.path(pdDir, paste0('PD_', filtration, '_', periodList[i], '.rds'))
    if (!file.exists(pd_file)) {
      PD_list[[i]] <- matrix(0, nrow = 0, ncol = 3)
      next
    }
    
    pd <- readRDS(pd_file)
    pd[, 3] <- pd[, 3] - pd[, 2]  # Convert to persistence
    
    if (sum(pd[, 1] == 0) > 0) {
      maxB0[i] <- max(pd[pd[, 1] == 0, 2])
      maxP0[i] <- max(pd[pd[, 1] == 0, 3])
    }
    if (sum(pd[, 1] == 1) > 0) {
      maxB1[i] <- max(pd[pd[, 1] == 1, 2])
      maxP1[i] <- max(pd[pd[, 1] == 1, 3])
    }
    PD_list[[i]] <- pd
  }
  
  # Second pass: compute PIs
  PI0 <- PI1 <- matrix(0, nrow = nDays, ncol = res^2)
  sigH0 <- 0.5 * max(maxP0, na.rm = TRUE) / res
  sigH1 <- 0.5 * max(maxP1, na.rm = TRUE) / res
  
  # Handle edge cases
  if (sigH0 == 0 || is.na(sigH0)) sigH0 <- 0.01
  if (sigH1 == 0 || is.na(sigH1)) sigH1 <- 0.01
  
  for (i in 1:nDays) {
    pd <- PD_list[[i]]
    if (nrow(pd) == 0) next
    
    PI0[i, ] <- PI_func(pd[pd[, 1] == 0, 2:3, drop = FALSE], res, sigH0, 
                        max(maxB0, na.rm = TRUE), max(maxP0, na.rm = TRUE))
    PI1[i, ] <- PI_func(pd[pd[, 1] == 1, 2:3, drop = FALSE], res, sigH1, 
                        max(maxB1, na.rm = TRUE), max(maxP1, na.rm = TRUE))
    
    if (i %% 30 == 0) cat(paste(i, 'out of', nDays, 'PI processed\n'))
  }
  
  saveRDS(data.frame(Time = paste0('pi0.', periodList), PI0), 
          file.path(piDir, 'allTokens_pi0.rds'))
  saveRDS(data.frame(Time = paste0('pi1.', periodList), PI1), 
          file.path(piDir, 'allTokens_pi1.rds'))
  
  cat("PI computation complete!\n")
}

########################################################################
# 6. HNAV (Hierarchical Normalized Averaged Velocity) - Unweighted
########################################################################

computeHNAV <- function(networkDF, pdDir, hnavDir, m, n_sub, filtration = "sublevel") {
  
  extractHNAV <- function(D, alpha = 0, beta = 1.01) {
    if (nrow(D) == 0) return(rep(0, m))
    
    births <- D[, 1]
    deaths <- D[, 2]
    n_features <- nrow(D)
    
    delta_s <- (beta - alpha) / m
    s_points <- seq(alpha, beta, by = delta_s)
    
    hnav_vec <- numeric(m)
    
    for (j in 1:m) {
      sj <- s_points[j]
      sj1 <- s_points[j + 1]
      
      delta_t <- (sj1 - sj) / n_sub
      t_points <- seq(sj, sj1, by = delta_t)
      
      V_j <- 0
      for (ell in 1:n_sub) {
        t_ell <- t_points[ell]
        t_ell1 <- t_points[ell + 1]
        
        N_b <- sum(births >= t_ell & births < t_ell1)
        N_d <- sum(deaths >= t_ell & deaths < t_ell1)
        
        delta_t_ell <- t_ell1 - t_ell
        if (delta_t_ell > 0) {
          V_ell <- (N_b + N_d) / (2 * delta_t_ell)
        } else {
          V_ell <- 0
        }
        V_j <- V_j + V_ell
      }
      
      if (n_features > 0) {
        hnav_vec[j] <- V_j / (n_sub * n_features)
      }
    }
    
    return(hnav_vec)
  }
  
  periodList <- unique(networkDF$time)
  nDays <- length(periodList)
  
  HNAV0 <- HNAV1 <- matrix(0, nrow = nDays, ncol = m)
  
  for (i in 1:nDays) {
    pd_file <- file.path(pdDir, paste0('PD_', filtration, '_', periodList[i], '.rds'))
    if (!file.exists(pd_file)) next
    
    PD <- readRDS(pd_file)
    HNAV0[i, ] <- extractHNAV(PD[PD[, 1] == 0, 2:3, drop = FALSE])
    HNAV1[i, ] <- extractHNAV(PD[PD[, 1] == 1, 2:3, drop = FALSE])
    
    if (i %% 30 == 0) cat(paste(i, 'out of', nDays, 'HNAV processed\n'))
  }
  
  colnames(HNAV0) <- paste0('int', 1:m)
  colnames(HNAV1) <- paste0('int', 1:m)
  
  saveRDS(data.frame(Time = paste0('hnav0.', periodList), HNAV0), 
          file.path(hnavDir, 'allTokens_hnav0.rds'))
  saveRDS(data.frame(Time = paste0('hnav1.', periodList), HNAV1), 
          file.path(hnavDir, 'allTokens_hnav1.rds'))
  
  cat("HNAV computation complete!\n")
}

########################################################################
# 7. HWNAV (Hierarchical Weighted Normalized Averaged Velocity)
########################################################################

computeHWNAV <- function(networkDF, pdDir, hwnavDir, m, n_sub, filtration = "sublevel") {
  
  extractHWNAV <- function(D, alpha = 0, beta = 1.01) {
    if (nrow(D) == 0) return(rep(0, m))
    
    births <- D[, 1]
    deaths <- D[, 2]
    persistence <- deaths - births
    
    P_total <- sum(persistence)
    if (P_total == 0) return(rep(0, m))
    
    delta_s <- (beta - alpha) / m
    s_points <- seq(alpha, beta, by = delta_s)
    
    hwnav_vec <- numeric(m)
    
    for (j in 1:m) {
      sj <- s_points[j]
      sj1 <- s_points[j + 1]
      
      delta_t <- (sj1 - sj) / n_sub
      t_points <- seq(sj, sj1, by = delta_t)
      
      V_j_w <- 0
      for (ell in 1:n_sub) {
        t_ell <- t_points[ell]
        t_ell1 <- t_points[ell + 1]
        
        W_b <- sum(persistence[births >= t_ell & births < t_ell1])
        W_d <- sum(persistence[deaths >= t_ell & deaths < t_ell1])
        
        delta_t_ell <- t_ell1 - t_ell
        if (delta_t_ell > 0) {
          V_ell_w <- (W_b + W_d) / (2 * delta_t_ell)
        } else {
          V_ell_w <- 0
        }
        V_j_w <- V_j_w + V_ell_w
      }
      
      hwnav_vec[j] <- V_j_w / (n_sub * P_total)
    }
    
    return(hwnav_vec)
  }
  
  periodList <- unique(networkDF$time)
  nDays <- length(periodList)
  
  HWNAV0 <- HWNAV1 <- matrix(0, nrow = nDays, ncol = m)
  
  for (i in 1:nDays) {
    pd_file <- file.path(pdDir, paste0('PD_', filtration, '_', periodList[i], '.rds'))
    if (!file.exists(pd_file)) next
    
    PD <- readRDS(pd_file)
    HWNAV0[i, ] <- extractHWNAV(PD[PD[, 1] == 0, 2:3, drop = FALSE])
    HWNAV1[i, ] <- extractHWNAV(PD[PD[, 1] == 1, 2:3, drop = FALSE])
    
    if (i %% 30 == 0) cat(paste(i, 'out of', nDays, 'HWNAV processed\n'))
  }
  
  colnames(HWNAV0) <- paste0('int', 1:m)
  colnames(HWNAV1) <- paste0('int', 1:m)
  
  saveRDS(data.frame(Time = paste0('hwnav0.', periodList), HWNAV0), 
          file.path(hwnavDir, 'allTokens_hwnav0.rds'))
  saveRDS(data.frame(Time = paste0('hwnav1.', periodList), HWNAV1), 
          file.path(hwnavDir, 'allTokens_hwnav1.rds'))
  
  cat("HWNAV computation complete!\n")
}

########################################################################
# 8. OW-HNPV (Overlap-Weighted Hierarchical Normalized Persistence Velocity)
########################################################################

computeOWHNPV <- function(networkDF, pdDir, owhnpvDir, m, n_sub, filtration = "sublevel") {
  
  extractOWHNPV <- function(D, alpha = 0, beta = 1.01) {
    if (nrow(D) == 0) return(rep(0, m))
    
    births <- D[, 1]
    deaths <- D[, 2]
    persistence <- deaths - births
    
    P_total <- sum(persistence)
    if (P_total == 0) return(rep(0, m))
    
    n_features <- nrow(D)
    
    delta_s <- (beta - alpha) / m
    s_points <- seq(alpha, beta, by = delta_s)
    
    owhnpv_vec <- numeric(m)
    
    for (j in 1:m) {
      sj <- s_points[j]
      sj1 <- s_points[j + 1]
      
      delta_t <- (sj1 - sj) / n_sub
      t_points <- seq(sj, sj1, by = delta_t)
      
      V_j_ow <- 0
      for (ell in 1:n_sub) {
        t_ell <- t_points[ell]
        t_ell1 <- t_points[ell + 1]
        delta_t_ell <- t_ell1 - t_ell
        
        # Compute overlap weights
        overlap_sum <- 0
        for (i in 1:n_features) {
          overlap_length <- max(0, min(deaths[i], t_ell1) - max(births[i], t_ell))
          overlap_sum <- overlap_sum + overlap_length
        }
        
        if (delta_t_ell > 0) {
          V_ell_ow <- overlap_sum / delta_t_ell
        } else {
          V_ell_ow <- 0
        }
        V_j_ow <- V_j_ow + V_ell_ow
      }
      
      owhnpv_vec[j] <- V_j_ow / (n_sub * P_total)
    }
    
    return(owhnpv_vec)
  }
  
  periodList <- unique(networkDF$time)
  nDays <- length(periodList)
  
  OWHNPV0 <- OWHNPV1 <- matrix(0, nrow = nDays, ncol = m)
  
  for (i in 1:nDays) {
    pd_file <- file.path(pdDir, paste0('PD_', filtration, '_', periodList[i], '.rds'))
    if (!file.exists(pd_file)) next
    
    PD <- readRDS(pd_file)
    OWHNPV0[i, ] <- extractOWHNPV(PD[PD[, 1] == 0, 2:3, drop = FALSE])
    OWHNPV1[i, ] <- extractOWHNPV(PD[PD[, 1] == 1, 2:3, drop = FALSE])
    
    if (i %% 30 == 0) cat(paste(i, 'out of', nDays, 'OW-HNPV processed\n'))
  }
  
  colnames(OWHNPV0) <- paste0('int', 1:m)
  colnames(OWHNPV1) <- paste0('int', 1:m)
  
  saveRDS(data.frame(Time = paste0('owhnpv0.', periodList), OWHNPV0), 
          file.path(owhnpvDir, 'allTokens_owhnpv0.rds'))
  saveRDS(data.frame(Time = paste0('owhnpv1.', periodList), OWHNPV1), 
          file.path(owhnpvDir, 'allTokens_owhnpv1.rds'))
  
  cat("OW-HNPV computation complete!\n")
}

########################################################################
# 9. ROLLING DEPTH
########################################################################

computeRollingDepth <- function(methodDir, methodName, depthDir, rollSize = 7) {
  
  files <- list.files(methodDir, pattern = "allTokens")
  
  for (f in files) {
    df <- readRDS(file.path(methodDir, f))
    
    # Normalize
    numeric_cols <- sapply(df, is.numeric)
    if (sum(numeric_cols) > 0) {
      df[, numeric_cols] <- df[, numeric_cols] / apply(df[, numeric_cols, drop = FALSE], 1, 
                                                        function(x) max(1, max(x, na.rm = TRUE)))
    }
    
    rollDepth_vec <- rep(NA, nrow(df))
    
    for (i in rollSize:nrow(df)) {
      sig_roll <- df[(i - rollSize + 1):i, -1, drop = FALSE]
      tryCatch({
        rollDepth_vec[i] <- MBD(df[(i - 1):i, -1, drop = FALSE], sig_roll, plotting = FALSE)$MBD[2]
      }, error = function(e) {
        rollDepth_vec[i] <- NA
      })
    }
    
    # Extract date from Time column (format: "method0.YYYY-MM-DD" or "method1.YYYY-MM-DD")
    time_str <- as.character(df$Time)
    date_part <- sub("^[^.]+\\.", "", time_str)  # Remove everything before and including first dot
    
    # Extract type (e.g., "vab0", "vab1")
    type_part <- sub("\\..*$", "", time_str)  # Remove everything after first dot
    
    res <- tibble(
      day = ymd(date_part),
      rollDepth = rollDepth_vec,
      type = type_part
    )
    
    output_file <- paste0("rd_", methodName, "_", f)
    saveRDS(res, file.path(depthDir, output_file))
  }
  
  cat(paste("Rolling depth for", methodName, "complete!\n"))
}

########################################################################
# 10. DATA MERGE
########################################################################

dataMerge <- function(graphDir, depthDir, priceDir, mergeDir, threshold = 0.05) {
  
  # Graph features
  ghFeature <- readRDS(file.path(graphDir, "graphFeature_allTokens.rds"))
  
  # Rolling depths - combine all
  depthFiles <- list.files(depthDir, pattern = "^rd_")
  
  if (length(depthFiles) == 0) {
    stop("No rolling depth files found in depthDir!")
  }
  
  rollDepth <- map_dfr(depthFiles, function(f) {
    readRDS(file.path(depthDir, f))
  })
  
  rollDepth <- rollDepth %>% 
    filter(!is.na(day)) %>%
    spread(type, rollDepth) %>% 
    arrange(day)
  
  # Price data
  priceFile <- list.files(priceDir, pattern = "Eth_price")
  if (length(priceFile) == 0) {
    stop("No price file found in priceDir!")
  }
  
  priceEth <- read_csv(file.path(priceDir, priceFile[1]), 
                       col_types = cols(Date = col_date(format = "%Y-%m-%d")),
                       show_col_types = FALSE)
  priceEth <- rename(priceEth, Open = `24h Open (USD)`)
  
  # Process price
  price <- priceEth %>% dplyr::select(Date, Open)
  priceVec <- price$Open
  priceDiff <- diff(priceVec)
  priceReturn <- priceDiff / priceVec[1:length(priceDiff)]
  price <- price %>% mutate(
    openNorm = Open / max(Open),
    priceReturn = c(NA, priceReturn)
  )
  
  # Flag periods
  flagDaysFun <- function(inputVector, ndays) {
    absReturn <- abs(inputVector)
    l <- length(inputVector)
    flagDays <- vector(length = l)
    for (i in 1:(l - 1)) {
      maxChange <- max(absReturn[(i + 1):min(i + ndays, l)], na.rm = TRUE)
      flagDays[i] <- ifelse(maxChange >= threshold, TRUE, FALSE)
    }
    flagDays[l] <- FALSE
    return(flagDays)
  }
  
  for (h in 1:7) {
    price <- price %>% mutate(!!paste0("flag", h) := flagDaysFun(priceReturn, h))
  }
  
  # Merge
  df <- inner_join(ghFeature, rollDepth, by = "day")
  df <- inner_join(df, price, by = c("day" = "Date"))
  
  saveRDS(df, file.path(mergeDir, paste0("df_merged_abs", 100 * threshold, ".rds")))
  
  cat(paste("Data merged! Shape:", nrow(df), "rows,", ncol(df), "columns\n"))
  cat(paste("Columns:", paste(names(df), collapse = ", "), "\n"))
  
  return(df)
}

########################################################################
# 11. RANDOM FOREST MODEL
########################################################################

fitRFModel <- function(mergeDir, modelDir, methodName, threshold = 0.05, repNum = 10) {
  
  df <- readRDS(file.path(mergeDir, paste0("df_merged_abs", 100 * threshold, ".rds")))
  
  flags <- df %>% dplyr::select(contains("flag")) %>% colnames()
  df <- df %>% mutate(across(starts_with("flag"), ~factor(., levels = c("TRUE", "FALSE"))))
  df$clusterCoef <- replace(df$clusterCoef, is.nan(df$clusterCoef), 0)
  
  # Train/test split
  dateRange <- range(df$day)
  spDate <- dateRange[1] + diff(dateRange) / 3 * 2
  dftrain <- df[df$day < spDate, ]
  dftest <- df[df$day >= spDate, ]
  
  # Define formulas based on method
  baseVars <- "vertexNum + edgeNum + clusterCoef + openNorm"
  
  # Check which columns exist (no trailing underscore)
  dim0_col <- paste0(methodName, "0")
  dim1_col <- paste0(methodName, "1")
  
  cat(paste("Looking for columns:", dim0_col, "and", dim1_col, "\n"))
  cat(paste("Available columns:", paste(names(df), collapse = ", "), "\n"))
  
  inputfmls <- list(
    M1 = paste0("~ ", baseVars)
  )
  
  if (dim0_col %in% names(df)) {
    inputfmls$M2 <- paste0("~ ", baseVars, " + ", dim0_col)
    cat(paste("M2 formula will use:", dim0_col, "\n"))
  } else {
    cat(paste("WARNING:", dim0_col, "not found in data!\n"))
  }
  
  if (dim0_col %in% names(df) && dim1_col %in% names(df)) {
    inputfmls$M3 <- paste0("~ ", baseVars, " + ", dim0_col, " + ", dim1_col)
    cat(paste("M3 formula will use:", dim0_col, "and", dim1_col, "\n"))
  } else if (!(dim1_col %in% names(df))) {
    cat(paste("WARNING:", dim1_col, "not found in data!\n"))
  }
  
  cat(paste("Formulas to fit:", paste(names(inputfmls), collapse = ", "), "\n"))
  
  # Run RF
  results <- map_dfr(flags, function(flagi) {
    map_dfr(1:repNum, function(repID) {
      map2_dfr(inputfmls, names(inputfmls), function(fmlstr, fmlname) {
        
        tryCatch({
          rf <- randomForest::randomForest(
            as.formula(paste0(flagi, fmlstr)),
            dftrain, importance = TRUE, na.action = na.omit
          )
          
          pred <- predict(rf, dftest)
          prob <- predict(rf, dftest, type = 'prob')[, "TRUE"]
          
          cf <- caret::confusionMatrix(pred, dftest %>% pull(flagi), positive = "TRUE")
          
          auc_val <- pROC::auc(
            as.numeric(dftest[[flagi]] == "TRUE"),
            prob, quiet = TRUE
          ) %>% as.numeric()
          
          tibble(
            horizon = flagi,
            modelType = fmlname,
            Accuracy = cf$overall["Accuracy"],
            Sensitivity = cf$byClass["Sensitivity"],
            Precision = cf$byClass["Precision"],
            AUC = auc_val,
            repID = repID
          )
        }, error = function(e) {
          cat(paste("Error in RF for", flagi, fmlname, ":", e$message, "\n"))
          return(NULL)
        })
      })
    })
  })
  
  saveRDS(results, file.path(modelDir, paste0("rf_results_", methodName, ".rds")))
  cat(paste("RF results saved for", methodName, "\n"))
  return(results)
}

########################################################################
# 12. GATHER AND SUMMARIZE RESULTS
########################################################################

summarizeResults <- function(modelDir, methodName) {
  
  results <- readRDS(file.path(modelDir, paste0("rf_results_", methodName, ".rds")))
  
  summary <- results %>%
    group_by(horizon, modelType) %>%
    summarise(
      Accuracy = mean(Accuracy, na.rm = TRUE),
      Sensitivity = mean(Sensitivity, na.rm = TRUE),
      Precision = mean(Precision, na.rm = TRUE),
      AUC = mean(AUC, na.rm = TRUE),
      .groups = 'drop'
    )
  
  return(summary)
}

########################################################################
# 13. COMPUTE AUC GAIN
########################################################################

computeAUCGain <- function(summary_df) {
  
  wide <- summary_df %>%
    dplyr::select(horizon, modelType, AUC) %>%
    pivot_wider(names_from = modelType, values_from = AUC)
  
  # Initialize gain columns
  wide$M2_gain <- NA

  wide$M3_gain <- NA
  
  if ("M2" %in% names(wide) && "M1" %in% names(wide)) {
    wide$M2_gain <- (wide$M2 - wide$M1) / wide$M1 * 100
  }
  if ("M3" %in% names(wide) && "M1" %in% names(wide)) {
    wide$M3_gain <- (wide$M3 - wide$M1) / wide$M1 * 100
  }
  
  wide$horizon_num <- as.numeric(gsub("flag", "", wide$horizon))
  
  return(wide)
}
