remove_outliers <- function(rawdata) {
  
  N <- ncol(rawdata)
  X <- rawdata[,2:N]
  
  median_X <- apply(X, 2, median, na.rm = TRUE)
  
  # Repeat median of each series over all data points in the series
  median_X_mat <- matrix(rep(median_X, nrow(X)), nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  
  # Calculate quartiles
  Q <- apply(X, 2, stats::quantile, probs = c(0.25, 0.75), na.rm = TRUE)
  
  # Calculate interquartile range (IQR) of each series
  IQR <- Q[2,] - Q[1,]
  
  # Repeat IQR of each series over all data points in the series
  IQR_mat <- matrix(rep(IQR, nrow(X)), nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  
  # Determine outliers
  Z <- abs(X - median_X_mat)
  outlier <- (Z > (10 * IQR_mat))
  
  # Replace outliers with NaN
  Y <- X
  Y[outlier] <- NA
  
  # Cleaned data
  outdata <- rawdata
  outdata[, 2:N] <- Y
  
  return(outdata)
  
}