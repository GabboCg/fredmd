prepare_missing <- function(rawdata, tcode, vardate = "sasdate") {
  
  # ===========================================================================
  # DESCRIPTION:
  #    This function transforms raw data based on each series' transformation
  #    code.
  # 
  # ---------------------------------------------------------------------------
  # INPUT:
  #     rawdata = raw data
  #     tcode = transformation codes for each series
  # 
  # OUTPUT:
  #     y_t = transformed data
  # 
  # ---------------------------------------------------------------------------
  # SUBFUNCTION:
  #     transxf: transforms a single series as specified by a given 
  #              transfromation code
  # 
  # ===========================================================================
  
  rawdata_date_id <- tibble::column_to_rownames(rawdata, var = vardate) 
  variables <- colnames(rawdata_date_id) 
  
  transformed_data <- list()
  
  for (i in seq_along(variables)) {
    
    x <- rawdata_date_id %>% select(all_of(variables[i]))
    transformed_data[[i]] <- transxf(x, as.numeric(tcode[i]))
    
  }
  
  y_t_all <- bind_cols(transformed_data) %>% 
    dplyr::mutate(yyyymm = rownames(.)) %>% 
    dplyr::as_tibble() %>% 
    dplyr::relocate(yyyymm) 
  
  return(y_t_all)  
  
}

transxf <- function(x, tcode = c("1", "2", "3", "4", "5", "6", "7")) {
  
  # ===========================================================================
  # DESCRIPTION:
  #  This function transforms a SINGLE SERIES (in a column vector) as specified
  #  by a given transformation code.
  # 
  # ---------------------------------------------------------------------------
  # INPUT:
  #     x = series (in a column vector) to be transformed
  #     tcode = transformation code (1-7)
  # 
  # OUTPUT:
  #     y_t = transformed series (as a column vector)
  # ---------------------------------------------------------------------------

  y_t <- switch(tcode,
                "1" = x,
                "2" = x - lag(x, 1),
                "3" = (x - lag(x, 1)) - (lag(x, 1) - lag(x, 2)),
                "4" = log(x),
                "5" = log(x) - lag(log(x), 1),
                "6" = (log(x) - lag(log(x), 1)) - (lag(log(x), 1) - lag(log(x), 2)), 
                "7" = (x / lag(x, 1) - 1)  - (lag(x, 1) / lag(x, 2) - 1))
  
  return(y_t)
  
}
