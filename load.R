#!/usr/bin/env Rscript
# ======================================================== #
#
#                   Data Cleaning FRED-MD
#
#                 Gabriel E. Cabrera-Guzmán
#                The University of Manchester
#
#                        Spring, 2026
#
#                https://gcabrerag.rbind.io
#
# ------------------------------ #
# email: gabriel.cabreraguzman@postgrad.manchester.ac.uk
# ======================================================== #

# Load packages
library(tidyverse)
library(lubridate)
library(Rcpp)

# Load auxiliary functions
source("R/prepare-missing.R")
source("R/remove-outliers.R")

Rcpp::sourceCpp("src/fred_factors.cpp")

# ==========================================
#               FRED MD Cleaning
# ------------------------------------------

# Read raw FRED MD dataset
raw_fred_md <- read_csv("data/FRED_MD_062024.csv")

# Extract columns names and transformation code
col_names <- colnames(raw_fred_md)[-1]
tcode <- raw_fred_md[1, -1]
raw_fred_md <- raw_fred_md |> slice(-1)

# Transform each series based on his code
transformed_fred_md <- prepare_missing(raw_fred_md, tcode, vardate = "sasdate")

# Change column date format
transformed_fred_md <- transformed_fred_md |>
  mutate(
    yyyymm = as.Date(yyyymm, format = "%m/%d/%Y")
  ) |>
  slice(-1:-2) # Remove the first two months

# Remove missing values
transformed_fred_md <- remove_outliers(rawdata = transformed_fred_md)

# Extract current date
current_date_fred_md <- transformed_fred_md$yyyymm

# Type of transformation performed on each series before factors are estimated
demean <- 2

# Information criterion used to select number of factors; for more details,
jj <- 2

# Maximum number of factors to be estimated; if set to 99, the number of
# factors selected is forced to equal 8
kmax <- 8

# Estimate factors — strip date column before calling C++
x_mat <- as.matrix(transformed_fred_md[, -1])   # T0 × N numeric matrix
result <- factors_em_cpp(x_mat, kmax, jj, demean)

# Add colnames
colnames(result[[5]]) <- colnames(x_mat)

fred_md_balanced <- result[[5]] |>
  as_tibble() |>
  mutate(date = current_date_fred_md) |>
  relocate(date, .before = 1) |>
  janitor::clean_names() |>
  drop_na()
