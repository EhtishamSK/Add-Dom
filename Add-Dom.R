# -------------------------------
# GBLUP with Additive + Dominance Marker Effects
# -------------------------------
# Author: Ehtisham Khokhar
# Email: ehtishamshakeel@gmail.com 
# Description: This script demonstrates how to implement GBLUP (Genomic Best Linear Unbiased Prediction)
#              using additive and dominance effects with the BGLR package in R. 
#              Input files: geno.csv (genotype matrix), pheno.txt (phenotype data).
#              Output files: marker_effects.csv (SNP effects), genetic_values.csv (GEBVs per individual).
# Note: This script is prepared for GitHub repository documentation purposes.
# -------------------------------

# ---- 0) Install/load required package ----
# The BGLR package is used to fit Bayesian genomic prediction models
# install.packages("BGLR")   # Uncomment and run if BGLR is not installed
library(BGLR)

# ---- 1) Set working directory ----
# Update path to the folder containing input files
setwd("C:/Users/ehtis/Downloads/simplemate")

# ---- 2) Load genotype data ----
# Input: geno.csv (rows = individuals, columns = SNPs coded as 0/1/2)
geno <- read.csv("geno.csv", row.names = 1)
# Preview first 5 rows and columns
geno[1:5, 1:5]

# ---- 3) Load phenotype data ----
# Input: pheno.txt (tab-delimited file with trait values)
pheno <- read.delim("pheno.txt", header = TRUE)
# Assign row names based on the first column (Taxa/individual ID)
rownames(pheno) <- pheno[[1]]  
# Inspect structure and first 5 entries
str(pheno)
head(pheno, 5)

# Extract phenotype vector (example: AUDPC trait)
y <- pheno$AUDPC

# ---- 4) Preprocess genotype data ----
# Step 1: Remove SNPs with all missing values
geno <- geno[, colSums(is.na(geno)) < nrow(geno)]

# Step 2: Impute missing genotypes with column mean
geno_imp <- geno
for(j in seq_len(ncol(geno_imp))){
  m <- mean(geno_imp[, j], na.rm = TRUE)
  geno_imp[is.na(geno_imp[, j]), j] <- m
}

# Step 3: Remove monomorphic SNPs (no variation across individuals)
is_monomorphic <- apply(geno_imp, 2, function(col) length(unique(col)) == 1)
geno_imp <- geno_imp[, !is_monomorphic]

# ---- 5) Build additive and dominance design matrices ----
# Calculate allele frequencies
p <- colMeans(geno_imp) / 2
q <- 1 - p

# Additive design matrix (centered by 2p)
Z_add <- scale(geno_imp, center = 2*p, scale = FALSE)

# Dominance design matrix using Vitezica et al. 2013 coding scheme
make_dom <- function(M){
  p <- colMeans(M)/2
  q <- 1 - p
  Z <- matrix(0, nrow(M), ncol(M), dimnames = dimnames(M))
  for(j in seq_len(ncol(M))){
    xj <- M[, j]
    pj <- p[j]; qj <- q[j]
    Z[, j] <- ifelse(xj == 0, -2*qj^2,
                     ifelse(xj == 1, 2*pj*qj,
                            ifelse(xj == 2, -2*pj^2, 0)))  # Missing values → 0
  }
  Z
}
Z_dom <- make_dom(geno_imp)

# Double-check for missing values
cat("NA in Z_add:", sum(is.na(Z_add)), "\n")
cat("NA in Z_dom:", sum(is.na(Z_dom)), "\n")

# ---- 6) Fit BGLR model with additive + dominance effects ----
# ETA list defines model components:
# A = additive effects, D = dominance effects
ETA <- list(
  A = list(X = Z_add, model = "BRR"),  # Bayesian Ridge Regression for additive
  D = list(X = Z_dom, model = "BRR")   # Bayesian Ridge Regression for dominance
)

# Run BGLR model
fit <- BGLR(y = y, ETA = ETA,
            nIter = 12000, burnIn = 2000,
            verbose = TRUE)

# ---- 7) Extract marker effects ----
# α = additive SNP effects, δ = dominance SNP effects
alpha <- fit$ETA$A$b   
delta <- fit$ETA$D$b   

# ---- 8) Compute genetic values per individual ----
# GEBV (additive genetic value), dominance deviation, and total genotypic value
gebv_add <- as.vector(Z_add %*% alpha)  # additive GEBV
dom_dev  <- as.vector(Z_dom %*% delta)  # dominance deviation
g_total  <- gebv_add + dom_dev          # total genetic value

# ---- 9) Save results ----
# (a) Marker effects
write.csv(data.frame(SNP = colnames(geno_imp),
                     Add_Effect = alpha,
                     Dom_Effect = delta),
          "marker_effects.csv", row.names = FALSE)

# (b) Individual-level genetic values
write.csv(data.frame(Taxa = rownames(geno_imp),
                     Add_GEBV = gebv_add,
                     Dom_Dev = dom_dev,
                     Total_GValue = g_total),
          "genetic_values.csv", row.names = FALSE)
