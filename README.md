## GBLUP with Additive and Dominance Effects

This R script implements **Genomic Best Linear Unbiased Prediction (GBLUP)** with both additive and dominance marker effects using the [BGLR](https://cran.r-project.org/package=BGLR) package.  
It demonstrates a full workflow from raw genotype/phenotype data to genomic estimated breeding values (GEBVs) and SNP effect estimates.  

---

### Workflow

1. **Load Libraries & Data**  
   - Load the **BGLR** package.  
   - Import genotype data (`geno.csv`) with SNPs coded as `0/1/2`.  
   - Import phenotype data (`pheno.txt`) containing trait measurements (e.g., `AUDPC`).  

2. **Preprocessing**  
   - Remove SNPs with all missing values.  
   - Impute missing genotypes using column means.  
   - Filter out monomorphic markers (no genetic variation).  

3. **Design Matrices**  
   - Compute allele frequencies (`p`, `q`).  
   - Build **additive design matrix (Z_add)** centered on `2p`.  
   - Build **dominance design matrix (Z_dom)** using the coding scheme of *Vitezica et al. 2013*.  

4. **Model Fitting (BGLR)**  
   - Specify model components (`ETA`) for additive and dominance effects.  
   - Fit the **Bayesian Ridge Regression (BRR)** model with MCMC iterations (`nIter = 12000`, `burnIn = 2000`).  

5. **Extract Marker Effects**  
   - Obtain **additive effects (α)** and **dominance effects (δ)** for each SNP.  

6. **Compute Genetic Values**  
   - Calculate additive GEBVs, dominance deviations, and total genotypic values for each individual.  

7. **Save Results**  
   - **marker_effects.csv** → additive and dominance SNP effects.  
   - **genetic_values.csv** → individual-level predictions (Add_GEBV, Dom_Dev, Total_GValue).  

---

### Outputs

- **Marker Effects**  
  Table of SNP effects with additive and dominance components.  
  | SNP | Add_Effect | Dom_Effect |

- **Genetic Values**  
  Genomic predictions for each individual.  
  | Taxa | Add_GEBV | Dom_Dev | Total_GValue |

---

### Summary

This workflow provides a comprehensive pipeline for:  
- Cleaning and preparing genotype and phenotype data  
- Partitioning genetic variance into additive and dominance components  
- Predicting total genetic values for breeding populations  
- Exporting results for downstream analyses (e.g., genomic selection, GWAS follow-up, or breeding decisions)  

By combining **additive and dominance models**, breeders can capture a more realistic representation of genetic architecture and improve prediction accuracy for complex traits.  

---

### Author
- **Ehtisham Khokhar**  
- 📧 ehtishamkhokhar@gmail.com  

---

### References
- Pérez, P., & de los Campos, G. (2014). BGLR: A statistical package for whole genome regression and prediction. *Genetics*, 198(2), 483–495.  
- Vitezica, Z. G., Varona, L., & Legarra, A. (2013). On the additive and dominant variance and covariance of individuals within the genomic selection scope. *Genetics*, 195(4), 1223–1230.  
