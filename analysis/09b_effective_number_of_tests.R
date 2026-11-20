# how many independent protein tests we actually have
# needs individual-level UK Biobank data, so it won't run from this repo
# gives 1,821 of 2,922 proteins, which is hard-coded in 09_proteome_mr.R

library(tidyverse)

# residualise each protein on the covariates
COVARIATES <- paste(
    "age + sex + age_squared + age_x_sex + age_squared_x_sex",
    "assessment_centre",
    "genotype_batch",
    "olink_batch",
    paste(sprintf("genetic_PC_%d", 1:20), collapse = " + "),
    sep = " + "
)

residualise <- function(protein, d) {
    f <- as.formula(paste(protein, "~", COVARIATES))
    keep <- d %>% select(ID, all_of(all.vars(f))) %>% drop_na()
    out <- data.frame(keep$ID, residuals(lm(f, data = keep)))
    setNames(out, c("ID", protein))
}

eur <- ukb_ppp %>% filter(genetic_ancestry_superpopulation == "EUR")
adjusted <- reduce(lapply(protein_cols, residualise, d = eur), full_join, by = "ID")

# simpleM (Gao) on the protein correlation matrix
eigenvalues <- eigen(
    cor(column_to_rownames(adjusted, "ID"), use = "pairwise.complete.obs")
)$values

gao <- function(ev, C) which(cumsum(sort(ev, decreasing = TRUE)) / sum(ev) > C)[1]

gao(eigenvalues, 0.95)   # 1821
