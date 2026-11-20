# same idea but with NLRP3 expression as the exposure
# uses expression-only cis-eQTLs, and the eight main instruments with their eQTL effects

suppressPackageStartupMessages(library(MendelianRandomization))
source(here::here("config.R"))
source(here::here("helpers.R"))

out_dir <- step_dir("11b_expression_only")

sel <- select_expression_instruments(READOUTS$NLRP3_expression, READOUTS$Neutrophil_count,
                                     CHR, INSTRUMENT_START, INSTRUMENT_END)
score_snps <- load_instruments()$SNP

raw <- lapply(READOUTS, read_region, chr = CHR, start = LOCUS_START, end = LOCUS_END)
readouts <- lapply(names(READOUTS), function(k) harmonise(raw[[k]], READOUTS[[k]]))
names(readouts) <- names(READOUTS)

# per SD lower expression
expression_exposure <- function(snps) {
    readouts$NLRP3_expression %>%
        filter(SNPid %in% snps) %>%
        transmute(SNP = SNPid, chr, pos_hg38 = pos, A1, A2, beta_exposure = -beta,
                  se_exposure = se) %>%
        arrange(pos_hg38)
}
SETS <- list(`Expression-only instrument` = expression_exposure(sel$instruments),
             `Activity score variants` = expression_exposure(score_snps))
ld_full <- interval_ld_matrix(unique(c(SETS[[1]]$SNP, SETS[[2]]$SNP)))

# joint F and R2
strength <- lapply(names(SETS), function(s) {
    e <- SETS[[s]]
    f <- joint_F_R2(e$beta_exposure, e$se_exposure, READOUTS$NLRP3_expression$n,
                    ld_full[e$SNP, e$SNP])
    tibble(set = s, n_snps = nrow(e), n_samples = READOUTS$NLRP3_expression$n,
           chi2 = f$chi2, F_joint = f$F_joint, F_mean = f$F_mean, F_min = f$F_min,
           R2_joint = f$R2_joint, R2_rho = f$R2_rho, R2_adj = f$R2_adj)
}) %>% bind_rows()


# CRP, GlycA, neutrophils
readout_mr <- function(e) {
    lapply(c("CRP", "GlycA", "Neutrophil_count"), function(k) {
        e %>%
            select(SNP, beta_exposure, se_exposure) %>%
            left_join(transmute(readouts[[k]], SNP = SNPid, beta_outcome = beta, se_outcome = se),
                      by = "SNP") %>%
            mutate(outcome = READOUTS[[k]]$label) %>%
            run_mr(ld_full)
    }) %>% bind_rows()
}


# IL1B, IL18, IL6
assays <- tibble(protein_id = readLines(ppp_manifest)) %>%
    mutate(gene_name = sub("_.*$", "", protein_id)) %>%
    filter(gene_name %in% c("IL1B", "IL18", "IL6"))
cytokine_mr <- function(e) {
    best <- assays %>%
        mutate(mr = map(protein_id, function(id) {
            mr_ppp_assay(file.path(ppp_dir, paste0(id, ".bgz")), e, ld_full)$res
        })) %>%
        unnest(mr) %>%
        group_by(gene_name) %>%
        slice_min(ivw_pval, n = 1, with_ties = FALSE) %>%
        ungroup()
    lapply(c("IL1B", "IL18", "IL6"), function(k) {
        r <- best[best$gene_name == k, ]
        tibble(outcome = k, protein_id = r$protein_id, nsnp = r$n_snps,
               method = c("IVW", "Weighted median"),
               estimate = c(r$ivw_beta, r$median_beta), se = c(r$ivw_se, r$median_se),
               p = c(r$ivw_pval, r$median_pval),
               ci_lower = estimate - 1.96 * se, ci_upper = estimate + 1.96 * se,
               egger_intercept = r$egger_intercept, egger_intercept_p = r$egger_intercept_pval,
               het_q = c(r$ivw_het_q, NA), het_p = c(r$ivw_het_p, NA))
    }) %>% bind_rows()
}


# CAD
cohorts <- CAD_STUDIES[c("aragam", "mvp", "allofus", "finngen")]
cad <- lapply(cohorts, read_cad)
names(cad) <- sapply(cohorts, `[[`, "label")
cad_mr <- function(e) {
    studies <- lapply(cad, cad_at, ex = transmute(e, SNP, pos_hg38, A1, A2,
                                                  bx = beta_exposure, bxse = se_exposure))
    fit <- function(d, label) {
        d %>%
            transmute(SNP, beta_exposure = bx, se_exposure = bxse,
                      beta_outcome = by, se_outcome = byse, outcome = "CAD") %>%
            run_mr(ld_full) %>%
            mutate(study = label)
    }
    bind_rows(lapply(names(studies), function(n) fit(studies[[n]], n)),
              fit(meta_snp(bind_rows(studies)), "Meta-analysis"))
}

results <- lapply(names(SETS), function(s) {
    bind_rows(readout_mr(SETS[[s]]), cytokine_mr(SETS[[s]]), cad_mr(SETS[[s]])) %>%
        mutate(set = s)
}) %>%
    bind_rows() %>%
    mutate(method = recode(method, "IVW (LD-corrected)" = "IVW")) %>%
    select(set, outcome, study, protein_id, method, n_snps = nsnp, estimate, se,
           ci_lower, ci_upper, p, egger_intercept, egger_intercept_p, het_q, het_p)


# per-variant table
# flip to the expression-lowering allele
inst <- SETS[[1]]
flip_snp <- setNames(inst$beta_exposure < 0, inst$SNP)
rsids <- lookup_rsids(unique(c(inst$SNP, score_snps)), raw[c("Neutrophil_count", "CRP", "GlycA")])

# is each variant an instrument, a proxy of one, or neither
r2_score <- ld_full[inst$SNP, score_snps]^2
tag <- score_snps[apply(r2_score, 1, which.max)]
status <- tibble(SNP = inst$SNP, r2_score = apply(r2_score, 1, max),
                 r2_score_rsid = rsids$rsid[match(tag, rsids$SNP)]) %>%
    mutate(score_status = case_when(SNP %in% score_snps ~ "Activity score instrument",
                                    r2_score >= PROXY_R2 ~ "In LD with an activity score signal",
                                    TRUE ~ "Neither"))

per_variant <- expand_grid(SNP = inst$SNP, trait = names(READOUTS)) %>%
    left_join(bind_rows(readouts, .id = "trait") %>% select(trait, SNP = SNPid, beta, se, p),
              by = c("SNP", "trait")) %>%
    left_join(select(inst, SNP, chr, pos_hg38, A1, A2), by = "SNP") %>%
    left_join(transmute(readouts$Neutrophil_count, SNP = SNPid, eaf), by = "SNP") %>%
    left_join(rsids, by = "SNP") %>%
    left_join(status, by = "SNP") %>%
    mutate(flipped = flip_snp[SNP]) %>%
    transmute(SNP, rsid, chr, pos_hg38,
              effect_allele = ifelse(flipped, A2, A1), other_allele = ifelse(flipped, A1, A2),
              eaf = ifelse(flipped, 1 - eaf, eaf),
              trait, beta = ifelse(flipped, -beta, beta), se, p,
              score_status, r2_score, r2_score_rsid)

print(as.data.frame(results %>% transmute(
    set, outcome, study, method, n_snps,
    estimate = ifelse(outcome == "CAD", sprintf("OR %.3f (%.3f, %.3f)", exp(estimate),
                                                exp(ci_lower), exp(ci_upper)),
                      sprintf("%.3f (%.3f, %.3f)", estimate, ci_lower, ci_upper)),
    p = signif(p, 3), het_p = signif(het_p, 3))), row.names = FALSE)

write_tsv(results, file.path(out_dir, "expression_only_mr.tsv"))
write_tsv(strength, file.path(out_dir, "expression_only_strength.tsv"))
write_tsv(per_variant, file.path(out_dir, "expression_only_instrument_table.tsv"))
write_tsv(sel$specificity, file.path(out_dir, "expression_only_specificity.tsv"))
