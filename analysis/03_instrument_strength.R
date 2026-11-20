# how strong are the instruments

source(here::here("config.R"))
source(here::here("helpers.R"))

out_dir <- step_dir("03_instrument_strength")

aligned <- fread(file.path(results_dir, "02_instrument_table", "instrument_table_aligned.tsv"),
                 data.table = FALSE)
K <- nrow(aligned)

ld <- interval_ld_matrix(aligned$SNP)
# joint F and R2
strength <- c(
    lapply(names(READOUTS), function(k) list(
        label = k, n = READOUTS[[k]]$n,
        beta = aligned[[paste0("beta_", k)]], se = aligned[[paste0("se_", k)]])),
    list(list(label = "cis_NLRP3_activity_score", n = NA_real_,
              beta = aligned$beta_exposure, se = aligned$se_exposure))
) %>% lapply(function(x) c(x, joint_F_R2(x$beta, x$se, x$n, ld)))

joint_tbl <- lapply(strength, function(s) data.frame(
    trait = s$label, n_snps = K, n_samples = s$n,
    chi2 = s$chi2, F_joint = s$F_joint, F_mean = s$F_mean, F_min = s$F_min,
    R2_joint = s$R2_joint, R2_rho = s$R2_rho, R2_adj = s$R2_adj
)) %>% bind_rows()
print(joint_tbl, row.names = FALSE, digits = 4)

per_snp_tbl <- lapply(strength, function(s) data.frame(
    SNP = aligned$SNP, trait = s$label, beta = s$beta, se = s$se, z = s$z,
    F = s$F_per_snp, R2 = s$R2_per_snp
)) %>% bind_rows()

# effective N from the SEs, should be close to the reported N
maf <- pmin(aligned$eaf_exposure, 1 - aligned$eaf_exposure)
eff_n <- lapply(names(READOUTS), function(k) {
    ne <- median(1 / (2 * maf * (1 - maf) * aligned[[paste0("se_", k)]]^2))
    data.frame(trait = k, n_stated = READOUTS[[k]]$n, n_effective = ne,
               ratio = ne / READOUTS[[k]]$n)
}) %>% bind_rows()

fwrite(joint_tbl, file.path(out_dir, "instrument_strength_joint.tsv"), sep = "\t")
fwrite(per_snp_tbl, file.path(out_dir, "instrument_strength_per_snp.tsv"), sep = "\t")
fwrite(eff_n, file.path(out_dir, "effective_n_check.tsv"), sep = "\t")
saveRDS(ld, file.path(out_dir, "instrument_ld_matrix.rds"))
