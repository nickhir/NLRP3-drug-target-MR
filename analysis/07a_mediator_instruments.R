# get instruments for the mediators (SBP, ApoB, T2D)

source(here::here("config.R"))
source(here::here("helpers.R"))

out_dir <- step_dir("07a_mediator_instruments")
MEDS <- names(MEDIATORS)

# plink can't break ties at P = 0 (top ApoB hits), so rank on -log10(p)
clump_ranked <- function(snps, nlog10, ...) {
    p <- rank(-nlog10, ties.method = "first") / (length(nlog10) + 1) * MED_CLUMP_P
    clumped <- ld_clump_local(data.frame(SNP = to_panel_id(snps), p = p),
                              ld_panel, MED_CLUMP_R2, MED_CLUMP_KB, ...)
    from_panel_id(clumped$ID)
}


# instruments per mediator
instruments <- lapply(MEDS, function(k) {
    cfg <- MEDIATORS[[k]]
    sig <- read_genome(cfg) %>%
        filter(as.numeric(p) < MED_CLUMP_P) %>%
        harmonise(cfg) %>%
        arrange(desc(nlog10))
    keep <- clump_ranked(sig$SNPid, sig$nlog10)
    sig %>% filter(SNPid %in% keep) %>% mutate(mediator = k, .before = 1)
}) %>% bind_rows()


# all instruments in all mediators and CAD
union_snps <- distinct(instruments, SNPid, chrom = chr, pos)

med_at_union <- lapply(MEDS, function(k) {
    lookup_at(MEDIATORS[[k]], union_snps) %>% transmute(SNPid, mediator = k, beta, se, eaf, n)
}) %>% bind_rows()

# SBP from mmHg to SD
sd_scales <- sapply(MEDS, function(k) {
    if (!isTRUE(MEDIATORS[[k]]$standardise_sd)) return(1)
    d <- filter(med_at_union, mediator == k)
    sd_from_se(d$se, d$eaf, d$n)
})
med_at_union <- mutate(med_at_union, beta_sd = beta / sd_scales[mediator],
                       se_sd = se / sd_scales[mediator])

# All of Us is an extract, missing variants are meta-analysed without it
cad <- lapply(names(CAD_STUDIES), function(k) {
    lookup_at(CAD_STUDIES[[k]], union_snps) %>% transmute(SNPid, beta, se, study = k)
}) %>% bind_rows()

# fixed effect meta-analysis
cad_meta <- cad %>%
    group_by(SNPid) %>%
    summarise(
        beta_cad = sum(beta / se^2) / sum(1 / se^2),
        se_cad = sqrt(1 / sum(1 / se^2)),
        n_studies = n(),
        Q = sum((beta - sum(beta / se^2) / sum(1 / se^2))^2 / se^2),
        .groups = "drop"
    ) %>%
    mutate(Q_df = n_studies - 1,
           Q_p = ifelse(Q_df > 0, pchisq(Q, Q_df, lower.tail = FALSE), NA_real_),
           I2 = ifelse(Q_df > 0, pmax(0, (Q - Q_df) / Q) * 100, NA_real_))


# MVMR design
design <- med_at_union %>%
    select(SNPid, mediator, beta_sd, se_sd) %>%
    pivot_wider(names_from = mediator, values_from = c(beta_sd, se_sd)) %>%
    inner_join(select(cad_meta, SNPid, beta_cad, se_cad, n_studies, I2), by = "SNPid") %>%
    drop_na(starts_with("beta_sd_"), starts_with("se_sd_"))

# re-clump the pooled instruments for each nested model (SBP, +ApoB, +T2D)
mvmr_instrument_sets <- lapply(seq_along(MEDS), function(i) {
    ins <- filter(instruments, mediator %in% MEDS[1:i], SNPid %in% design$SNPid)
    priority <- tapply(ins$nlog10, ins$SNPid, max)
    clump_ranked(names(priority), priority, extra = c("--threads", 1, "--memory", 2048))
})
names(mvmr_instrument_sets) <- sapply(seq_along(MEDS), function(i) paste(MEDS[1:i], collapse = "+"))


saveRDS(mvmr_instrument_sets, file.path(out_dir, "mvmr_instrument_sets.rds"))
fwrite(instruments, file.path(out_dir, "mediator_instruments.tsv"), sep = "\t")
fwrite(med_at_union, file.path(out_dir, "mediators_at_union.tsv"), sep = "\t")
fwrite(design, file.path(out_dir, "mvmr_design.tsv"), sep = "\t")
fwrite(cad_meta, file.path(out_dir, "cad_meta_at_instruments.tsv"), sep = "\t")
saveRDS(list(sd_scales = sd_scales,
             n_instruments = table(instruments$mediator),
             clump = c(p = MED_CLUMP_P, r2 = MED_CLUMP_R2, kb = MED_CLUMP_KB),
             ld_panel = ld_panel),
        file.path(out_dir, "selection_summary.rds"))
