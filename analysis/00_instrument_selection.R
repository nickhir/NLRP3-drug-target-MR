# pick the instruments and build the NLRP3 activity score

source(here::here("config.R"))
source(here::here("helpers.R"))

out_dir <- step_dir("00_instrument_selection")

sel <- select_instruments(READOUTS, CHR, INSTRUMENT_START, INSTRUMENT_END)

# eaf from the neutrophil GWAS, eaf_panel from INTERVAL
out <- sel$score %>%
    transmute(SNP, chr = CHR, pos_hg38 = as.integer(snp_field(SNP, 2)),
              A1 = snp_field(SNP, 3), A2 = snp_field(SNP, 4),
              beta_exposure = beta, se_exposure = se) %>%
    left_join(select(sel$stats$neutro, SNP, eaf = A1_freq), by = "SNP") %>%
    left_join(sel$panel_eaf, by = "SNP") %>%
    mutate(r2_threshold = 0.1, pc1_var_explained = sel$pc1_var,
           scaling_constant_k = sel$k, score_sign_flipped = sel$flipped) %>%
    select(SNP, chr, pos_hg38, A1, A2, eaf, beta_exposure, se_exposure, eaf_panel,
           r2_threshold, pc1_var_explained, scaling_constant_k, score_sign_flipped) %>%
    arrange(pos_hg38)

write_tsv(out, file.path(out_dir, "nlrp3_instruments.tsv"))
write_tsv(sel$proxies, file.path(out_dir, "nlrp3_proxy_replacements.tsv"))
write_tsv(sel$specificity, file.path(out_dir, "nlrp3_specificity.tsv"))
write_tsv(filter(sel$freq_cmp, delta > 0.05) %>% arrange(desc(delta)),
          file.path(out_dir, "nlrp3_frequency_discordant_variants.tsv"))
# signals per readout (Sup Fig 1)
write_tsv(sel$components, file.path(out_dir, "nlrp3_signal_components.tsv"))

print(as.data.frame(select(out, SNP, pos_hg38, A1, A2, eaf, beta_exposure, se_exposure)),
      row.names = FALSE)
