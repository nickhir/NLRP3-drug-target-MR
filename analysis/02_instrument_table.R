# make the instrument table

source(here::here("config.R"))
source(here::here("helpers.R"))

out_dir <- step_dir("02_instrument_table")

# keep the original allele orientation
exposure <- load_instruments(negate = FALSE)
raw <- lapply(READOUTS, read_region, chr = CHR, start = LOCUS_START, end = LOCUS_END)

readouts <- lapply(names(READOUTS), function(k) {
    harmonise(raw[[k]], READOUTS[[k]]) %>%
        filter(SNPid %in% exposure$SNP) %>%
        transmute(SNP = SNPid, trait = k, beta, se, p, n = READOUTS[[k]]$n) %>%
        arrange(match(SNP, exposure$SNP))
}) %>% bind_rows()

aligned <- readouts %>%
    select(SNP, trait, beta, se) %>%
    pivot_wider(names_from = trait, values_from = c(beta, se)) %>%
    inner_join(select(exposure, SNP, chr, pos_hg38, pos_hg19, A1, A2,
                      eaf_exposure, beta_exposure, se_exposure),
               by = "SNP") %>%
    arrange(match(SNP, exposure$SNP))

# rsIDs, the two rarest variants are not in the mapping file
rsids <- lookup_rsids(exposure$SNP, raw[c("Neutrophil_count", "CRP", "GlycA")])

fwrite(readouts, file.path(out_dir, "instrument_readouts_long.tsv"), sep = "\t")
fwrite(aligned, file.path(out_dir, "instrument_table_aligned.tsv"), sep = "\t")
fwrite(rsids, file.path(out_dir, "instrument_rsids.tsv"), sep = "\t")
