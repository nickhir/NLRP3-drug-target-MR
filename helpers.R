# functions used by all the scripts, source config.R first
suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
})


# summary statistics

# read a sumstats file and give it our standard column names
read_sumstats <- function(cfg, p_as_text = FALSE) {
    cols <- c(chrom = cfg$chr_col, pos = cfg$pos_col, ea = cfg$ea_col,
              oa = cfg$oa_col, marker = cfg$marker_col, beta = cfg$effect_col,
              se = cfg$se_col, eaf = cfg$eaf_col, p = cfg$p_col,
              rsid = cfg$rsid_col, n = cfg$n_col, gene = cfg$gene_col,
              nlog10 = cfg$nlog10_col, ci_lower = cfg$ci_lower_col,
              ci_upper = cfg$ci_upper_col)
    # chrom as text ("1" and "chr1" both occur)
    df <- fread(cfg$file, sep = if (is.null(cfg$sep)) "\t" else cfg$sep,
                select = unname(cols), data.table = FALSE, showProgress = FALSE,
                colClasses = list(character = c(cfg$chr_col, if (p_as_text) cfg$p_col),
                                  numeric = cfg$pos_col)) %>%
        select(all_of(cols)) %>%
        mutate(chrom = sub("^chr", "", chrom))

    if (!is.null(cfg$marker_col)) {
        df$oa <- sub("^.*_([^/]+)/.*$", "\\1", df$marker)
        df$ea <- sub("^.*/", "", df$marker)
    }
    df
}

# one region (and one gene for eQTL files)
read_region <- function(cfg, chr, start, end) {
    df <- read_sumstats(cfg) %>% filter(chrom == chr, pos >= start, pos <= end)
    if (!is.null(cfg$gene)) df <- filter(df, gene == cfg$gene)
    mutate(df, pos = as.integer(pos))
}

# like read_region() but always on GRCh38, GRCh37 files get lifted with the UCSC chain
# (only keeps variants that map uniquely to the same chromosome)
read_region_grch38 <- function(cfg, chr, start, end) {
    if (!identical(cfg$build, "GRCh37")) return(read_region(cfg, chr, start, end))
    df <- read_region(cfg, chr, 1, Inf)
    lifted <- rtracklayer::liftOver(
        GenomicRanges::GRanges(paste0("chr", chr), IRanges::IRanges(df$pos, width = 1)),
        rtracklayer::import.chain(hg19_to_hg38_chain))
    one <- lengths(lifted) == 1
    hit <- as.data.frame(unlist(lifted[one]))
    same_chr <- hit$seqnames == paste0("chr", chr)
    df <- df[one, ][same_chr, ]
    df$pos <- hit$start[same_chr]
    filter(df, pos >= start, pos <= end)
}

read_genome <- function(cfg) mutate(read_sumstats(cfg, p_as_text = TRUE), pos = as.integer(pos))

# pull variants out of a genome-wide file (joined on position)
lookup_at <- function(cfg, want) {
    read_genome(cfg) %>%
        inner_join(distinct(want, chrom = as.character(chrom), pos = as.integer(pos)),
                   by = c("chrom", "pos")) %>%
        harmonise(cfg) %>%
        filter(SNPid %in% want$SNPid)
}

se_from_ci <- function(ci_lower, ci_upper) (log(ci_upper) - log(ci_lower)) / (2 * qnorm(0.975))

# SE from the p-value (deCODE)
se_from_p <- function(beta, p) abs(beta) / qnorm(p / 2, lower.tail = FALSE)

# -log10(p) from the text so tiny p-values don't underflow
nlog10_from_p <- function(x) {
    x <- trimws(as.character(x))
    out <- -log10(as.numeric(x))
    m <- regmatches(x, regexec("^([0-9.]+)[eE]([+-]?[0-9]+)$", x))
    hit <- lengths(m) == 3
    out[hit] <- -log10(as.numeric(vapply(m[hit], `[`, "", 2))) -
        as.numeric(vapply(m[hit], `[`, "", 3))
    out
}

# project IDs (CHR_POS_A1_A2, A1 = ASCII-first allele), beta and eaf on A1
harmonise <- function(df, cfg) {
    chr <- as.integer(df$chrom)
    ea <- toupper(df$ea)
    oa <- toupper(df$oa)
    A1 <- pmin(ea, oa)
    A2 <- pmax(ea, oa)
    beta <- as.numeric(df$beta)
    if (identical(cfg$effect_type, "OR")) beta <- log(beta)
    p <- as.numeric(df$p)
    se <- switch(if (is.null(cfg$se_source)) "column" else cfg$se_source,
                 column = as.numeric(df$se),
                 ci = se_from_ci(as.numeric(df$ci_lower), as.numeric(df$ci_upper)),
                 p = se_from_p(beta, p))
    nlog10 <- if (isTRUE(cfg$neglog10_p)) p else nlog10_from_p(df$p)
    if (!is.null(cfg$nlog10_col)) {
        alt <- as.numeric(df$nlog10)
        nlog10 <- ifelse(is.finite(alt) & (!is.finite(nlog10) | alt > nlog10), alt, nlog10)
    }
    tibble(
        SNPid = paste(chr, df$pos, A1, A2, sep = "_"),
        chr = chr,
        pos = as.integer(df$pos),
        A1 = A1,
        A2 = A2,
        beta = ifelse(ea == A1, beta, -beta),
        se = se,
        p = if (isTRUE(cfg$neglog10_p)) 10^(-p) else p,
        nlog10 = nlog10,
        eaf = if (is.null(cfg$eaf_col)) NA_real_ else
            ifelse(ea == A1, as.numeric(df$eaf), 1 - as.numeric(df$eaf)),
        n = if (is.null(cfg$n_col)) NA_real_ else as.numeric(df$n)
    ) %>%
        filter(!is.na(beta), !is.na(se), se > 0,
               A1 %in% c("A", "C", "G", "T"), A2 %in% c("A", "C", "G", "T")) %>%
        distinct(SNPid, .keep_all = TRUE)
}

# phenotype SD from the SEs (like coloc does), to put a trait on an SD scale
sd_from_se <- function(se, eaf, n) {
    ok <- !is.na(se) & !is.na(eaf) & !is.na(n)
    coloc:::sdY.est(vbeta = se[ok]^2, maf = pmin(eaf[ok], 1 - eaf[ok]),
                    n = round(median(n[ok])))
}


# instruments

# load the instruments and score from step 00, per unit decrease unless negate = FALSE
load_instruments <- function(negate = TRUE) {
    hg19 <- tribble(
        ~SNP              , ~pos_hg19 ,
        "1_247406019_C_T" , 247569321 ,
        "1_247432548_C_T" , 247595850 ,
        "1_247433558_A_C" , 247596860 ,
        "1_247438293_C_T" , 247601595 ,
        "1_247442302_A_G" , 247605604 ,
        "1_247452478_A_G" , 247615780 ,
        "1_247459572_C_T" , 247622874 ,
        "1_247460342_C_G" , 247623644
    )
    fread(instrument_file, data.table = FALSE) %>%
        transmute(SNP, chr, pos_hg38, A1, A2, eaf_exposure = eaf,
                  beta_exposure = if (negate) -beta_exposure else beta_exposure,
                  se_exposure) %>%
        left_join(hg19, by = "SNP")
}

# rsIDs from the mapping file, otherwise by position from the readout files
lookup_rsids <- function(snps, raw) {
    rsid_map <- fread(rsid_map_file, data.table = FALSE)
    fallback <- lapply(raw, transmute, pos_hg38 = pos, rsid2 = rsid) %>%
        bind_rows() %>%
        filter(grepl("^rs", rsid2)) %>%
        distinct(pos_hg38, .keep_all = TRUE)
    tibble(SNP = snps, pos_hg38 = as.integer(snp_field(snps, 2))) %>%
        mutate(rsid = sapply(SNP, function(s) {
            f <- strsplit(s, "_")[[1]]
            ids <- sprintf("chr%s:%s:%s:%s", f[1], f[2], c(f[3], f[4]), c(f[4], f[3]))
            rsid_map$rsid[rsid_map$variant_id %in% ids][1]
        })) %>%
        left_join(fallback, by = "pos_hg38") %>%
        transmute(SNP, rsid = coalesce(rsid, rsid2))
}

# joint F and R2 for a set of instruments (NA for the score, it has no n)
joint_F_R2 <- function(beta, se, n, ld) {
    K <- length(beta)
    ld_inv <- solve(ld)
    z <- beta / se
    chi2 <- as.numeric(t(z) %*% ld_inv %*% z)
    F_joint <- chi2 / K
    rho <- sign(z) * sqrt(z^2 / (z^2 + n - 2))
    list(
        z = z, chi2 = chi2, F_joint = F_joint, F_mean = mean(z^2), F_min = min(z^2),
        F_per_snp = z^2,
        R2_joint = F_joint / (F_joint + (n - K - 1) / K),
        R2_rho = as.numeric(t(rho) %*% ld_inv %*% rho),
        R2_adj = (chi2 - K) / (chi2 - K + n - K - 1),
        R2_per_snp = z^2 / (z^2 + n - 2)
    )
}


# LD from the INTERVAL panel

# panel IDs are chr1:POS:A1:A2, ours 1_POS_A1_A2
to_panel_id <- function(x) paste0("chr", gsub("_", ":", x))
from_panel_id <- function(x) gsub(":", "_", sub("^chr", "", x))

# field of a variant ID: 2 = pos, 3 = A1, 4 = A2
snp_field <- function(snp, i) sapply(strsplit(snp, "_"), `[`, i)

# signed LD (r relative to each variant's A1)
interval_ld_matrix <- function(snp_ids) {
    ids <- to_panel_id(snp_ids)
    tmp <- tempfile()
    writeLines(ids, paste0(tmp, ".extract"))
    write.table(data.frame(ids, sub("^chr[^:]+:[0-9]+:([^:]+):.*$", "\\1", ids)),
                paste0(tmp, ".ref"), sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    system2(plink2_bin, c("--bfile", ld_panel, "--extract", paste0(tmp, ".extract"),
                          "--ref-allele", "force", paste0(tmp, ".ref"), "2", "1",
                          "--make-pgen", "--out", tmp, "--threads", 2), stdout = FALSE)
    system2(plink2_bin, c("--pfile", tmp, "--r-unphased", "square", "ref-based",
                          "--out", tmp, "--threads", 2), stdout = FALSE)
    ld <- as.matrix(read.table(paste0(tmp, ".unphased.vcor1")))
    rownames(ld) <- colnames(ld) <-
        from_panel_id(read.table(paste0(tmp, ".unphased.vcor1.vars"))[, 1])
    ld[snp_ids, snp_ids]
}

# plink2 clumping
ld_clump_local <- function(variants, bfile, r2, kb, extra = NULL) {
    fn <- tempfile()
    write.table(data.frame(SNP = variants$SNP, P = variants$p), fn,
                row.names = FALSE, quote = FALSE)
    system2(plink2_bin, c(if (file.exists(paste0(bfile, ".pgen"))) "--pfile" else "--bfile",
                          bfile, "--clump", fn, "--clump-p1", 5e-8, "--clump-r2", r2,
                          "--clump-kb", kb, "--out", fn, extra),
            stdout = FALSE, stderr = FALSE)
    fread(paste0(fn, ".clumps"))
}

# all variants within kb of a lead at r2 >= r2
get_high_ld_snps <- function(leads, reference, r2, kb) {
    out <- tempfile()
    writeLines(leads, paste0(out, ".snps"))
    system2(plink2_bin, c("--pfile", reference, "--r2-unphased", "--ld-snp-list",
                          paste0(out, ".snps"), "--ld-window-r2", r2,
                          "--ld-window-kb", kb, "--threads", 4, "--out", out),
            stdout = FALSE)
    fread(paste0(out, ".vcor"), data.table = FALSE)
}

# LD proxies within PROXY_KB at r2 >= PROXY_R2, best first
# (r is signed relative to A1, so sign(r) aligns the proxy's effect)
proxy_candidates <- function(snps) {
    out <- tempfile()
    writeLines(to_panel_id(snps), paste0(out, ".snps"))
    pos <- as.integer(snp_field(snps, 2))
    system2(plink2_bin, c("--bfile", ld_panel, "--chr", snp_field(snps[1], 1),
                          "--from-bp", min(pos) - PROXY_KB * 1000L,
                          "--to-bp", max(pos) + PROXY_KB * 1000L,
                          "--r-unphased", "ref-based", "cols=chrom,pos,id,ref,alt",
                          "--ld-snp-list", paste0(out, ".snps"),
                          "--ld-window-r2", PROXY_R2, "--ld-window-kb", PROXY_KB,
                          "--threads", 2, "--out", out), stdout = FALSE)
    fread(paste0(out, ".vcor"), data.table = FALSE) %>%
        transmute(SNP = from_panel_id(ID_A), proxy = from_panel_id(ID_B),
                  # flip r where REF is A2
                  r = UNPHASED_R * ifelse(REF_A == snp_field(SNP, 3), 1, -1) *
                      ifelse(REF_B == snp_field(proxy, 3), 1, -1),
                  r2 = r^2,
                  kb = abs(as.integer(snp_field(proxy, 2)) - as.integer(snp_field(SNP, 2))) / 1e3) %>%
        filter(!proxy %in% snps, nchar(snp_field(proxy, 3)) == 1, nchar(snp_field(proxy, 4)) == 1) %>%
        arrange(SNP, desc(r2), kb)
}


# instrument selection

# one readout over the instrument window, on A1 (indels kept)
prepare_readout <- function(cfg, chr, start, end) {
    df <- read_region(cfg, chr, start, end)
    A1 <- ifelse(df$oa < df$ea, df$oa, df$ea)
    A2 <- ifelse(df$oa < df$ea, df$ea, df$oa)
    flip <- A1 != df$ea
    p <- as.numeric(df$p)
    tibble(
        SNP = paste(df$chrom, df$pos, A1, A2, sep = "_"),
        position_hg38 = df$pos,
        A1 = A1,
        A2 = A2,
        A1_freq = if (is.null(cfg$eaf_col)) NA_real_ else
            ifelse(flip, 1 - as.numeric(df$eaf), as.numeric(df$eaf)),
        beta = ifelse(flip, -df$beta, df$beta),
        se = as.numeric(df$se),
        p = if (isTRUE(cfg$neglog10_p)) 10^(-p) else p,
        n = cfg$n,
        rsid_source = if (is.null(cfg$rsid_col)) NA_character_ else as.character(df$rsid)
    ) %>%
        filter(!grepl("D|I", SNP), !is.na(beta), !is.na(se), se > 0, !is.na(p)) %>%
        distinct(SNP, .keep_all = TRUE)
}

# drop variants associated with another gene's expression, before and
# after conditioning on that gene's lead eQTL (COJO)
eqtl_specificity <- function(snps, eqtl, chr) {
    rows <- fread(eqtl$file, sep = "\t", data.table = FALSE,
                  select = c("phenotype_id", "chr", "pos_b38", "effect_allele",
                             "other_allele", "slope", "slope_se", "pval_nominal"))
    # one gene at a time
    by_gene <- function(r) {
        lapply(split(r, r$phenotype_id), function(g) {
            harmonise(data.frame(chrom = g$chr, pos = g$pos_b38, ea = g$effect_allele,
                                 oa = g$other_allele, beta = g$slope,
                                 se = g$slope_se, p = g$pval_nominal), list()) %>%
                mutate(gene = g$phenotype_id[1])
        }) %>% bind_rows()
    }
    hits <- by_gene(rows[rows$pos_b38 %in% as.integer(snp_field(snps, 2)), ]) %>%
        filter(SNPid %in% snps, gene != eqtl$gene, p < EQTL_SPECIFICITY_P)
    gene_rows <- by_gene(rows[rows$phenotype_id %in% hits$gene, ])

    detail <- lapply(unique(hits$gene), function(g) {
        d <- gene_rows[gene_rows$gene == g, ]
        tag <- tempfile(g)
        # PLINK1 files for COJO
        system2(plink2_bin, c("--bfile", ld_panel, "--chr", chr,
                              "--from-bp", min(d$pos), "--to-bp", max(d$pos),
                              "--make-bed", "--out", tag, "--threads", 4),
                stdout = FALSE, stderr = FALSE)
        system2(plink2_bin, c("--bfile", tag, "--freq", "--out", tag, "--threads", 4),
                stdout = FALSE, stderr = FALSE)
        afreq <- fread(paste0(tag, ".afreq"), data.table = FALSE)
        # A1 (ASCII-first allele) is REF in the panel
        a1 <- sub("^chr[^:]+:[0-9]+:([^:]+):.*$", "\\1", afreq$ID)
        freq <- data.frame(SNPid = from_panel_id(afreq$ID),
                           freq = ifelse(a1 != afreq$ALT, 1 - afreq$ALT_FREQS, afreq$ALT_FREQS))
        ma <- inner_join(d, freq, by = "SNPid") %>%
            filter(freq > 0, freq < 1) %>%
            # z-score instead of p, p underflows to 0 for the strongest eQTLs
            mutate(p = pmax(2 * pnorm(-abs(beta / se)), 1e-300))
        lead <- ma$SNPid[which.max(abs(ma$beta / ma$se))]
        write.table(transmute(ma, SNP = to_panel_id(SNPid), A1, A2, freq, b = beta,
                              se, p, N = eqtl$n),
                    paste0(tag, ".ma"), sep = "\t", quote = FALSE, row.names = FALSE)
        writeLines(to_panel_id(lead), paste0(tag, ".cond"))
        system2(gcta_bin, c("--bfile", tag, "--cojo-file", paste0(tag, ".ma"),
                            "--cojo-cond", paste0(tag, ".cond"), "--out", tag,
                            "--thread-num", 4),
                stdout = FALSE, stderr = FALSE)
        cma <- fread(paste0(tag, ".cma.cojo"), data.table = FALSE)
        h <- hits[hits$gene == g, ]
        tibble(SNP = h$SNPid, gene = g, p = h$p, lead = lead,
               p_lead = ma$p[ma$SNPid == lead],
               p_conditional = cma$pC[match(to_panel_id(h$SNPid), cma$SNP)],
               excluded = p_conditional < EQTL_SPECIFICITY_P)
    }) %>% bind_rows()
    list(excluded = unique(detail$SNP[detail$excluded]), detail = detail)
}

# cut the INTERVAL panel to the window + 300 kb, with our variant IDs
# (needs shQuote, otherwise the shell reads # as a comment), and drop variants
# whose frequency is off from the reference GWAS by more than 0.05
instrument_panel <- function(freq_readout, chr, start, end) {
    work <- tempfile()
    system2(plink2_bin, c("--bfile", ld_panel, "--chr", chr,
                          "--from-bp", as.integer(start - 300e3),
                          "--to-bp", as.integer(end + 300e3),
                          "--set-all-var-ids", shQuote("@_#_$1_$2"),
                          "--new-id-max-allele-len", 200, "--rm-dup", "force-first",
                          "--make-pgen", "--out", paste0(work, "_raw"), "--threads", 4),
            stdout = FALSE, stderr = FALSE)

    system2(plink2_bin, c("--pfile", paste0(work, "_raw"), "--freq", "--out", work,
                          "--threads", 4), stdout = FALSE, stderr = FALSE)
    panel_freq <- read.table(paste0(work, ".afreq"), header = TRUE, comment.char = "") %>%
        transmute(ID, freq_panel = ifelse(snp_field(ID, 3) != ALT, 1 - ALT_FREQS, ALT_FREQS))
    freq_cmp <- inner_join(
        panel_freq,
        prepare_readout(freq_readout, chr, start, end) %>% transmute(ID = SNP, freq_gwas = A1_freq),
        by = "ID"
    ) %>% mutate(delta = abs(freq_panel - freq_gwas))
    ids <- filter(freq_cmp, delta <= 0.05)$ID
    writeLines(ids, paste0(work, ".keep"))
    ref <- paste0(work, "_ref")
    system2(plink2_bin, c("--pfile", paste0(work, "_raw"), "--extract", paste0(work, ".keep"),
                          "--make-pgen", "--out", ref, "--threads", 4),
            stdout = FALSE, stderr = FALSE)

    list(ref = ref, ids = ids, freq_cmp = freq_cmp,
         eaf = panel_freq %>% filter(ID %in% ids) %>% transmute(SNP = ID, eaf_panel = freq_panel))
}

# pick the instruments and build the activity score at one locus
# (readouts in the order eQTL, CRP, GlycA, neutrophil count)
select_instruments <- function(readouts, chr, start, end, r2 = 0.1) {
    names(readouts) <- c("eQTLs", "CRP", "GlycA", "neutro")

    # neutrophil GWAS as the frequency reference
    panel <- instrument_panel(readouts$neutro, chr, start, end)
    ref <- panel$ref

    stats <- lapply(readouts, function(cfg) {
        filter(prepare_readout(cfg, chr, start, end), SNP %in% panel$ids)
    })
    availability <- function(snp) sum(vapply(stats, function(x) snp %in% x$SNP, TRUE))

    # clump each readout, then group each lead with its r2 > 0.95 neighbours within 40 kb
    blocks <- list()
    for (tn in names(stats)) {
        lead <- ld_clump_local(stats[[tn]], ref, r2, 250)$ID
        friends <- get_high_ld_snps(lead, ref, 0.95, 40) %>%
            select(ID_A, ID_B, UNPHASED_R2) %>%
            rbind(data.frame(ID_A = lead, ID_B = lead, UNPHASED_R2 = 1))
        b <- friends %>% group_by(ID_A) %>% group_split() %>% map(~ pull(.x, ID_B))
        blocks[paste0(tn, "_LD_block", seq_along(b))] <- b
    }

    # merge blocks that share a variant
    block_snps <- lapply(blocks, unique)
    adj <- sapply(block_snps, function(i) {
        sapply(block_snps, function(j) length(intersect(i, j)) > 0)
    })
    diag(adj) <- FALSE
    membership <- igraph::components(igraph::graph_from_adjacency_matrix(
        adj, mode = "undirected", diag = FALSE))$membership
    comp_list <- split(names(membership), membership)
    components <- tibble(comp_id = names(comp_list))
    for (tn in names(stats)) {
        components[[tn]] <- vapply(comp_list, function(b) any(grepl(paste0("^", tn), b)), TRUE)
    }

    # one instrument per signal seen in 2+ readouts: prefer variants in all four, then lowest CRP p
    shared <- comp_list[rowSums(components[names(stats)]) >= 2]
    instruments <- unlist(lapply(shared, function(b) {
        snps <- unique(unlist(block_snps[b]))
        av <- vapply(snps, availability, 1)
        stats$CRP %>%
            filter(SNP %in% snps[av == max(av)]) %>%
            slice_min(p, n = 1, with_ties = FALSE) %>%
            pull(SNP)
    }), use.names = FALSE)

    # replace a missing instrument with the best proxy present in all four
    # (r2 within 0.05 counts as a tie, broken on CRP p)
    need <- instruments[vapply(instruments, availability, 1) != 4]
    proxies <- tibble(original = character(), replacement = character(), r2 = numeric())
    if (length(need) > 0) {
        hl <- get_high_ld_snps(need, ref, PROXY_R2, PROXY_KB)
        for (s in unique(hl$ID_A)) {
            cand <- hl[hl$ID_A == s, ]
            cand <- cand[vapply(cand$ID_B, availability, 1) == 4, ]
            cand <- cand[cand$UNPHASED_R2 >= max(cand$UNPHASED_R2) - 0.05, ]
            best <- stats$CRP %>%
                filter(SNP %in% cand$ID_B) %>%
                slice_min(p, n = 1, with_ties = FALSE) %>%
                pull(SNP)
            proxies <- bind_rows(proxies, tibble(
                original = s, replacement = best,
                r2 = cand$UNPHASED_R2[match(best, cand$ID_B)]))
        }
        swap <- match(instruments, proxies$original)
        instruments <- ifelse(is.na(swap), instruments, proxies$replacement[swap])
    }
    instruments <- unique(instruments)

    specificity <- eqtl_specificity(instruments, readouts$eQTLs, chr)
    instruments <- setdiff(instruments, specificity$excluded)

    # PC1 of the standardised effects
    effect_matrix <- function(field) {
        lapply(names(stats), function(tn) {
            stats[[tn]] %>% filter(SNP %in% instruments) %>% mutate(trait = tn)
        }) %>%
            bind_rows() %>%
            select(SNP, trait, all_of(field)) %>%
            pivot_wider(names_from = trait, values_from = all_of(field)) %>%
            column_to_rownames("SNP") %>%
            drop_na() %>%
            as.matrix()
    }
    beta_m <- effect_matrix("beta")
    se_m <- effect_matrix("se")[rownames(beta_m), colnames(beta_m), drop = FALSE]
    pca <- prcomp(beta_m, center = TRUE, scale. = TRUE)
    loadings <- pca$rotation[, 1]
    if (all(loadings < 0)) loadings <- -loadings
    beta_latent <- as.numeric(beta_m %*% loadings)

    # conservative SE, the three blood readouts share participants so we treat
    # their errors as perfectly correlated
    w <- loadings[c("eQTLs", "CRP", "GlycA", "neutro")]
    se_latent <- sqrt(w[1]^2 * se_m[, "eQTLs"]^2 +
        (w[2] * se_m[, "CRP"] + w[3] * se_m[, "GlycA"] + w[4] * se_m[, "neutro"])^2)

    # scaled to CRP units, higher score = higher expression
    anchor <- tibble(SNP = rownames(beta_m), beta_latent, se_latent) %>%
        inner_join(select(stats$CRP, SNP, beta_CRP = beta, se_CRP = se), by = "SNP")
    k <- coef(lm(beta_CRP ~ 0 + beta_latent, data = anchor, weights = 1 / se_CRP^2))[["beta_latent"]]
    score <- tibble(SNP = rownames(beta_m), beta = beta_latent * k, se = se_latent * abs(k))
    flipped <- cor(score$beta, stats$eQTLs$beta[match(score$SNP, stats$eQTLs$SNP)]) < 0
    if (flipped) {
        score$beta <- -score$beta
        k <- -k
    }

    list(score = score, stats = stats, freq_cmp = panel$freq_cmp, panel_eaf = panel$eaf,
         components = components, proxies = proxies,
         specificity = specificity$detail, loadings = loadings,
         pc1_var = summary(pca)$importance[2, 1], k = k, flipped = flipped)
}

# expression-only instruments for step 11b, the eQTL is the only readout
select_expression_instruments <- function(eqtl, freq_readout, chr, start, end, r2 = 0.1) {
    panel <- instrument_panel(freq_readout, chr, start, end)
    stats <- filter(prepare_readout(eqtl, chr, start, end), SNP %in% panel$ids)
    leads <- ld_clump_local(stats, panel$ref, r2, 250)$ID
    specificity <- eqtl_specificity(leads, eqtl, chr)
    list(instruments = setdiff(leads, specificity$excluded), specificity = specificity$detail)
}


# MR

# correlated IVW, weighted median, Egger intercept and Cochran's Q
run_mr <- function(dat, ld_full) {
    d <- filter(dat, !is.na(beta_outcome), !is.na(se_outcome), se_outcome > 0)
    mr_in <- MendelianRandomization::mr_input(
        bx = d$beta_exposure, bxse = d$se_exposure,
        by = d$beta_outcome, byse = d$se_outcome,
        snps = d$SNP, correlation = ld_full[d$SNP, d$SNP]
    )
    ivw <- MendelianRandomization::mr_ivw(mr_in)
    med <- MendelianRandomization::mr_median(mr_in, weighting = "weighted")
    egger <- MendelianRandomization::mr_egger(mr_in)
    tibble(
        outcome = dat$outcome[1],
        nsnp = nrow(d),
        method = c("IVW (LD-corrected)", "Weighted median"),
        estimate = c(ivw$Estimate, med$Estimate),
        se = c(ivw$StdError, med$StdError),
        ci_lower = c(ivw$CILower, med$CILower),
        ci_upper = c(ivw$CIUpper, med$CIUpper),
        p = c(ivw$Pvalue, med$Pvalue),
        # Egger intercept
        egger_intercept = egger$Intercept,
        egger_intercept_p = egger$Pvalue.Int,
        het_q = c(ivw$Heter.Stat[1], NA),
        het_p = c(ivw$Heter.Stat[2], NA)
    )
}

# one UKB-PPP assay at the instrument positions (GENPOS is GRCh38)
read_ppp <- function(path, ex) {
    txt <- system2(tabix_bin, c(shQuote(path), sprintf("%d:%d-%d", ex$chr, ex$pos_hg38, ex$pos_hg38)),
                   stdout = TRUE, stderr = FALSE)
    ppp <- fread(text = paste(txt, collapse = "\n"), sep = "\t", header = FALSE,
                 colClasses = "character",
                 col.names = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO",
                               "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA"))
    ex %>%
        left_join(transmute(ppp, pos_hg38 = as.integer(GENPOS), ea = toupper(ALLELE1),
                            oa = toupper(ALLELE0), b = as.numeric(BETA),
                            se_outcome = as.numeric(SE), eaf = as.numeric(A1FREQ),
                            info = as.numeric(INFO), n_outcome = as.integer(N),
                            p_outcome = 10^(-as.numeric(LOG10P))),
                  by = "pos_hg38", relationship = "one-to-many") %>%
        mutate(beta_outcome = case_when(ea == A1 & oa == A2 ~ b, ea == A2 & oa == A1 ~ -b),
               eaf_outcome = case_when(ea == A1 & oa == A2 ~ eaf, ea == A2 & oa == A1 ~ 1 - eaf)) %>%
        filter(!is.na(beta_outcome), !is.na(se_outcome), se_outcome > 0) %>%
        distinct(SNP, .keep_all = TRUE)
}

# MR on one UKB-PPP assay, the weighted median ignores LD
mr_ppp_assay <- function(path, ex, ld_full) {
    dat <- read_ppp(path, ex)
    mi <- MendelianRandomization::mr_input(
        bx = dat$beta_exposure, bxse = dat$se_exposure,
        by = dat$beta_outcome, byse = dat$se_outcome,
        snps = dat$SNP, correlation = ld_full[dat$SNP, dat$SNP]
    )
    ivw <- MendelianRandomization::mr_ivw(mi, correl = TRUE)
    egger <- MendelianRandomization::mr_egger(mi, correl = TRUE)
    wm <- MendelianRandomization::mr_median(mi, weighting = "weighted")
    list(
        dat = dat,
        res = tibble(
            n_snps = nrow(dat), n_ppp = max(dat$n_outcome), min_info = min(dat$info),
            ivw_beta = ivw$Estimate, ivw_se = ivw$StdError, ivw_pval = ivw$Pvalue,
            ivw_het_q = ivw$Heter.Stat[1], ivw_het_p = ivw$Heter.Stat[2],
            egger_beta = egger$Estimate, egger_se = egger$StdError.Est,
            egger_pval = egger$Pvalue.Est, egger_intercept = egger$Intercept,
            egger_intercept_se = egger$StdError.Int, egger_intercept_pval = egger$Pvalue.Int,
            median_beta = wm$Estimate, median_se = wm$StdError, median_pval = wm$Pvalue
        )
    )
}

# one CAD study over the locus window
read_cad <- function(cfg) {
    read_region(cfg, CHR, LOCUS_START, LOCUS_END) %>%
        transmute(pos, ea = toupper(ea), oa = toupper(oa),
                  b = if (identical(cfg$effect_type, "OR")) log(as.numeric(beta)) else as.numeric(beta),
                  se = if (identical(cfg$se_source, "ci")) se_from_ci(as.numeric(ci_lower), as.numeric(ci_upper))
                       else as.numeric(se))
}

# one CAD study at the instruments, on A1
cad_at <- function(study, ex) {
    ex %>%
        left_join(study, by = c(pos_hg38 = "pos")) %>%
        mutate(by = case_when(ea == A1 & oa == A2 ~ b, ea == A2 & oa == A1 ~ -b)) %>%
        filter(!is.na(by), !is.na(se), se > 0) %>%
        transmute(SNP, bx, bxse, by, byse = se)
}

# fixed-effect meta-analysis per variant
meta_snp <- function(per_study) {
    per_study %>%
        group_by(SNP) %>%
        summarise(bx = bx[1], bxse = bxse[1],
                  by = sum(by / byse^2) / sum(1 / byse^2),
                  byse = sqrt(1 / sum(1 / byse^2)),
                  n_studies = n(), .groups = "drop")
}

# MVMR of the mediators on CAD (07c, 07d), multiplicative random effects
fit_mvmr <- function(design, snps, meds) {
    d <- design[design$SNPid %in% snps, ]
    x <- setNames(d[paste0("beta_sd_", meds)], meds)
    x$by <- d$beta_cad
    f <- lm(reformulate(meds, response = "by", intercept = FALSE),
            data = x, weights = 1 / d$se_cad^2)
    sig <- summary(f)$sigma
    list(beta = coef(f)[meds], V = vcov(f)[meds, meds, drop = FALSE] * (max(sig, 1) / sig)^2,
         sigma = sig, n = nrow(d), Q = sum(residuals(f)^2 / d$se_cad^2), snps = d$SNPid)
}


# locuszoom figures

# r2 with the index variant
add_LD <- function(loc, index_snp) {
    tmp <- tempfile()
    writeLines(loc$data$panel_id, paste0(tmp, ".ids"))
    system2(plink2_bin, c("--chr", CHR, "--extract", paste0(tmp, ".ids"), "--bfile", ld_panel,
                          "--make-pgen", "--out", tmp), stdout = FALSE)
    system2(plink2_bin, c("--pfile", tmp, "--r2-unphased", "--ld-snp", to_panel_id(index_snp),
                          "--ld-window-r2", 0, "--ld-window-kb", 99999,
                          "--ld-window", 99999, "--out", tmp), stdout = FALSE)
    ld <- fread(paste0(tmp, ".vcor"), data.table = FALSE)
    loc$data <- left_join(loc$data, transmute(ld, panel_id = ID_B, ld = UNPHASED_R2),
                          by = "panel_id")
    loc$index_snp <- index_snp
    loc
}

# recombination rates from UCSC's European mirror, the main site's
# redirect breaks rtracklayer 1.62 ("httpss://")
link_recomb_euro <- function(loc) {
    session <- rtracklayer::browserSession("UCSC", url = "https://genome-euro.ucsc.edu/cgi-bin/")
    GenomeInfoDb::genome(session) <- "hg38"
    region <- GenomicRanges::GRanges(sub("^(chr)?", "chr", loc$seqname),
                                     IRanges::IRanges(loc$xrange[1], loc$xrange[2]))
    loc$recomb <- rtracklayer::getTable(
        rtracklayer::ucscTableQuery(session, table = "recomb1000GAvg", range = region))
    loc
}

# thin grey recombination line on a secondary axis (cM/Mb)
add_recomb_line <- function(p, loc, scale_max = 100,
                            ylab = expression(-log[10](P - value))) {
    rec <- loc$recomb
    f <- max(loc$data[[loc$yvar]], na.rm = TRUE) / scale_max
    step <- data.frame(x = c(rbind(rec$start, rec$end)) / 1e6,
                       y = rep(rec$value, each = 2) * f)
    p$layers <- c(list(geom_line(data = step, aes(x = x, y = y), inherit.aes = FALSE,
                                 colour = "grey65", linewidth = 0.25)), p$layers)
    p +
        scale_y_continuous(
            name = ylab, limits = c(0, NA),
            expand = expansion(mult = c(0, 0.18)),
            sec.axis = sec_axis(~ . / f, name = "Recombination rate (cM/Mb)")
        ) +
        theme(
            axis.title.y.right = element_text(colour = "grey45"),
            axis.text.y.right = element_text(colour = "grey45")
        )
}
