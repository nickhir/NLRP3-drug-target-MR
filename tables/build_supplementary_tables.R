#!/usr/bin/env Rscript
# build the supplementary tables workbook
# run from the project root: Rscript tables/build_supplementary_tables.R

suppressPackageStartupMessages({
    library(openxlsx)
    library(MendelianRandomization)
})


# paths
template <- here::here("NLRP3_supplementary_tables_sources.xlsx")
output <- here::here("NLRP3_supplementary_tables.xlsx")
results <- here::here("results")

read_result <- function(path) {
    read.delim(file.path(results, path), check.names = FALSE, stringsAsFactors = FALSE,
               na.strings = c("", "NA"), quote = "\"", comment.char = "")
}
read_result_rds <- function(path) readRDS(file.path(results, path))


# formatting
# CI as "(1.02, 1.45)", 3 significant figures
ci <- function(lo, hi) {
    # keep trailing zeros, drop the bare decimal point on values >= 100
    f <- function(v) trimws(sub("\\.$", "", formatC(signif(v, 3), format = "fg", digits = 3, flag = "#")))
    ifelse(is.na(lo) | is.na(hi), NA_character_, sprintf("(%s, %s)", f(lo), f(hi)))
}

# all IVW estimates are LD-corrected, just make the method label consistent
method_label <- function(x) {
    x <- sub("^IVW \\(LD-corrected\\)$", "IVW", x)
    sub("^Weighted Median$", "Weighted median", x)
}

# Egger intercept on the IVW row only
egger_on_ivw <- function(d, column) {
    ifelse(method_label(d$method) == "IVW", d[[column]], NA)
}

# scientific notation below 0.001
P_FMT <- "[<0.001]0.00E+00;[<0.01]0.0000;0.000"
INT_FMT <- "#,##0"


# sheets
sheets <- list()
add <- function(id, contents, table, fmt = character()) {
    sheets[[id]] <<- list(contents = contents, table = table, fmt = fmt)
    invisible(NULL)
}


# ST01-ST04, ST14: manual tables
# the template has a title in row 1 and the header in row 3
MANUAL <- c(
    ST01 = "GWAS data sources for the traits used to construct and validate the NLRP3 activity score and the IL1Ra positive control.",
    ST02 = "GWAS data sources for coronary artery disease and imaging-defined atherosclerosis.",
    ST03 = "GWAS data sources for the diseases under clinical investigation for NLRP3 inhibition.",
    ST04 = "GWAS data sources for the cardiometabolic and lifestyle traits.",
    ST14 = "Clinical trials of selective NLRP3 inhibitors."
)
for (id in names(MANUAL)) {
    # colNames = FALSE so read.xlsx doesn't mangle the headers
    raw <- read.xlsx(template, sheet = id, startRow = 3,
                     colNames = FALSE, skipEmptyRows = TRUE, skipEmptyCols = TRUE,
                     detectDates = FALSE)
    header_row <- trimws(as.character(unlist(raw[1, ], use.names = FALSE)))
    named <- !is.na(header_row) & nzchar(header_row)
    x <- raw[-1, named, drop = FALSE]
    names(x) <- header_row[named]
    rownames(x) <- NULL

    # restore numbers, thousands separators except for PMIDs
    fmt <- character()
    for (j in names(x)) {
        text <- trimws(as.character(x[[j]]))
        text[text %in% c("", "NA")] <- NA_character_
        as_number <- suppressWarnings(as.numeric(text))
        # only numeric if every entry parses (some PMID cells hold a DOI)
        numeric_column <- any(!is.na(text)) && all(!is.na(as_number[!is.na(text)]))
        x[[j]] <- if (numeric_column) as_number else text
        if (numeric_column && all(is.na(as_number) | as_number == floor(as_number))) {
            fmt[j] <- if (grepl("PMID", j)) "0" else INT_FMT
        }
    }
    if (id == "ST04") {
        x <- x[order(tolower(x[[1]])), , drop = FALSE]
        rownames(x) <- NULL
    }
    add(id, MANUAL[[id]], x, fmt)
}


# shared inputs
inst <- read_result("00_instrument_selection/nlrp3_instruments.tsv")
readouts_long <- read_result("02_instrument_table/instrument_readouts_long.tsv")
rsid <- read_result("02_instrument_table/instrument_rsids.tsv")
proteome <- read_result("09_proteome_mr/ukb_ppp_proteome_mr_results.tsv")
leave_one_out <- read_result("11_mr_sensitivity/leave_one_out.tsv")
stopifnot(setdiff(leave_one_out$dropped, "none") %in% rsid$SNP)  # ST13 labels use rsIDs

ld_inst <-read_result_rds("03_instrument_strength/instrument_ld_matrix.rds")[inst$SNP, inst$SNP]

TRAIT <- c(NLRP3_expression = "NLRP3 expression", CRP = "CRP",
           GlycA = "GlycA", Neutrophil_count = "Neutrophil count")

# step 12 has the imaging traits (ST07) and the other diseases (ST15)
indications <- read_result("12_mr_indications/additional_indications_mr_results.tsv")
coverage <- read_result("12_mr_indications/additional_indications_coverage.tsv")
indication_key <- coverage$key[match(indications$outcome, coverage$outcome)]
IMAGING <- c("cac", "sis", "carotid")


# ST05 instruments
# flipped to the allele that lowers the activity score, as in Fig 2B
flip_snp <- setNames(inst$beta_exposure > 0, inst$SNP)

# Wald P for the score rows
score_rows <- data.frame(SNP = inst$SNP, trait = "NLRP3 activity score",
                         beta = inst$beta_exposure, se = inst$se_exposure,
                         p = 2 * pnorm(-abs(inst$beta_exposure / inst$se_exposure)),
                         stringsAsFactors = FALSE)
st05_rows <- rbind(readouts_long[c("SNP", "trait", "beta", "se", "p")], score_rows)
st05_rows$trait <- ifelse(st05_rows$trait %in% names(TRAIT),
                          TRAIT[st05_rows$trait], st05_rows$trait)

variant <- match(st05_rows$SNP, inst$SNP)
flipped <- unname(flip_snp[st05_rows$SNP])

st05 <- data.frame(
    SNP = st05_rows$SNP,
    rsID = rsid$rsid[match(st05_rows$SNP, rsid$SNP)],
    Chromosome = inst$chr[variant],
    `Position (GRCh38)` = inst$pos_hg38[variant],
    `Effect allele` = ifelse(flipped, inst$A2[variant], inst$A1[variant]),
    `Other allele` = ifelse(flipped, inst$A1[variant], inst$A2[variant]),
    `Effect allele frequency` = ifelse(flipped, 1 - inst$eaf[variant], inst$eaf[variant]),
    Trait = st05_rows$trait,
    Beta = ifelse(flipped, -st05_rows$beta, st05_rows$beta),
    SE = st05_rows$se, `P value` = st05_rows$p,
    check.names = FALSE, stringsAsFactors = FALSE)

TRAIT_ORDER <- c(unname(TRAIT), "NLRP3 activity score")
st05 <- st05[order(st05$`Position (GRCh38)`, match(st05$Trait, TRAIT_ORDER)), ]
add("ST05",
    "The eight cis-NLRP3 instruments and their per-trait associations, oriented to the allele that lowers NLRP3 activity.",
    st05, c(Chromosome = "0", `Position (GRCh38)` = "0", `P value` = P_FMT))


# ST06 activity score validation
biomarkers <- read_result("04_mr_biomarkers/mr_biomarkers.tsv")
# protein names instead of gene names
CYTOKINE <- c(IL1B = "IL-1β", IL18 = "IL-18", IL6 = "IL-6")

add("ST06",
    "Mendelian randomization estimates for the NLRP3 activity score on the inflammatory biomarkers and cytokines.",
    data.frame(
        Exposure = "NLRP3 activity score (lower)",
        Outcome = ifelse(biomarkers$outcome %in% names(CYTOKINE),
                         CYTOKINE[biomarkers$outcome], biomarkers$outcome),
        `# of Instruments` = biomarkers$nsnp,
        Method = method_label(biomarkers$method),
        `Effect estimate (beta)` = biomarkers$estimate,
        SE = biomarkers$se,
        `95% CI` = ci(biomarkers$ci_lower, biomarkers$ci_upper),
        `P value` = biomarkers$p,
        `Egger intercept` = egger_on_ivw(biomarkers, "egger_intercept"),
        `Egger intercept P` = egger_on_ivw(biomarkers, "egger_intercept_p"),
        `Outcome GWAS N` = biomarkers$n_total,
        check.names = FALSE, stringsAsFactors = FALSE),
    c(`# of Instruments` = INT_FMT, `Outcome GWAS N` = INT_FMT, `P value` = P_FMT,
      `Egger intercept P` = P_FMT))


# ST07 CAD and imaging
cad <- read_result("05_mr_cad/cad_meta_studies.tsv")
# the meta-analysis without 1_247460342_C_G is already in ST13
cad <- cad[!grepl("^Meta-analysis \\(excl", cad$study), ]

OUTCOME_OF <- c(`Aragam et al.` = "Coronary artery disease",
                MVP = "Coronary atherosclerosis", `All of Us` = "Coronary atherosclerosis",
                FinnGen = "Coronary atherosclerosis", `Meta-analysis` = "Coronary artery disease")
stopifnot(all(cad$study %in% names(OUTCOME_OF)))

# Egger intercepts, fitted here from the saved per-SNP effects and LD
per_snp <- read_result("05_mr_cad/cad_meta_studies_per_snp.tsv")
egger_of <- t(vapply(unique(cad$study), function(s) {
    d <- per_snp[per_snp$study == s & per_snp$SNP %in% inst$SNP, ]
    fit <- mr_egger(mr_input(bx = d$bx, bxse = d$bxse, by = d$by, byse = d$byse,
                             snps = d$SNP, correlation = ld_inst[d$SNP, d$SNP]))
    c(fit$Intercept, fit$Pvalue.Int)
}, numeric(2)))
cad_ivw <- method_label(cad$method) == "IVW"

cad_rows <- data.frame(
    Outcome = unname(OUTCOME_OF[cad$study]),
    `GWAS source` = cad$study, Method = method_label(cad$method),
    `# of Instruments` = cad$n_snps,
    Scale = "OR",
    `Effect estimate (beta or OR)` = exp(cad$beta),
    `95% CI` = ci(exp(cad$ci_lower), exp(cad$ci_upper)),
    `SE (beta or log OR)` = cad$se, `P value` = cad$pval,
    `Egger intercept` = ifelse(cad_ivw, egger_of[cad$study, 1], NA),
    `Egger intercept P` = ifelse(cad_ivw, egger_of[cad$study, 2], NA),
    check.names = FALSE, stringsAsFactors = FALSE)

imaging <- indications[indication_key %in% IMAGING, ]
# CAC is continuous (inverse-normal), the plaque scores are ordinal (proportional OR)
proportional_or <- indication_key[indication_key %in% IMAGING] != "cac"
imaging_rows <- data.frame(
    Outcome = imaging$outcome, `GWAS source` = "SCAPIS",
    Method = method_label(imaging$method),
    `# of Instruments` = imaging$nsnp,
    Scale = ifelse(proportional_or, "Proportional OR", "Beta (SD)"),
    `Effect estimate (beta or OR)` = ifelse(proportional_or, exp(imaging$estimate),
                                            imaging$estimate),
    `95% CI` = ifelse(proportional_or, ci(exp(imaging$ci_lower), exp(imaging$ci_upper)),
                      ci(imaging$ci_lower, imaging$ci_upper)),
    `SE (beta or log OR)` = imaging$se, `P value` = imaging$p,
    `Egger intercept` = egger_on_ivw(imaging, "egger_intercept"),
    `Egger intercept P` = egger_on_ivw(imaging, "egger_intercept_p"),
    check.names = FALSE, stringsAsFactors = FALSE)

add("ST07",
    "Mendelian randomization estimates for the NLRP3 activity score on coronary artery disease and imaging-defined atherosclerosis.",
    rbind(cad_rows, imaging_rows),
    c(`# of Instruments` = INT_FMT, `P value` = P_FMT, `Egger intercept P` = P_FMT))


# ST08 cardiometabolic traits
cardio <- read_result("06_mr_cardiometabolic/cardiometabolic_mr.tsv")
# only log-odds betas become ORs (GSCAN smoking is in SD units)
log_odds <- cardio$log_odds %in% c(TRUE, "TRUE")
add("ST08",
    "Mendelian randomization estimates for the NLRP3 activity score on cardiometabolic and lifestyle traits.",
    data.frame(
        Outcome = cardio$outcome, Method = method_label(cardio$method),
        `# of Instruments` = cardio$nsnp,
        `Effect estimate (beta or OR)` = ifelse(log_odds, exp(cardio$estimate),
                                                cardio$estimate),
        `95% CI` = ifelse(log_odds, ci(exp(cardio$ci_lower), exp(cardio$ci_upper)),
                          ci(cardio$ci_lower, cardio$ci_upper)),
        `SE (beta or log OR)` = cardio$se, `P value` = cardio$p,
        `Egger intercept` = egger_on_ivw(cardio, "egger_intercept"),
        `Egger intercept P` = egger_on_ivw(cardio, "egger_intercept_p"),
        check.names = FALSE, stringsAsFactors = FALSE),
    c(`# of Instruments` = INT_FMT, `P value` = P_FMT, `Egger intercept P` = P_FMT))


# ST09 mediation
mediation <- read_result("07c_mediation/mediation_summary.tsv")
per_mediator <- read_result("07c_mediation/mediation_per_mediator.tsv")
proportion <- read_result("07c_mediation/proportion_mediated.tsv")
conditional_f <- read_result_rds("07c_mediation/mediation.rds")$cond_F

total <- mediation[mediation$quantity == "Total effect", ]
joint <- mediation[mediation$quantity == "Joint indirect", ]
direct <- mediation[mediation$quantity == "Direct effect", ]
n_mediators <- nrow(per_mediator)

pct_ci <- sprintf("(%s, %s)",
                  formatC(100 * proportion$ci_lower, format = "f", digits = 1),
                  formatC(100 * proportion$ci_upper, format = "f", digits = 1))

add("ST09",
    "Mediation of the coronary artery disease association by systolic blood pressure, apolipoprotein B and type 2 diabetes.",
    data.frame(
        # lower case mid-sentence ("via apolipoprotein B")
        `Effect component` = c("Total effect",
                               paste("Indirect via", sub("^(.)", "\\L\\1", per_mediator$label,
                                                         perl = TRUE)),
                               "Joint indirect effect", "Residual direct effect"),
        # mediator -> CAD columns only for the per-mediator rows
        `Mediator to CAD (log OR)` = c(NA, per_mediator$beta, NA, NA),
        SE = c(NA, per_mediator$se_beta, NA, NA),
        `Conditional F` = c(NA, conditional_f, NA, NA),
        `Effect (log OR)` = c(total$estimate, per_mediator$indirect,
                              joint$estimate, direct$estimate),
        `SE (log OR)` = c(total$se, per_mediator$se_indirect, joint$se, direct$se),
        OR = c(total$or, per_mediator$or, joint$or, direct$or),
        `95% CI (OR)` = c(ci(total$or_lower, total$or_upper),
                          ci(per_mediator$or_lower, per_mediator$or_upper),
                          ci(joint$or_lower, joint$or_upper),
                          ci(direct$or_lower, direct$or_upper)),
        `Proportion mediated (%)` = 100 * c(NA, per_mediator$pm_alone,
                                            proportion$proportion_mediated, NA),
        # Monte Carlo CI only on the joint row
        `Proportion mediated 95% CI (%)` = c(rep(NA_character_, 1 + n_mediators),
                                             pct_ci, NA),
        `P value` = c(total$p, per_mediator$p, joint$p, direct$p),
        check.names = FALSE, stringsAsFactors = FALSE),
    c(`P value` = P_FMT))


# ST10 IL1RN positive control
il1rn <- read_result("08_il1rn_positive_control/il1rn_mr_results.tsv")
add("ST10",
    "Mendelian randomization estimates for the IL1Ra activity score.",
    data.frame(
        Exposure = "IL1Ra activity score (higher)",
        Outcome = il1rn$outcome, `# of Instruments` = il1rn$n_snps,
        Method = method_label(il1rn$method),
        Scale = il1rn$scale,
        # already ORs
        `Effect estimate (beta or OR)` = il1rn$estimate,
        `95% CI` = ci(il1rn$ci_lower, il1rn$ci_upper),
        `SE (beta or log OR)` = il1rn$se, `P value` = il1rn$p,
        check.names = FALSE, stringsAsFactors = FALSE),
    c(`# of Instruments` = INT_FMT, `P value` = P_FMT))


# ST11 proteome-wide MR
ranked <- proteome[order(proteome$ivw_pval), ]
add("ST11",
    "Proteome-wide Mendelian randomization estimates across all tested UK Biobank plasma proteins.",
    data.frame(
        `Protein ID` = ranked$protein_id, `Gene name` = ranked$gene_name,
        `# of Instruments` = ranked$n_snps,
        `IVW beta` = ranked$ivw_beta, `IVW SE` = ranked$ivw_se,
        `IVW P` = ranked$ivw_pval, `IVW P adjusted` = ranked$ivw_pval_bonf,
        `Weighted median beta` = ranked$median_beta,
        `Weighted median SE` = ranked$median_se,
        `Weighted median P` = ranked$median_pval,
        # Egger intercept only, the slope isn't used
        `Egger intercept` = ranked$egger_intercept,
        `Egger intercept P` = ranked$egger_intercept_pval,
        check.names = FALSE, stringsAsFactors = FALSE),
    c(`# of Instruments` = INT_FMT, `IVW P` = P_FMT, `IVW P adjusted` = P_FMT,
      `Weighted median P` = P_FMT, `Egger intercept P` = P_FMT))


# ST12 over-representation analysis
ora <- read_result("10_proteome_ora/ukb_ppp_proteome_ora.tsv")
DIRECTION <- c(lowered_by_inhibition = "Lower with lower NLRP3 activity",
               raised_by_inhibition  = "Higher with lower NLRP3 activity")
ora <- ora[order(match(ora$direction, names(DIRECTION)), ora$p.adjust), ]
add("ST12",
    "Over-representation analysis of the proteins altered by lower NLRP3 activity.",
    data.frame(
        Direction = unname(DIRECTION[ora$direction]),
        `Term ID` = ora$ID, Description = ora$Description,
        # text, otherwise Excel turns the ratios into dates
        `Gene ratio` = ora$GeneRatio, `Background ratio` = ora$BgRatio,
        `# Genes` = ora$Count,
        `P value` = ora$pvalue, `Adjusted P value` = ora$p.adjust,
        Genes = gsub("/", ", ", ora$geneID),
        check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL),
    c(`# Genes` = INT_FMT, `P value` = P_FMT, `Adjusted P value` = P_FMT))


# ST13 sensitivity analyses
sweep <- read_result("11_mr_sensitivity/r2_sweep.tsv")
expression_only <- read_result("11b_expression_only/expression_only_mr.tsv")
expression_only <- expression_only[is.na(expression_only$study) |
                                       expression_only$study == "Meta-analysis", ]
EXPRESSION_SET <- c(`Expression-only instrument` = "Expression-only instrument: per SD lower NLRP3 expression",
                    `Activity score variants` = "Activity score variants: per SD lower NLRP3 expression")
# CAD first
expression_only <- expression_only[order(match(expression_only$set, names(EXPRESSION_SET)),
                                         expression_only$outcome != "CAD"), ]
expression_outcome <- ifelse(expression_only$outcome == "CAD", "Coronary artery disease",
                             ifelse(expression_only$outcome %in% names(CYTOKINE),
                                    CYTOKINE[expression_only$outcome], expression_only$outcome))

sensitivity_block <- function(analysis, outcome, d) {
    # one outcome per row, otherwise ifelse() returns a single value
    outcome <- rep_len(outcome, nrow(d))
    log_odds <- outcome == "Coronary artery disease"
    data.frame(
        Analysis = analysis, Outcome = outcome, Method = method_label(d$method),
        `# of Instruments` = d$n_snps,
        `Effect estimate (beta or OR)` = ifelse(log_odds, exp(d$estimate), d$estimate),
        `95% CI` = ifelse(log_odds, ci(exp(d$ci_lower), exp(d$ci_upper)),
                          ci(d$ci_lower, d$ci_upper)),
        `SE (beta or log OR)` = d$se, `P value` = d$p,
        `Egger intercept` = egger_on_ivw(d, "egger_intercept"),
        `Egger intercept P` = egger_on_ivw(d, "egger_intercept_p"),
        `Cochran's Q` = d$het_q, `Cochran's Q P` = d$het_p,
        check.names = FALSE, stringsAsFactors = FALSE)
}

add("ST13",
    "Sensitivity analyses for the coronary artery disease association, and estimates with NLRP3 expression as the exposure.",
    rbind(
        sensitivity_block(sprintf("LD clumping r2 < %.1f", sweep$r2_threshold),
                          "Coronary artery disease", sweep),
        # rsIDs
        sensitivity_block(ifelse(leave_one_out$dropped == "none", "All instruments",
                                 paste("Leave-one-out: excluding",
                                       rsid$rsid[match(leave_one_out$dropped, rsid$SNP)])),
                          "Coronary artery disease", leave_one_out),
        sensitivity_block(unname(EXPRESSION_SET[expression_only$set]), unname(expression_outcome),
                          expression_only)),
    c(`# of Instruments` = INT_FMT, `P value` = P_FMT, `Egger intercept P` = P_FMT,
      `Cochran's Q P` = P_FMT))


# ST15 exploratory indications
exploratory <- indications[!indication_key %in% IMAGING, ]
add("ST15",
    "Mendelian randomization estimates for the NLRP3 activity score across the exploratory disease indications.",
    data.frame(
        Outcome = exploratory$outcome, Method = method_label(exploratory$method),
        `# of Instruments` = exploratory$nsnp,
        OR = exp(exploratory$estimate),
        `95% CI (OR)` = ci(exp(exploratory$ci_lower), exp(exploratory$ci_upper)),
        `SE (log OR)` = exploratory$se, `P value` = exploratory$p,
        `Egger intercept` = egger_on_ivw(exploratory, "egger_intercept"),
        `Egger intercept P` = egger_on_ivw(exploratory, "egger_intercept_p"),
        check.names = FALSE, stringsAsFactors = FALSE),
    c(`# of Instruments` = INT_FMT, `P value` = P_FMT, `Egger intercept P` = P_FMT))


# write
# zero-padded ids sort in order
ids <- sort(names(sheets))
contents <- data.frame(
    Sheet = ids,
    Contents = vapply(sheets[ids], `[[`, character(1), "contents"),
    check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)

header_style <- createStyle(textDecoration = "bold", valign = "top", wrapText = TRUE)
body_style <- createStyle(valign = "top")

# default formats: 4 decimals for numbers
column_format <- function(name, values, fmt) {
    if (name %in% names(fmt)) fmt[[name]] else if (is.numeric(values)) "0.0000" else NULL
}
# column widths
column_width <- function(name, values) {
    longest_value <- if (is.character(values)) {
        suppressWarnings(max(nchar(values), na.rm = TRUE)) + 2
    } else {
        12
    }
    min(46, max(11, nchar(name) + 2, longest_value))
}

# merge repeated values into one cell (not for the long sheets ST11 and ST12)
MERGE <- list(
    ST05 = c("SNP", "rsID", "Chromosome", "Position (GRCh38)", "Effect allele",
             "Other allele", "Effect allele frequency"),
    ST06 = c("Exposure", "Outcome", "# of Instruments", "Outcome GWAS N"),
    ST07 = c("Outcome", "GWAS source", "# of Instruments", "Scale"),
    ST08 = c("Outcome", "# of Instruments"),
    ST10 = "Exposure",
    ST13 = c("Analysis", "Outcome", "# of Instruments"),
    ST15 = c("Outcome", "# of Instruments")
)
centred <- createStyle(halign = "center", valign = "center", wrapText = TRUE)
merge_runs <- function(id, x, cols) {
    group <- rep("", nrow(x))
    for (name in cols) {
        j <- match(name, names(x))
        group <- paste(group, x[[j]], sep = "\r")
        first <- which(c(TRUE, group[-1] != group[-nrow(x)]))
        last <- c(first[-1] - 1L, nrow(x))
        addStyle(wb, id, centred, rows = seq_len(nrow(x)) + 1L, cols = j, stack = TRUE)
        for (k in which(last > first)) {
            rows <- first[k]:last[k] + 1L
            mergeCells(wb, id, cols = j, rows = rows)
            # only the top cell of a merged range keeps its value
            deleteData(wb, id, cols = j, rows = rows[-1], gridExpand = TRUE)
        }
    }
}

wb <- createWorkbook()
write_sheet <- function(id, x, fmt = character(), widths = NULL) {
    addWorksheet(wb, id)
    writeData(wb, id, x, headerStyle = header_style, keepNA = TRUE, na.string = "NA",
              withFilter = FALSE)
    data_rows <- seq_len(nrow(x)) + 1L
    addStyle(wb, id, body_style, rows = data_rows, cols = seq_along(x), gridExpand = TRUE)

    for (j in seq_along(x)) {
        name <- names(x)[j]
        values <- x[[j]]

        # "NA" as text, right-aligned
        if (is.numeric(values)) {
            addStyle(wb, id, createStyle(halign = "right"), data_rows, j,
                     gridExpand = TRUE, stack = TRUE)
        }

        f <- column_format(name, values, fmt)
        if (!is.null(f)) {
            addStyle(wb, id, createStyle(numFmt = f), data_rows, j,
                     gridExpand = TRUE, stack = TRUE)
            # don't show small nonzero values as 0
            if (f == "0.0000") {
                tiny <- which(!is.na(values) & values != 0 & abs(values) < 5e-5)
                if (length(tiny)) {
                    addStyle(wb, id, createStyle(numFmt = P_FMT), data_rows[tiny], j,
                             stack = TRUE)
                }
            }
        }
        setColWidths(wb, id, j,
                     if (!is.null(widths)) widths[j] else column_width(name, values))
    }
    freezePane(wb, id, firstActiveRow = 2)
}

write_sheet("Contents", contents, widths = c(10, 120))
for (id in ids) {
    write_sheet(id, sheets[[id]]$table, sheets[[id]]$fmt)
    if (id %in% names(MERGE)) merge_runs(id, sheets[[id]]$table, MERGE[[id]])
}

saveWorkbook(wb, output, overwrite = TRUE)

# openxlsx 4.2.8 adds links to drawings that don't exist, which breaks
# openpyxl/pandas, so remove them
dir <- tempfile("xlsxfix_")
utils::unzip(output, exdir = dir)
for (f in list.files(file.path(dir, "xl", "worksheets", "_rels"), full.names = TRUE)) {
    txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
    writeLines(gsub('<Relationship[^>]*Target="\\.\\./drawings/[^"]*"[^>]*/>', "", txt), f)
}
# italic gene names and "cis" (openxlsx can't do rich text, so we patch the XML)
RICH <- list(
    c("", "NLRP3", " expression"),
    c("", "NLRP3", " and ", "IL1RN", " eQTLs"),
    c("The eight ", "cis", "-", "NLRP3",
      " instruments and their per-trait associations, oriented to the allele that lowers NLRP3 activity."),
    c("Sensitivity analyses for the coronary artery disease association, and estimates with ", "NLRP3",
      " expression as the exposure."),
    c("Expression-only instrument: per SD lower ", "NLRP3", " expression"),
    c("Activity score variants: per SD lower ", "NLRP3", " expression")
)
run_xml <- function(text, italic) {
    if (!nzchar(text)) return("")
    sprintf('<r><rPr>%s<sz val="11"/><color rgb="FF000000"/><rFont val="Calibri"/><family val="2"/></rPr><t xml:space="preserve">%s</t></r>',
            if (italic) "<i/>" else "", text)
}
ss_file <- file.path(dir, "xl", "sharedStrings.xml")
ss <- paste(readLines(ss_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
for (segments in RICH) {
    plain <- paste0(segments, collapse = "")
    runs <- paste0(mapply(run_xml, segments, seq_along(segments) %% 2 == 0), collapse = "")
    # <t> or <t xml:space="preserve">
    targets <- sprintf(c("<si><t>%s</t></si>", "<si><t xml:space=\"preserve\">%s</t></si>"), plain)
    hit <- targets[vapply(targets, grepl, logical(1), x = ss, fixed = TRUE)]
    stopifnot(length(hit) == 1)
    ss <- gsub(hit, paste0("<si>", runs, "</si>"), ss, fixed = TRUE)
}
writeLines(ss, ss_file, useBytes = TRUE)

unlink(output)
setwd(dir)
utils::zip(output, list.files(".", recursive = TRUE, all.files = TRUE), flags = "-qX")
