# define paths and parameters for all the datasets

analysis_dir <- here::here()
dataset_dir <- file.path(analysis_dir, "datasets")
results_dir <- file.path(analysis_dir, "results")

step_dir <- function(step) {
    d <- file.path(results_dir, step)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    d
}


# external tools and data
plink2_bin <- "/rds/user/nh608/hpc-work/software/plink2/plink2"
gcta_bin <- "/rds/user/nh608/hpc-work/software/gcta/gcta-1.94.1"
tabix_bin <- "/rds/user/nh608/hpc-work/software/micromamba/envs/sambcfenv/bin/tabix"

# hg19 -> hg38 chain
hg19_to_hg38_chain <- "/rds/user/nh608/hpc-work/software/UCSC_liftOver/hg19ToHg38.over.chain"

# UKB-PPP summary statistics, one tabix file per assay
ppp_dir <- paste0(
    "/rds/project/rds-C1Ph08tkaOA/public/proteomics/UKB-PPP/sun23/",
    "UKB-PPP pGWAS summary statistics (reformatted)/Combined_European"
)
ppp_manifest <- paste0(ppp_dir, ".lst")

# INTERVAL WGS LD panel (GRCh38), variants named chr1:POS:A1:A2
ld_panel <- paste0(
    "/rds/user/nh608/hpc-work/oxLDL/data/oxLDL_data/INTERVAL_reference/",
    "WGS_reference/ld_panels/INTERVAL_allchr.GRCh38.alpha_sorted_alleles"
)


# NLRP3 locus (GRCh38)
CHR <- 1
GENE_START <- 247416156 # NLRP3 gene body
GENE_END <- 247449108

# +/-150 kb instrument window
INSTRUMENT_START <- as.integer(GENE_START - 150e3)
INSTRUMENT_END <- as.integer(GENE_END + 150e3)

# +/-200 kb for colocalisation and outcome lookups
LOCUS_START <- as.integer(GENE_START - 200e3)
LOCUS_END <- as.integer(GENE_END + 200e3)

# +/-1 Mb for the flank colocalisation
FLANK_START <- as.integer(GENE_START - 1e6)
FLANK_END <- as.integer(GENE_END + 1e6)

NLRP3_ENSG <- "ENSG00000162711"

# IL1RN (Ensembl 111)
IL1RN_CHR <- 2
IL1RN_START <- 113099315
IL1RN_END <- 113134016
IL1RN_ENSG <- "ENSG00000136689"

# drop variants associated with another gene's expression at this P,
# before and after conditioning on that gene's lead eQTL
EQTL_SPECIFICITY_P <- 1e-3

# LD proxies: r2 >= PROXY_R2 within PROXY_KB kb
PROXY_R2 <- 0.9
PROXY_KB <- 50


# instruments
instrument_file <- file.path(
    results_dir,
    "00_instrument_selection",
    "nlrp3_instruments.tsv"
)

# variant ID -> rsID
rsid_map_file <- file.path(dataset_dir, "rsid_mapping.tsv.gz")

# everything is per unit lower activity score (i.e. inhibition)
EXPOSURE_DIRECTION <- "1-unit DECREASE in cis-NLRP3 activity score"


# GWAS Catalog harmonised format
GWAS_CATALOG <- list(
    build = "GRCh38",
    chr_col = "chromosome",
    pos_col = "base_pair_location",
    ea_col = "effect_allele",
    oa_col = "other_allele",
    effect_col = "beta",
    effect_type = "beta",
    se_source = "column",
    se_col = "standard_error",
    eaf_col = "effect_allele_frequency",
    p_col = "p_value"
)
# same, with OR and 95% CI (SE from the CI)
GWAS_CATALOG_OR <- modifyList(GWAS_CATALOG, list(
    effect_col = "odds_ratio",
    effect_type = "OR",
    se_source = "ci",
    se_col = NULL,
    ci_lower_col = "ci_lower",
    ci_upper_col = "ci_upper"
))
# older hm_ column names
GWAS_CATALOG_HM <- modifyList(GWAS_CATALOG, list(
    chr_col = "hm_chrom",
    pos_col = "hm_pos",
    ea_col = "hm_effect_allele",
    oa_col = "hm_other_allele",
    effect_col = "hm_beta",
    eaf_col = "hm_effect_allele_frequency"
))

T2DGGI <- list(
    chr_col = "Chromsome",
    pos_col = "Position",
    ea_col = "EffectAllele",
    oa_col = "NonEffectAllele",
    effect_col = "Beta",
    se_col = "SE",
    eaf_col = "EAF",
    p_col = "Pval"
)
# GSCAN: ALT is the effect allele, no allele frequencies
GSCAN <- list(
    build = "GRCh37",
    chr_col = "CHROM",
    pos_col = "POS",
    ea_col = "ALT",
    oa_col = "REF",
    effect_col = "BETA",
    se_col = "SE",
    p_col = "PVALUE"
)


# the four readouts used to build the score
READOUTS <- list(
    NLRP3_expression = list(
        label = "NLRP3 expression",
        file = file.path(dataset_dir, "nlrp3_eqtl_interval_chr1.tsv"),
        n = 4732,
        chr_col = "chr",
        pos_col = "pos_b38",
        ea_col = "effect_allele",
        oa_col = "other_allele",
        effect_col = "slope",
        se_col = "slope_se",
        p_col = "pval_nominal",
        eaf_col = "af",
        # all genes in the window
        gene_col = "phenotype_id",
        gene = NLRP3_ENSG
    ),
    CRP = modifyList(GWAS_CATALOG_HM, list(
        label = "CRP concentration",
        file = file.path(dataset_dir, "crp_GCST90029070.h.tsv.gz"),
        n = 575531,
        eaf_col = NULL,
        rsid_col = "hm_rsid"
    )),
    GlycA = modifyList(GWAS_CATALOG, list(
        label = "GlycA concentration",
        file = file.path(dataset_dir, "glyca_GCST90497330.h.tsv.gz"),
        n = 434646,
        p_col = "neg_log_10_p_value",
        neglog10_p = TRUE,
        rsid_col = "rsid"
    )),
    # frequency reference for the panel check
    Neutrophil_count = modifyList(GWAS_CATALOG_HM, list(
        label = "Neutrophil count",
        file = file.path(dataset_dir, "neutrophil_count_GCST90002351.h.tsv.gz"),
        n = 519288,
        rsid_col = "hm_rsid"
    ))
)

# same readouts at IL1RN
IL1RN_READOUTS <- list(
    modifyList(READOUTS$NLRP3_expression, list(
        label = "IL1RN expression",
        file = file.path(dataset_dir, "il1rn_eqtl_interval_chr2.tsv"),
        gene = IL1RN_ENSG,
        rsid_col = "variant_id"
    )),
    READOUTS$CRP,
    modifyList(READOUTS$GlycA, list(p_col = "p_value", neglog10_p = NULL)),
    READOUTS$Neutrophil_count
)

# IL1Ra protein
ppp_il1rn_assay <- file.path(ppp_dir, "IL1RN_P18510_OID20700_v1_Inflammation.bgz")


# disease and trait outcomes
outcome_entry <- function(template, label, file, ...) {
    modifyList(template, list(label = label, file = file.path(dataset_dir, file), ...))
}
OUTCOMES <- list(
    cac = outcome_entry(GWAS_CATALOG, "Coronary artery calcification (Agatston, CT)",
                  "coronary_artery_calcification_GCST90503074.h.tsv.gz"),
    sis = outcome_entry(GWAS_CATALOG, "Coronary plaque burden (SIS, CTA)",
                  "coronary_plaque_burden_GCST90503075.h.tsv.gz"),
    carotid = outcome_entry(GWAS_CATALOG, "Carotid plaque (ultrasound)",
                      "carotid_plaque_burden_GCST90503076.h.tsv.gz"),
    allergic_rhinitis = outcome_entry(GWAS_CATALOG_OR, "Allergic rhinitis",
                                "allergic_rhinitis_GCST90476021.h.tsv.gz",
                                n_cases = 92310, n_controls = 311377),
    als = outcome_entry(GWAS_CATALOG_HM, "Amyotrophic lateral sclerosis",
                  "amyotrophic_lateral_sclerosis_GCST90027164.h.tsv.gz",
                  n_cases = 27205, n_controls = 110881),
    asthma = outcome_entry(GWAS_CATALOG, "Asthma", "asthma_GCST90399686.h.tsv.gz",
                     n_cases = 121940, n_controls = 1254131),
    chronic_airway_obstruction = outcome_entry(GWAS_CATALOG_OR, "Chronic airway obstruction",
                                         "chronic_airway_obstruction_GCST90476027.h.tsv.gz",
                                         n_cases = 103054, n_controls = 315450),
    gout = outcome_entry(GWAS_CATALOG_OR, "Gout", "gout_GCST90475735.h.tsv.gz",
                   n_cases = 42034, n_controls = 397989),
    hfref = outcome_entry(GWAS_CATALOG_OR, "Heart failure (reduced EF)",
                    "heart_failure_reduced_ef_GCST90475983.h.tsv.gz",
                    n_cases = 31220, n_controls = 407213),
    # b38_* columns are GRCh38 despite the file name
    covid19 = outcome_entry(GWAS_CATALOG, "Hospitalized COVID-19",
                      "hospitalized_covid19_HGI_B2_release7_GRCh37.tsv.gz",
                      chr_col = "b38_chr", pos_col = "b38_pos", ea_col = "b38_alt",
                      oa_col = "b38_ref", effect_col = "all_inv_var_meta_beta",
                      se_col = "all_inv_var_meta_sebeta", eaf_col = "all_meta_AF",
                      p_col = "all_inv_var_meta_p",
                      n_cases = 32519, n_controls = 2062805),
    # UK Biobank WGS, ICD-10 M17, non-Finnish Europeans
    # two instruments are missing and get LD proxies
    knee_oa = outcome_entry(GWAS_CATALOG, "Knee osteoarthritis",
                      "knee_osteoarthritis_GCST90474022.h.tsv.gz",
                      n_cases = 44190, n_controls = 414250),
    mi = outcome_entry(GWAS_CATALOG_OR, "Myocardial infarction",
                 "myocardial_infarction_GCST90475932.h.tsv.gz",
                 n_cases = 39074, n_controls = 392979),
    obesity = outcome_entry(GWAS_CATALOG_OR, "Obesity", "obesity_GCST90475762.h.tsv.gz",
                      n_cases = 169600, n_controls = 244254),
    parkinsons = outcome_entry(GWAS_CATALOG, "Parkinson's disease", "parkinsons_disease_GCST009325.tsv",
                         build = "GRCh37", n_cases = 33674, n_controls = 449056),
    # no SE or CI, so we get the SE from the p-value
    pericarditis = outcome_entry(GWAS_CATALOG_OR, "Pericarditis", "pericarditis_decode.txt.gz",
                           chr_col = "Chrom", pos_col = "Pos", ea_col = "Effect_Allele",
                           oa_col = "Other_allele", # lower-case a
                           effect_col = "comb_Effect", se_source = "p",
                           ci_lower_col = NULL, ci_upper_col = NULL, eaf_col = NULL,
                           p_col = "comb_Pval", n_cases = 4894, n_controls = 1457822),
    ra_ishigaki = outcome_entry(GWAS_CATALOG, "Rheumatoid arthritis",
                          "rheumatoid_arthritis_ishigaki_GCST90132223.h.tsv.gz",
                          n_cases = 22350, n_controls = 74823),
    t2d = outcome_entry(T2DGGI, "Type 2 diabetes", "type2_diabetes_T2DGGI_EUR_GRCh37.sumstats.zip",
                  build = "GRCh37", effect_type = "beta", se_source = "column",
                  n_cases = 242283, n_controls = 1569734),
    uc = outcome_entry(GWAS_CATALOG_HM, "Ulcerative colitis", "ulcerative_colitis_GCST004133.h.tsv.gz",
                 n_cases = 12366, n_controls = 33609)
)


# mediators
MEDIATORS <- list(
    SBP = modifyList(GWAS_CATALOG, list(
        label = "Systolic blood pressure",
        file = file.path(dataset_dir, "systolic_bp_GCST90310294.h.tsv.gz"),
        n_col = "n",
        # mmHg -> SD
        standardise_sd = TRUE
    )),
    ApoB = modifyList(GWAS_CATALOG, list(
        label = "Apolipoprotein B",
        file = file.path(dataset_dir, "apolipoprotein_b_GCST90497142.h.tsv.gz"),
        n_col = "n",
        # p_value is just "0.0" at APOE, this column keeps the actual size for clumping
        nlog10_col = "neg_log_10_p_value"
    )),
    # lifted to GRCh38 with datasets/liftover_t2d_to_grch38.R
    T2D = modifyList(T2DGGI, list(
        label = "Type 2 diabetes",
        file = file.path(dataset_dir, "type2_diabetes_T2DGGI_EUR_GRCh38.tsv.gz"),
        n_col = "Neff"
    ))
)

# TwoSampleMR default clumping, against INTERVAL
MED_CLUMP_P <- 5e-8
MED_CLUMP_R2 <- 0.001
MED_CLUMP_KB <- 10000


# cardiometabolic traits (Fig 3C)
# blood pressure in mmHg -> SD, T2D betas are log odds
# GSCAN smoking initiation is case/control but its betas are in SD units
# GSCAN Ns exclude 23andMe (Liu et al. 2019, Suppl. Table 7)
trait_entry <- function(template, label, group, n_label, file, ...) {
    modifyList(template, list(label = label, group = group, n_label = n_label,
                              file = file.path(dataset_dir, file), ...))
}
CARDIOMETABOLIC <- list(
    Triglycerides = trait_entry(GWAS_CATALOG, "Triglycerides", "Lipids", "1,320,016",
                          "triglycerides_GCST90239664.h.tsv.gz"),
    Lpa = trait_entry(GWAS_CATALOG, "Lp(a)", "Lipids", "394,642",
                "lipoprotein_a_GCST90468081.h.tsv.gz"),
    ApoB = trait_entry(GWAS_CATALOG, "ApoB", "Lipids", "434,646",
                 "apolipoprotein_b_GCST90497142.h.tsv.gz"),
    non_HDL_C = trait_entry(GWAS_CATALOG, "Non-HDL-C", "Lipids", "1,320,016",
                      "non_hdl_cholesterol_GCST90239670.h.tsv.gz"),
    LDL_C = trait_entry(GWAS_CATALOG, "LDL-C", "Lipids", "1,320,016",
                  "ldl_cholesterol_GCST90239658.h.tsv.gz"),
    DBP = trait_entry(GWAS_CATALOG, "DBP", "Blood pressure", "1,028,980",
                "diastolic_bp_GCST90310295.h.tsv.gz", n_col = "n", standardise_sd = TRUE),
    SBP = trait_entry(GWAS_CATALOG, "SBP", "Blood pressure", "1,028,980",
                "systolic_bp_GCST90310294.h.tsv.gz", n_col = "n", standardise_sd = TRUE),
    Alcohol = trait_entry(GSCAN, "Alcohol", "Lifestyle", "537,349",
                    "alcohol_consumption_GCST007461.NLRP3region.tsv.gz"),
    Smoking = trait_entry(GSCAN, "Smoking", "Lifestyle", "311,629 / 321,173",
                    "smoking_initiation_GCST007474.NLRP3region.tsv.gz"),
    BMI = trait_entry(GSCAN, "BMI", "Metabolic", "806,834", "bmi_GCST009004.NLRP3region.tsv.gz",
                chr_col = "CHR", ea_col = "Tested_Allele", oa_col = "Other_Allele",
                p_col = "P", sep = " "),
    T2D = trait_entry(T2DGGI, "T2D", "Metabolic", "242,283 / 1,569,734",
                "type2_diabetes_T2DGGI_EUR_GRCh38.tsv.gz", build = "GRCh38", log_odds = TRUE)
)


# CAD studies for the meta-analysis
CAD_STUDIES <- list(
    aragam = modifyList(GWAS_CATALOG, list(
        label = "Aragam et al.",
        file = file.path(dataset_dir, "coronary_artery_disease_aragam_GCST90132314.h.tsv.gz")
    )),
    # no SE, taken from the OR CI
    mvp = modifyList(GWAS_CATALOG_OR, list(
        label = "MVP",
        file = file.path(dataset_dir, "coronary_atherosclerosis_mvp_GCST90475936.h.tsv.gz")
    )),
    finngen = list(
        label = "FinnGen",
        file = file.path(dataset_dir, "coronary_atherosclerosis_finngen_R12.gz"),
        chr_col = "#chrom",
        pos_col = "pos",
        ea_col = "alt",
        oa_col = "ref",
        effect_col = "beta",
        se_col = "sebeta",
        eaf_col = "af_alt",
        p_col = "pval"
    ),
    # All of Us isn't available genome-wide, a collaborator pulled out NLRP3 +/-1 Mb
    # and the mediator instruments for us (SAIGE, phecode CV_404.2, EUR)
    # MarkerID is chr:pos_Allele1/Allele2, BETA and AF_Allele2 are for Allele2
    allofus = list(
        label = "All of Us",
        file = file.path(dataset_dir, "coronary_atherosclerosis_allofus_CV_404_2.instruments.tsv"),
        chr_col = "CHR",
        pos_col = "POS",
        marker_col = "MarkerID",
        effect_col = "BETA",
        se_col = "SE",
        eaf_col = "AF_Allele2",
        p_col = "Pvalue"
    )
)


# analysis constants
PP_THRESHOLD <- 0.8 # HyPrColoc posterior probability
HYPR_PRIOR_1 <- 1e-4 # variant-specific prior, trait association
HYPR_PRIOR_C <- 0.02 # conditional colocalisation prior
