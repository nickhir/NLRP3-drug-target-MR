# how much of the CAD effect goes through SBP, ApoB and T2D

source(here::here("config.R"))
source(here::here("helpers.R"))

out_dir <- step_dir("07c_mediation")
med_dir <- file.path(results_dir, "07a_mediator_instruments")
set.seed(20260907)
N_DRAWS <- 100000
MEDS <- names(MEDIATORS)

ex <- load_instruments()
R <- interval_ld_matrix(ex$SNP)
bx <- ex$beta_exposure

# correlated IVW with multiplicative random effects
correlated_ivw <- function(by, se_y) {
    S <- diag(se_y, nrow = length(se_y))
    Oi <- solve(S %*% R %*% S)
    info <- as.numeric(t(bx) %*% Oi %*% bx)
    est <- as.numeric(t(bx) %*% Oi %*% by) / info
    resid <- by - est * bx
    sigma2 <- as.numeric(t(resid) %*% Oi %*% resid) / (length(bx) - 1)
    list(est = est, var = max(sigma2, 1) / info, info = info, Sinv = solve(S))
}


# total effect
cad <- fread(file.path(results_dir, "05_mr_cad", "cad_meta_studies_per_snp.tsv"),
             data.table = FALSE) %>%
    filter(study == "Meta-analysis")
cad <- cad[match(ex$SNP, cad$SNP), ]
tau_fit <- correlated_ivw(cad$by, cad$byse)
tau <- tau_fit$est
var_tau <- tau_fit$var


# NLRP3 -> mediators
sd_scales <- readRDS(file.path(med_dir, "selection_summary.rds"))$sd_scales
want <- transmute(ex, SNPid = SNP, chrom = chr, pos = pos_hg38)
med_rows <- lapply(MEDS, function(k) {
    d <- lookup_at(MEDIATORS[[k]], want)
    d[match(ex$SNP, d$SNPid), ]
})
names(med_rows) <- MEDS
alpha_fits <- lapply(MEDS, function(k) {
    correlated_ivw(med_rows[[k]]$beta / sd_scales[[k]], med_rows[[k]]$se / sd_scales[[k]])
})
names(alpha_fits) <- MEDS
alpha <- sapply(alpha_fits, `[[`, "est")
var_alpha <- sapply(alpha_fits, `[[`, "var")
info <- sapply(alpha_fits, `[[`, "info")


# mediators -> CAD (MVMR)
design <- fread(file.path(med_dir, "mvmr_design.tsv"), data.table = FALSE)
keep <- readRDS(file.path(med_dir, "mvmr_instrument_sets.rds"))[[paste(MEDS, collapse = "+")]]
mvmr <- fit_mvmr(design, keep, MEDS)
beta <- mvmr$beta
Sigma_bb <- mvmr$V

# conditional F
rho <- readRDS(file.path(results_dir, "07b_sample_overlap", "rho.rds"))[MEDS, MEDS]
d <- design[design$SNPid %in% keep, ]
SEX <- as.matrix(d[paste0("se_sd_", MEDS)])
fmt <- MVMR::format_mvmr(BXGs = as.matrix(d[paste0("beta_sd_", MEDS)]), BYG = d$beta_cad,
                         seBXGs = SEX, seBYG = d$se_cad, RSID = d$SNPid)
cond_F <- as.numeric(MVMR::strength_mvmr(
    r_input = fmt, gencov = MVMR::phenocov_mvmr(pcor = rho, seBXGs = SEX)))


# covariance of the alphas
# the alphas are correlated through sample overlap (rho from 07b)
Ri <- solve(R)
K <- length(MEDS)
iA <- seq_len(K)
iB <- K + iA

build_sigma <- function(rho_mat) {
    S_a <- diag(var_alpha, nrow = K)
    dimnames(S_a) <- list(MEDS, MEDS)
    for (i in seq_len(K)) for (j in seq_len(K)) {
        if (j <= i) next
        C <- as.numeric(t(bx) %*% alpha_fits[[i]]$Sinv %*% Ri %*% alpha_fits[[j]]$Sinv %*% bx)
        S_a[i, j] <- S_a[j, i] <- rho_mat[i, j] * C / (info[i] * info[j])
    }
    nm <- c(paste0("a_", MEDS), paste0("b_", MEDS))
    S <- matrix(0, 2 * K, 2 * K, dimnames = list(nm, nm))
    S[iA, iA] <- S_a
    S[iB, iB] <- Sigma_bb
    S
}
Sigma_theta <- build_sigma(rho)


# indirect effects and proportion mediated
per_mediator <- tibble(
    mediator = MEDS,
    label = sapply(MEDIATORS[MEDS], `[[`, "label"),
    alpha = alpha, se_alpha = sqrt(var_alpha),
    beta = beta, se_beta = sqrt(diag(Sigma_bb)),
    indirect = alpha * beta,
    se_indirect = sqrt(beta^2 * var_alpha + alpha^2 * diag(Sigma_bb))
) %>%
    mutate(or = exp(indirect),
           or_lower = exp(indirect - 1.96 * se_indirect),
           or_upper = exp(indirect + 1.96 * se_indirect),
           p = 2 * pnorm(-abs(indirect / se_indirect)),
           pm_alone = indirect / tau)

# delta method SE
var_joint <- function(Sigma) {
    grad <- c(beta, alpha)
    as.numeric(t(grad) %*% Sigma %*% grad)
}
ie <- sum(alpha * beta)
var_ie <- var_joint(Sigma_theta)
se_ie <- sqrt(var_ie)
direct <- tau - ie
se_dir <- sqrt(var_tau + var_ie)

# CI for the proportion mediated by Monte Carlo
draw_theta <- MASS::mvrnorm(N_DRAWS, mu = c(alpha, beta), Sigma = Sigma_theta)
draw_tau <- rnorm(N_DRAWS, tau, sqrt(var_tau))
draw_pm <- rowSums(draw_theta[, iA, drop = FALSE] * draw_theta[, iB, drop = FALSE]) / draw_tau
pm_ci <- quantile(draw_pm, c(0.025, 0.975), names = FALSE)

# SE assuming rho = 0 and rho = 1
bound <- sapply(list(zero = diag(1, K), est = rho, one = matrix(1, K, K)), function(m) {
    sqrt(var_joint(build_sigma(m)))
})

# difference method as a check
med_at_cis <- sapply(MEDS, function(k) med_rows[[k]]$beta / sd_scales[[k]])
direct2 <- correlated_ivw(cad$by - as.numeric(med_at_cis %*% beta), cad$byse)$est


summary_tbl <- tibble(
    quantity = c("Total effect", "Joint indirect", "Direct effect"),
    estimate = c(tau, ie, direct),
    se = c(sqrt(var_tau), se_ie, se_dir)
) %>%
    mutate(or = exp(estimate),
           or_lower = exp(estimate - 1.96 * se),
           or_upper = exp(estimate + 1.96 * se),
           p = 2 * pnorm(-abs(estimate / se)))

pm_tbl <- tibble(proportion_mediated = ie / tau, ci_lower = pm_ci[1], ci_upper = pm_ci[2],
                 method = "Monte Carlo", n_draws = N_DRAWS)

print(as.data.frame(summary_tbl), row.names = FALSE, digits = 3)
print(as.data.frame(pm_tbl), row.names = FALSE, digits = 3)

fwrite(summary_tbl, file.path(out_dir, "mediation_summary.tsv"), sep = "\t")
fwrite(per_mediator, file.path(out_dir, "mediation_per_mediator.tsv"), sep = "\t")
fwrite(pm_tbl, file.path(out_dir, "proportion_mediated.tsv"), sep = "\t")
saveRDS(list(tau = tau, var_tau = var_tau, alpha = alpha, var_alpha = var_alpha,
             beta = beta, Sigma_bb = Sigma_bb, rho = rho, Sigma_theta = Sigma_theta,
             ie = ie, var_ie = var_ie, direct = direct, se_direct = se_dir,
             pm = ie / tau, pm_ci = pm_ci, pm_draws = draw_pm,
             cond_F = cond_F, Q = mvmr$Q, n_mvmr = mvmr$n,
             mvmr_method = "IVW after pooled LD clumping, multiplicative random effects",
             mvmr_snps = mvmr$snps,
             se_bounds = bound, direct_difference_method = direct2,
             exposure_direction = EXPOSURE_DIRECTION),
        file.path(out_dir, "mediation.rds"))
