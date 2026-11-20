# which GO processes are over-represented among the proteome MR hits

source(here::here("config.R"))
source(here::here("helpers.R"))
suppressPackageStartupMessages({
    library(clusterProfiler)
    library(msigdbr)
})

out_dir <- step_dir("10_proteome_ora")

mr <- fread(file.path(results_dir, "09_proteome_mr", "ukb_ppp_proteome_mr_results.tsv"),
            data.table = FALSE)
sig <- filter(mr, ivw_pval_bonf < 0.05)
raised <- sig$gene_name[sig$ivw_beta > 0]
lowered <- sig$gene_name[sig$ivw_beta < 0]
background <- unique(mr$gene_name)

# msigdbr 25.1.1 (MSigDB 2025.1.Hs)
go_bp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
    select(gs_name, gene_symbol)

run_ora <- function(genes, term2gene, direction, collection) {
    enricher(gene = genes, pvalueCutoff = 1, qvalueCutoff = 1,
             universe = background, TERM2GENE = term2gene)@result %>%
        as_tibble() %>%
        mutate(direction = direction, collection = collection, .before = 1)
}
ora <- bind_rows(
    run_ora(raised, go_bp, "raised_by_inhibition", "GO:BP"),
    run_ora(lowered, go_bp, "lowered_by_inhibition", "GO:BP")
)

print(ora %>%
          filter(p.adjust < 0.05) %>%
          select(direction, Description, GeneRatio, p.adjust),
      n = 50)

write_tsv(filter(ora, p.adjust < 0.5), file.path(out_dir, "ukb_ppp_proteome_ora.tsv"))
write_tsv(
    sig %>%
        mutate(direction = if_else(ivw_beta > 0, "raised_by_inhibition", "lowered_by_inhibition")) %>%
        select(gene_name, direction, n_snps, ivw_beta, ivw_se, ivw_pval, ivw_pval_bonf) %>%
        arrange(ivw_pval),
    file.path(out_dir, "ukb_ppp_proteome_significant.tsv")
)
