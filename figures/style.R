# shared look for all the R figures, same as style.py

FONT <- "Open Sans" # bold
FONT_BODY <- "Open Sans Semibold" # body text

# also set on the device (plotmath)
cairo_pdf_font <- function(filename, ...) cairo_pdf(filename, ..., family = FONT_BODY)

figures_dir <- here::here("figures_out")
dir.create(figures_dir, showWarnings = FALSE)

# LD colours: no r2, five r2 bins, index variant
LD_SCHEME <- c(
    "grey62", # NA - not in the LD panel
    "#1B36C4", # 0.0 - 0.2
    "#00BEEF", # 0.2 - 0.4
    "#00B81F", # 0.4 - 0.6
    "#FF8A00", # 0.6 - 0.8
    "#EB0000", # 0.8 - 1.0
    "#9A00C8" # index SNP
)
