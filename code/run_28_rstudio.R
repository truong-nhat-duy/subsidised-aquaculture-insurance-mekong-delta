# Chay trong RStudio: mo tep nay, bam Source (Ctrl+Shift+S). Mat ~5 phut (che do FAST); dat .R28_FAST <- FALSE de chay day du ~25 phut.
# Ket xuat 5 hinh PDF bang cairo_pdf, font Times New Roman, roi chep vao ...\01. Latex\figures.
LATEX <- "D:/01. HOC HANH/09. FULLBRIGHT/09. MON HOC/17. Environmental Sustainability Policy/10. De cuong luan van_Bao hiem nong nghiep/02. de cuong nghien cuu_new/01. Latex"
RP    <- "D:/01. HOC HANH/09. FULLBRIGHT/09. MON HOC/17. Environmental Sustainability Policy/10. De cuong luan van_Bao hiem nong nghiep/03. data & code/04. R_pipeline"
need <- c("data.table", "sandwich", "ggplot2", "quadprog", "patchwork", "showtext", "sysfonts")
miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(miss)) install.packages(miss)
stopifnot(dir.exists(RP), file.exists(file.path(LATEX, "code", "28_robustness_paper.R")))
OUT <- file.path(LATEX, "robust_out")
setwd(RP)                                   # du lieu goc chi duoc doc, khong ghi
.R28_ARGS <- c(RP, OUT)
.R28_FAST <- TRUE      # TRUE: bo qua hoan vi cap ho (chay ~5 phut, hinh khong doi); FALSE: chay day du ~25 phut
Sys.setenv(R28_RP = RP, R28_OUT = OUT)
source(file.path(LATEX, "code", "28_robustness_paper.R"), encoding = "UTF-8", echo = FALSE)
figs <- list.files(file.path(OUT, "figures"), pattern = "^R[A-D]_.*[.]pdf$", full.names = TRUE)
ok <- file.copy(figs, file.path(LATEX, "figures"), overwrite = TRUE)
cat("\nDa chep", sum(ok), "hinh vao", file.path(LATEX, "figures"), "\n")
