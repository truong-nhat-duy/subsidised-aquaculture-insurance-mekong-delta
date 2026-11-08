# =============================================================================
# 28_robustness_paper.R (18/10/2026) -- bo kiem dinh do vung cho paper_vi.tex
# Chay: Rscript 28_robustness_paper.R <R_pipeline> <thu muc dau ra>
# Doc  : data_proc/{panel_aqua,hazard_tc_py,hazard_emdat_py,vhlss_model}.rds|csv (chi doc, khong sua)
# Ghi  : <dau ra>/tables/R*.csv, <dau ra>/figures/R*.pdf, <dau ra>/R28.log
# Nhom : A tinh (TWFE bao -> san luong), B chinh sach (DiD/SDID), C ho (hazard -> san luong),
#        D chi so luoi (AUC). Moi nhom: tai lap ket qua goc truoc, roi moi mo rong.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(sandwich); library(quadprog)})
args <- if (exists(".R28_ARGS")) .R28_ARGS else if (nzchar(Sys.getenv("R28_OUT"))) c(Sys.getenv("R28_RP"), Sys.getenv("R28_OUT")) else commandArgs(TRUE); RP <- if (length(args) >= 1) args[1] else "."; OUT <- if (length(args) >= 2) args[2] else "robust_out"
dir.create(file.path(OUT, "tables"), recursive = TRUE, showWarnings = FALSE); dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, "R28.log"); cat("", file = LOG)
lg <- function(...) { s <- paste0(...); cat(s, "\n"); cat(s, "\n", file = LOG, append = TRUE) }
set.seed(20261018)
BLUE <- "#00196E"; GOLD <- "#C99700"; GREY <- "#8C8C8C"; RED <- "#C0392B"
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a
fam <- (function() {   # font hoc thuat: Times New Roman nhung bang showtext (Windows); TeX Gyre Termes (Linux) neu co
  if (.Platform$OS.type == "windows") {
    wf <- file.path(Sys.getenv("WINDIR", "C:/Windows"), "Fonts")
    if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE) && file.exists(file.path(wf, "times.ttf"))) {
      sysfonts::font_add("TNR", regular = file.path(wf, "times.ttf"), bold = file.path(wf, "timesbd.ttf"),
                         italic = file.path(wf, "timesi.ttf"), bolditalic = file.path(wf, "timesbi.ttf"))
      showtext::showtext_auto(); return("TNR") }
    return("Times New Roman") }
  fl <- tryCatch(system("fc-list | grep -i 'TeX Gyre Termes'", intern = TRUE, ignore.stderr = TRUE), error = function(e) character(0))
  if (length(fl) && nzchar(fl[1])) "TeX Gyre Termes" else "serif" })()
th <- function() theme_minimal(base_size = 9, base_family = fam) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
  strip.text = element_text(face = "bold", hjust = 0), plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
savef <- function(p, name, w, h) ggsave(file.path(OUT, "figures", paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
wt <- function(x, name) fwrite(x, file.path(OUT, "tables", paste0(name, ".csv")))
fmt <- function(x, d = 3) formatC(x, digits = d, format = "f")

# ---- helpers: hai chieu FE bang chieu xen ke (MAP), SE cum CR1 ---------------------------------
demean2 <- function(M, g1, g2 = NULL, tol = 1e-10, maxit = 500) {
  M <- as.matrix(M); n1 <- tabulate(g1); n2 <- if (!is.null(g2)) tabulate(g2)
  for (it in seq_len(maxit)) { old <- M
    M <- M - (rowsum(M, g1) / n1)[g1, , drop = FALSE]
    if (!is.null(g2)) M <- M - (rowsum(M, g2) / n2)[g2, , drop = FALSE]
    if (is.null(g2) || max(abs(M - old)) < tol) break }
  M }
fe_ols <- function(y, X, g1, g2 = NULL, cl, kfe2 = NULL, ret_resid = FALSE) {
  X <- as.matrix(X); Z <- demean2(cbind(y, X), g1, g2); yt <- Z[, 1]; Xt <- Z[, -1, drop = FALSE]
  keep <- apply(Xt, 2, function(v) any(abs(v) > 1e-12)); Xt <- Xt[, keep, drop = FALSE]
  XtX <- crossprod(Xt); b <- solve(XtX, crossprod(Xt, yt)); e <- as.vector(yt - Xt %*% b)
  G <- length(unique(cl)); N <- length(y); K <- ncol(Xt) + (if (is.null(kfe2)) 0 else kfe2)
  sc <- rowsum(Xt * e, cl); meat <- crossprod(sc); XtXi <- solve(XtX)
  V <- XtXi %*% meat %*% XtXi * (G / (G - 1)) * ((N - 1) / (N - K))
  out <- data.frame(term = colnames(Xt), est = as.vector(b), se = sqrt(diag(V)), G = G, N = N)
  out$t <- out$est / out$se; out$p <- 2 * pt(-abs(out$t), df = G - 1)
  if (ret_resid) attr(out, "resid") <- list(e = e, Xt = Xt, yt = yt); out }
idx <- function(x) as.integer(factor(x))
lm_cl <- function(f, d, cl) { m <- lm(f, d); V <- vcovCL(m, cluster = d[[cl]], type = "HC1"); list(m = m, V = V) }

# =============================================================================
# A. PANEL TINH: bao -> log san luong (tai lap Bang A1 roi mo rong)
# =============================================================================
pan <- fread(file.path(RP, "data_proc/panel_aqua.csv")); tc <- fread(file.path(RP, "data_proc/hazard_tc_py.csv"))
setorder(tc, key, year); tc[, `:=`(tc_di_l1 = shift(tc_di), tc_di_l2 = shift(tc_di, 2), tc_ts_l1 = shift(tc_share_ts)), by = key]
d <- merge(pan, tc, by = c("key", "year"), all.x = TRUE); d <- d[year <= min(2023, max(tc$year)) & (coastal | mekong)]
d[, `:=`(tc_ts_pos = as.integer(tc_share_ts > 0), t_ = year - 2005)]
lg("A. mau tinh: ", uniqueN(d$key), " tinh, ", min(d$year), "-", max(d$year))
twfe_p <- function(dep, rhs, dat, extra = "") {
  dat <- dat[complete.cases(dat[, c(dep, rhs), with = FALSE])]
  dat[, `:=`(g1 = idx(key), g2 = idx(year))]
  f <- as.formula(paste(dep, "~", paste(rhs, collapse = "+"), "+ factor(key) + factor(year)", extra))
  m <- lm(f, dat); V <- vcovCL(m, cluster = dat$key, type = "HC1"); b <- coef(m)[rhs]
  data.frame(term = rhs, est = unname(b), se = sqrt(diag(V)[rhs]), N = nobs(m), G = uniqueN(dat$key), vc = I(replicate(length(rhs), V[rhs, rhs, drop = FALSE], simplify = FALSE)))
}
base_shr <- d[q_shrimp >= 100 & !is.na(q_shrimp)]
A0 <- rbind(cbind(model = "Tong", twfe_p("ln_q", c("tc_di", "tc_di_l1"), d)), cbind(model = "Tom", twfe_p("ln_shrimp", c("tc_di", "tc_di_l1"), base_shr)))
lg("A0 tai lap: tom tc_di ", fmt(A0$est[3]), " (SE ", fmt(A0$se[3]), "), lag ", fmt(A0$est[4]), " (SE ", fmt(A0$se[4]), "); tong ", fmt(A0$est[1]))
lg("KIEM TRA TAI LAP A0: tom ", fmt(A0$est[3]), " vs T5 -0.5308; tong ", fmt(A0$est[1]), " vs T5 -0.1586 -> ", ifelse(abs(A0$est[3] + 0.5308) < 1e-3 && abs(A0$est[1] + 0.1586) < 1e-3, "KHOP", "LECH"))

# A1. leave-one-province-out (tom) -> do on dinh cua he so
prov <- sort(unique(base_shr$key))
loo <- rbindlist(lapply(prov, function(p) { r <- twfe_p("ln_shrimp", c("tc_di", "tc_di_l1"), base_shr[key != p]); r$drop <- p; as.data.table(r)[, !"vc"] }))
loo <- merge(loo, unique(pan[, .(key = key, mekong = mekong)]), by.x = "drop", by.y = "key")
wt(loo, "RA1_leave_one_province")
lg("A1 LOO tom tc_di: min ", fmt(min(loo[term == "tc_di", est])), ", max ", fmt(max(loo[term == "tc_di", est])), "; tong |t| lon nhat ",
   fmt(max(abs(loo[term == "tc_di", est / se]))))

# A2. suy luan hoan vi: hoan vi bien phoi nhiem trong tung tinh qua cac nam (B = 5000)
perm_p <- function(dat, dep, B = 5000) {
  dat <- copy(dat[complete.cases(dat[, c(dep, "tc_di", "tc_di_l1"), with = FALSE])]); setorder(dat, key, year)
  g1 <- idx(dat$key); g2 <- idx(dat$year); y <- dat[[dep]]
  stat <- function(X) { r <- fe_ols(y, X, g1, g2, cl = g1, kfe2 = max(g2) - 1); r }
  obs <- stat(cbind(tc_di = dat$tc_di, tc_di_l1 = dat$tc_di_l1))
  # hoan vi tc_di trong tinh (giu tuong quan tc_di va tc_di_l1: hoan vi chuoi (t, t-1) theo nam)
  yrs_by <- split(seq_len(nrow(dat)), g1)
  draws <- replicate(B, { ord <- unlist(lapply(yrs_by, function(ix) if (length(ix) > 1) sample(ix) else ix), use.names = FALSE)
    X <- cbind(tc_di = dat$tc_di, tc_di_l1 = dat$tc_di_l1); Xp <- X; Xp[unlist(yrs_by, use.names = FALSE), ] <- X[ord, ]
    stat(Xp)$est }, simplify = TRUE)
  list(obs = obs, draws = t(draws)) }
pr_s <- perm_p(base_shr, "ln_shrimp"); pr_t <- perm_p(d, "ln_q")
ri_tab <- rbind(data.frame(model = "Tom", term = c("tc_di", "tc_di_l1"), est = pr_s$obs$est, p_cl = pr_s$obs$p,
                           p_perm = c(mean(abs(pr_s$draws[, 1]) >= abs(pr_s$obs$est[1])), mean(abs(pr_s$draws[, 2]) >= abs(pr_s$obs$est[2])))),
                data.frame(model = "Tong", term = c("tc_di", "tc_di_l1"), est = pr_t$obs$est, p_cl = pr_t$obs$p,
                           p_perm = c(mean(abs(pr_t$draws[, 1]) >= abs(pr_t$obs$est[1])), mean(abs(pr_t$draws[, 2]) >= abs(pr_t$obs$est[2])))))
wt(ri_tab, "RA2_permutation"); lg("A2 hoan vi (B=5000): "); for (i in seq_len(nrow(ri_tab))) lg("   ", ri_tab$model[i], " ", ri_tab$term[i], " est ", fmt(ri_tab$est[i]), " p_CR1 ", fmt(ri_tab$p_cl[i]), " p_perm ", fmt(ri_tab$p_perm[i]))

# A3. cac dac ta thay the (tom)
alt <- rbindlist(list(
  cbind(spec = "Co so", twfe_p("ln_shrimp", c("tc_di", "tc_di_l1"), base_shr)),
  cbind(spec = "Them tre t-2", twfe_p("ln_shrimp", c("tc_di", "tc_di_l1", "tc_di_l2"), base_shr)),
  cbind(spec = "Xu the tuyen tinh rieng", twfe_p("ln_shrimp", c("tc_di", "tc_di_l1"), base_shr, extra = "+ factor(key):t_")),
  cbind(spec = "Bien nhi phan (co gio >=34 kt)", twfe_p("ln_shrimp", c("tc_ts_pos"), base_shr)),
  cbind(spec = "Ngo ai tru Mekong", twfe_p("ln_shrimp", c("tc_di", "tc_di_l1"), base_shr[mekong == FALSE])),
  cbind(spec = "Bo tinh co Q tom < 1.000 t TB", twfe_p("ln_shrimp", c("tc_di", "tc_di_l1"), base_shr[key %in% base_shr[, .(m = mean(q_shrimp)), by = key][m >= 1000, key]])),
  cbind(spec = "Bo 2020-2021", twfe_p("ln_shrimp", c("tc_di", "tc_di_l1"), base_shr[!year %in% 2020:2021]))), fill = TRUE)
alt[, vc := NULL]; alt[, p := 2 * pt(-abs(est / se), df = G - 1)]
# tong tac dong tich luy t + t-1 (co so)
V0 <- A0$vc[[3]]; cum <- sum(A0$est[3:4]); cum_se <- sqrt(sum(V0)); alt <- rbind(alt, data.table(spec = "Tich luy t + t-1 (co so)", term = "tong", est = cum, se = cum_se, N = A0$N[3], G = A0$G[3], p = 2 * pt(-abs(cum / cum_se), df = A0$G[3] - 1)), fill = TRUE)
wt(alt, "RA3_alt_specs_province"); lg("A3 dac ta thay the: "); for (i in seq_len(nrow(alt))) lg("   ", alt$spec[i], " | ", alt$term[i], " ", fmt(alt$est[i]), " (", fmt(alt$se[i]), ") p=", fmt(alt$p[i]))

pA <- ggplot(loo[term == "tc_di"], aes(est, reorder(drop, est), colour = mekong)) + geom_vline(xintercept = A0$est[3], colour = GREY, linetype = 2) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), height = 0, linewidth = 0.3, alpha = 0.6) + geom_point(size = 1.2) +
  scale_colour_manual(values = c(`TRUE` = GOLD, `FALSE` = BLUE), labels = c(`TRUE` = "Tỉnh ĐBSCL", `FALSE` = "Tỉnh khác"), name = NULL) +
  labs(x = "Hệ số của chỉ số bão (log sản lượng tôm) khi bỏ lần lượt từng tỉnh; KTC 95% (CR1)", y = NULL,
       title = "(a) Bỏ một tỉnh: hệ số tôm cùng năm") + th() + theme(axis.text.y = element_text(size = 5.5))
pAb <- ggplot(data.frame(b = pr_s$draws[, 1]), aes(b)) + geom_histogram(bins = 50, fill = "grey80", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = pr_s$obs$est[1], colour = RED, linewidth = 0.6) +
  labs(x = "Hệ số khi hoán vị phơi nhiễm trong từng tỉnh (5.000 lần)", y = "Số lần", title = "(b) Suy luận hoán vị, tôm, bão cùng năm") + th()
if (requireNamespace("patchwork", quietly = TRUE)) { savef(pA + pAb + patchwork::plot_layout(widths = c(1.15, 1)), "RA_province_robust", 7.2, 3.9)
} else { savef(pA, "RA_province_robust_a", 4, 4); savef(pAb, "RA_province_robust_b", 4, 4) }

# =============================================================================
# B. CHINH SACH: DiD / SDID (13 tinh ĐBSCL, 5 thi diem) -- tai lap roi mo rong
# =============================================================================
pm <- pan[mekong == TRUE]; pm[, `:=`(D = as.integer(treated), post = as.integer(year >= 2019))]
twfe_att <- function(dd, y, post_year = 2019, trend = FALSE, drop_years = NULL) {
  dd <- dd[!is.na(get(y))]; if (!is.null(drop_years)) dd <- dd[!year %in% drop_years]
  dd[, pst := as.integer(year >= post_year)]
  f <- as.formula(paste(y, "~ I(D*pst) + factor(key) + factor(year)", if (trend) "+ factor(key):I(year-2005)" else ""))
  m <- lm(f, dd); V <- vcovCL(m, cluster = dd$key, type = "HC1"); c(est = unname(coef(m)[2]), se = sqrt(V[2, 2]), N = nobs(m)) }
OUTS <- c(ln_q = "Tổng sản lượng", ln_shrimp = "Tôm", ln_fish = "Cá")
keys13 <- sort(unique(pm$key)); tk <- unique(pm[D == 1, key])
okmat <- function(dd, y, shrimp = FALSE) { ok <- dd[, .(ok = all(is.finite(get(y)))), by = key][ok == TRUE, key]; db <- dd[key %in% ok]
  if (shrimp) { big <- db[year < 2019, .(m = mean(q_shrimp)), by = key][m >= 1000 | key %in% tk, key]; db <- db[key %in% big] }; db }
B1 <- rbindlist(lapply(names(OUTS), function(y) {
  db <- okmat(pm, y, y == "ln_shrimp"); r <- twfe_att(copy(db), y)
  r_tr <- twfe_att(copy(db), y, trend = TRUE); r_nc <- twfe_att(copy(db), y, drop_years = 2011:2013); r_20 <- twfe_att(copy(db), y, post_year = 2020)
  data.table(outcome = OUTS[y], n_tr = uniqueN(db[D == 1, key]), n_co = uniqueN(db[D == 0, key]),
             did = r["est"], did_se = r["se"], did_trend = r_tr["est"], did_trend_se = r_tr["se"],
             did_no1113 = r_nc["est"], did_no1113_se = r_nc["se"], did_post2020 = r_20["est"], did_post2020_se = r_20["se"]) }))
wt(B1, "RB1_did_variants"); lg("B1 DiD (13 tinh) -> kiem tra tai lap:"); print(B1)
lg("KIEM TRA TAI LAP B1 DiD tong ", fmt(B1$did[1]), " vs T7b 0.00514")

# B2. placebo theo thoi gian: gan nam ap dung gia 2010..2017, du lieu <= 2018
pb <- rbindlist(lapply(names(OUTS), function(y) {
  db <- okmat(pm[year <= 2018], y, y == "ln_shrimp"); rbindlist(lapply(2010:2017, function(py) { r <- twfe_att(copy(db), y, post_year = py)
    data.table(outcome = OUTS[y], placebo_year = py, est = r["est"], se = r["se"]) })) }))
act <- B1[, .(outcome, est = did, se = did_se)]
wt(pb, "RB2_placebo_in_time")
lg("B2 placebo thoi gian: |est| lon nhat theo ket qua: "); print(pb[, .(max_abs = max(abs(est)), mean = mean(est)), by = outcome])

# B3. placebo khong gian: hoan vi tat ca C(13,5) = 1287 cach chon 5 tinh (Tong, Ca; tom dung mau con rieng)
ri <- rbindlist(lapply(names(OUTS), function(y) { db <- okmat(pm, y, y == "ln_shrimp"); ks <- sort(unique(db$key)); n_tr <- sum(ks %in% tk)
  cm <- combn(ks, n_tr); est <- apply(cm, 2, function(t5) { dd <- copy(db); dd[, D := as.integer(key %in% t5)]; twfe_att(dd, y)["est"] })
  obs <- twfe_att(copy(db), y)["est"]; data.table(outcome = OUTS[y], est = est, obs = obs, p_ri = mean(abs(est) >= abs(obs) - 1e-12), n_assign = length(est)) }))
wt(unique(ri[, .(outcome, obs, p_ri, n_assign)]), "RB3_randomisation_inference"); lg("B3 RI: "); print(unique(ri[, .(outcome, obs, p_ri, n_assign)]))

# B4. SDID (cai dat cua nhom: quadprog) -- tre-dich nam ap dung, bo tung tinh thi diem, tap tinh cho thay the
sdid <- function(Y, n_tr, t_pre) { N <- nrow(Y); TT <- ncol(Y); n_co <- N - n_tr; t_post <- TT - t_pre
  Yco <- Y[1:n_co, , drop = FALSE]; Ytr <- Y[(n_co + 1):N, , drop = FALSE]
  dif <- t(apply(Yco[, 1:t_pre, drop = FALSE], 1, diff)); zeta <- (n_tr * t_post)^(1/4) * sd(as.vector(dif))
  sw <- function(A, b, pen) { k <- ncol(A); X <- cbind(1, A); Dm <- crossprod(X) + diag(c(1e-6 * (1 + sum(X^2)), rep(pen + 1e-8, k))); dv <- crossprod(X, b)
    Am <- cbind(c(0, rep(1, k)), rbind(0, diag(k))); bv <- c(1, rep(0, k)); sol <- solve.QP(Dm, dv, Am, bv, meq = 1)$solution; pmax(sol[-1], 0) / sum(pmax(sol[-1], 0)) }
  w <- sw(t(Yco[, 1:t_pre, drop = FALSE]), colMeans(Ytr[, 1:t_pre, drop = FALSE]), zeta^2 * t_pre)
  lam <- sw(Yco[, 1:t_pre, drop = FALSE], rowMeans(Yco[, (t_pre + 1):TT, drop = FALSE]), 1e-6 * sd(as.vector(Yco[, 1:t_pre]))^2 * n_co)
  tr_post <- mean(Ytr[, (t_pre + 1):TT]); tr_pre <- sum(lam * colMeans(Ytr[, 1:t_pre, drop = FALSE]))
  co_post <- sum(w * rowMeans(Yco[, (t_pre + 1):TT, drop = FALSE])); co_pre <- sum(w * (Yco[, 1:t_pre, drop = FALSE] %*% lam))
  list(tau = (tr_post - tr_pre) - (co_post - co_pre), w = w) }
to_mat <- function(dd, y, tkeys) { w <- dcast(dd[, c("key", "year", y), with = FALSE], key ~ year, value.var = y); w <- w[order(key %in% tkeys)]
  M <- as.matrix(w[, -1]); rownames(M) <- w$key; M }
sd_run <- function(db, y, tkeys, y0 = 2019, B = 300, se = TRUE) { M <- to_mat(db, y, tkeys); t_pre <- sum(as.integer(colnames(M)) < y0); n_tr <- sum(rownames(M) %in% tkeys); n_co <- nrow(M) - n_tr
  fit <- sdid(M, n_tr, t_pre); s <- NA_real_
  if (se && n_co > n_tr) { pl <- replicate(B, { ix <- sample(seq_len(n_co), min(n_tr, n_co - 1)); sdid(M[c(setdiff(seq_len(n_co), ix), ix), ], length(ix), t_pre)$tau }); s <- sd(pl) }
  c(tau = fit$tau, se = s) }
B4 <- rbindlist(lapply(names(OUTS), function(y) { db <- okmat(pm, y, y == "ln_shrimp"); base <- sd_run(db, y, tk, 2019)
  loo5 <- rbindlist(lapply(tk, function(p) { r <- sd_run(db[key != p], y, setdiff(tk, p), 2019, se = FALSE); data.table(drop = p, tau = r["tau"]) }))
  y20 <- sd_run(db, y, tk, 2020, se = FALSE); y18 <- sd_run(db, y, tk, 2018, se = FALSE)
  data.table(outcome = OUTS[y], sdid = base["tau"], sdid_se = base["se"], sdid_post2020 = y20["tau"], sdid_post2018 = y18["tau"],
             loo_min = min(loo5$tau), loo_max = max(loo5$tau)) }))
set.seed(20261018)
# tap tinh cho rong: 13 tinh ĐBSCL + ven bien (mau rong) cho tong san luong va ca
pw <- pan[(mekong | coastal)]; pw[, D := as.integer(treated)]
B4w <- rbindlist(lapply(c("ln_q", "ln_fish"), function(y) { db <- okmat(pw, y); r <- sd_run(db, y, tk, 2019, B = 300); data.table(outcome = paste0(OUTS[y], " (tap cho rong)"), sdid = r["tau"], sdid_se = r["se"]) }), fill = TRUE)
B4 <- rbind(B4, B4w, fill = TRUE); wt(B4, "RB4_sdid_variants"); lg("B4 SDID: "); print(B4)

pB1 <- {x <- rbind(act[, .(outcome, src = "Thực tế (2019)", year = 2019, est, se)], pb[, .(outcome, src = "Placebo", year = placebo_year, est, se)])
  ggplot(x, aes(factor(year), est, colour = src)) + geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se), width = 0, linewidth = 0.4, position = position_dodge(0.4)) + geom_point(size = 1.6) +
  facet_wrap(~outcome, nrow = 1, scales = "free_y") + scale_colour_manual(values = c(`Thực tế (2019)` = RED, Placebo = BLUE), name = NULL) +
  labs(x = "Năm áp dụng (placebo: dữ liệu đến 2018)", y = "Hệ số DiD (log); KTC 95%", title = "(a) DiD thực tế và placebo theo thời gian") + th() }
pB2 <- ggplot(ri, aes(est)) + geom_histogram(bins = 40, fill = "grey80", colour = "white", linewidth = 0.15) + geom_vline(aes(xintercept = obs), colour = RED, linewidth = 0.6) +
  facet_wrap(~outcome, nrow = 1, scales = "free") + labs(x = "Ước lượng DiD khi gán ngẫu nhiên 5 tỉnh thí điểm (mọi cách gán)", y = "Số cách gán", title = "(b) Suy luận hoán vị theo không gian") + th()
if (requireNamespace("patchwork", quietly = TRUE)) savef(pB1 / pB2 + patchwork::plot_layout(heights = c(1, 0.9)), "RB_policy_robust", 7.2, 5.6) else { savef(pB1, "RB_policy_a", 7, 3); savef(pB2, "RB_policy_b", 7, 3) }

# =============================================================================
# C. HO: hazard -> san luong ho (tai lap Bang 5 roi mo rong)
# =============================================================================
s <- as.data.table(readRDS(file.path(RP, "data_proc/vhlss_model.rds")))
s[, region_year := paste(region, year, sep = "_")]
ctrl <- intersect(c("head_age", "head_male", "head_grade_f", "grade_na", "n_mem", "dep_ratio"), names(s))
sf <- s[aq_farm %in% 1]
hh_fit <- function(dep, rhs, dat, fe2 = "prov_year", cl = "commune_id", fe1 = "commune_id") {
  dat <- dat[complete.cases(dat[, c(dep, rhs, fe1, fe2, cl), with = FALSE])]; g1 <- idx(dat[[fe1]]); g2 <- idx(dat[[fe2]]); g <- idx(dat[[cl]])
  r <- fe_ols(dat[[dep]], as.matrix(dat[, rhs, with = FALSE]), g1, g2, cl = g, kfe2 = max(g2) - 1); r$n_cl <- uniqueN(dat[[cl]]); r }
C0 <- hh_fit("y_farm", c("H_wet_z", "H_dry_z", ctrl), sf)
lg("C0 tai lap Bang 5 (y_farm ~ wet+dry | commune + prov_year): H_wet ", fmt(C0$est[1]), " SE ", fmt(C0$se[1]), " N ", C0$N[1], "; H_dry ", fmt(C0$est[2]), " SE ", fmt(C0$se[2]))
lg("KIEM TRA TAI LAP C0 H_wet ", fmt(C0$est[1]), " vs bai 0.055; N=", C0$N[1], " vs 10442")
sh <- s[aq_shrimp %in% 1]
C0s <- hh_fit("y_shrimp", c("H_wet_z", "H_dry_z", ctrl), sh); lg("   tom: H_wet ", fmt(C0s$est[1]), " (", fmt(C0s$se[1]), ") N ", C0s$N[1])

# C1. duong cong dac ta (specification curve)
specs <- CJ(dep = c("y_farm", "y_shrimp", "y_out"), haz = c("z", "pca"), samp = c("Tat ca", "ĐBSCL"), fe2 = c("prov_year", "region_year", "year"), ctl = c(TRUE, FALSE))
res <- rbindlist(lapply(seq_len(nrow(specs)), function(i) { sp <- specs[i]
  dat <- if (sp$dep == "y_shrimp") sh else if (sp$dep == "y_out") s[aq_any %in% 1 | aq_farm %in% 1] else sf
  if (sp$samp == "ĐBSCL") dat <- dat[mekong %in% TRUE]
  hz <- if (sp$haz == "z") c("H_wet_z", "H_dry_z") else c("H_wet_pca", "H_dry_pca")
  r <- tryCatch(hh_fit(sp$dep, c(hz, if (sp$ctl) ctrl), dat, fe2 = sp$fe2), error = function(e) NULL); if (is.null(r)) return(NULL)
  r <- r[r$term %in% hz, ]; cbind(sp[rep(1, nrow(r))], r[, c("term", "est", "se", "p", "N")]) }))
res[, `:=`(block = ifelse(grepl("wet", term), "Ướt (lũ, bão, mưa cực trị)", "Khô–nóng (hạn, nhiệt)"), lo = est - 1.96 * se, hi = est + 1.96 * se, sig = p < 0.05)]
wt(res, "RC1_spec_curve"); n_sig <- sum(res$sig); lg("C1 duong cong dac ta: ", nrow(res), " he so; ", n_sig, " co p<0,05 (", fmt(100 * n_sig / nrow(res), 1), "%), trung vi he so ", fmt(median(res$est)),
   "; khoang ", fmt(min(res$est)), " den ", fmt(max(res$est)), "; |t| lon nhat ", fmt(max(abs(res$est / res$se)), 2))
res[, rk := rank(est, ties.method = "first"), by = block]
pC1 <- ggplot(res, aes(rk, est, colour = sig)) + geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.25, alpha = 0.5) + geom_point(size = 0.9) + facet_wrap(~block, nrow = 1) +
  scale_colour_manual(values = c(`FALSE` = BLUE, `TRUE` = RED), labels = c(`FALSE` = "p ≥ 0,05", `TRUE` = "p < 0,05"), name = NULL) +
  labs(x = "Đặc tả (xếp theo hệ số)", y = "Hệ số trên 1 độ lệch chuẩn hiểm họa (asinh); KTC 95%", title = "Đường cong đặc tả: mọi tổ hợp kết quả, chỉ số, mẫu, hiệu ứng cố định, biến kiểm soát") + th()
savef(pC1, "RC_spec_curve", 7.2, 3.4)

# C2. cap cum SE khac nhau (huyen, tinh) cho dac ta co so
sf[, dist_id := paste(tinh, huyen)]; sh[, dist_id := paste(tinh, huyen)]
cl_tab <- rbindlist(lapply(list(list("Tong nuoi trong", sf, "y_farm"), list("Tom", sh, "y_shrimp")), function(a) rbindlist(lapply(c("commune_id", "dist_id", "tinh"), function(cc) {
  r <- hh_fit(a[[3]], c("H_wet_z", "H_dry_z", ctrl), a[[2]], cl = cc); r <- r[r$term %in% c("H_wet_z", "H_dry_z"), ]; data.table(sample = a[[1]], cluster = cc, term = r$term, est = r$est, se = r$se, p = r$p, G = r$G) }))))
wt(cl_tab, "RC2_cluster_levels"); lg("C2 cum SE:"); print(cl_tab)

# C3. suy luan hoan vi: hoan vi (H_wet, H_dry) giua cac xa trong tung tinh-nam (B = 2000)
perm_hh <- function(dep, dat, B = 2000) {
  dat <- dat[complete.cases(dat[, c(dep, "H_wet_z", "H_dry_z", ctrl, "commune_id", "prov_year"), with = FALSE])]
  g1 <- idx(dat$commune_id); g2 <- idx(dat$prov_year); y <- dat[[dep]]; Cm <- as.matrix(dat[, ctrl, with = FALSE]); H <- as.matrix(dat[, .(H_wet_z, H_dry_z)])
  Zy <- demean2(cbind(y, Cm), g1, g2); yt <- Zy[, 1]; Ct <- Zy[, -1, drop = FALSE]
  stat <- function(Hh) { Hd <- demean2(Hh, g1, g2); X <- cbind(Hd, Ct); b <- solve(crossprod(X), crossprod(X, yt)); b[1:2] }
  obs <- stat(H); cells <- split(seq_len(nrow(dat)), g2)
  dr <- replicate(B, { Hp <- H; for (ix in cells) if (length(ix) > 1) Hp[ix, ] <- H[sample(ix), ]; stat(Hp) })
  data.frame(term = c("H_wet_z", "H_dry_z"), est = obs, p_perm = c(mean(abs(dr[1, ]) >= abs(obs[1])), mean(abs(dr[2, ]) >= abs(obs[2])))) }
FAST <- (exists(".R28_FAST") && isTRUE(.R28_FAST)) || nzchar(Sys.getenv("R28_FAST"))   # FAST: bo qua hoan vi cap ho (21 phut), chi anh huong bang RC3, khong anh huong hinh
if (!FAST) { t0 <- Sys.time(); RI_f <- perm_hh("y_farm", sf); RI_s <- perm_hh("y_shrimp", sh); lg("C3 hoan vi ho (", fmt(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " phut):")
RI_h <- rbind(cbind(sample = "Tong nuoi trong", RI_f), cbind(sample = "Tom", RI_s)); wt(RI_h, "RC3_permutation_household"); print(RI_h) } else { lg("C3 bo qua (che do FAST); giu nguyen RC3 cu"); RI_h <- NULL }

# C4. do manh va kiem dinh tuong duong: MDE (80%, 5%) = 2,8 SE; TOST voi |delta| = 0,20 (phan doi ung 20% gia tri)
rob <- function(lbl, dep, dat, trans = identity) { dd <- copy(dat); dd[, yy := trans(get(dep))]; r <- hh_fit("yy", c("H_wet_z", "H_dry_z", ctrl), dd)[1:2, ]
  data.table(outcome = lbl, term = r$term, est = r$est, se = r$se, p = r$p, N = r$N) }
w99 <- function(v) { q <- quantile(v, c(0.01, 0.99), na.rm = TRUE); pmin(pmax(v, q[1]), q[2]) }
C4 <- rbindlist(list(rob("Tổng nuôi trồng (cơ sở)", "y_farm", sf), rob("Tổng nuôi trồng, winsor 1–99", "y_farm", sf, w99),
                     rob("Tổng nuôi trồng, bỏ 1% cao nhất", "y_farm", sf[y_farm <= quantile(y_farm, 0.99, na.rm = TRUE)]),
                     rob("Tôm (cơ sở)", "y_shrimp", sh), rob("Tôm, winsor 1–99", "y_shrimp", sh, w99),
                     rob("Tôm, log(1+x) thay asinh", "aq_v_shrimp", sh, function(v) log1p(pmax(v, 0) / 1000))))
DEL <- 0.20; C4[, `:=`(lo = est - 1.96 * se, hi = est + 1.96 * se, mde = 2.8 * se, tost_p = pmax(pt((est - DEL) / se, 1e6), pt(-(est + DEL) / se, 1e6)) )]
C4[, tost_p := pmax(pnorm((est - DEL) / se), pnorm((-DEL - est) / se))]   # p cua TOST (hai mot phia): can ca hai < 0,05 => KTC 90% nam trong (-delta, delta)
C4[, equiv := (est - 1.645 * se > -DEL) & (est + 1.645 * se < DEL)]
wt(C4, "RC4_outliers_equivalence"); lg("C4 outlier/tuong duong:"); print(C4[, .(outcome, term, est = round(est, 3), se = round(se, 3), mde = round(mde, 2), equiv)])
pC4 <- ggplot(C4, aes(est, factor(outcome, rev(unique(outcome))), colour = term)) + annotate("rect", xmin = -DEL, xmax = DEL, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_vline(xintercept = 0, linewidth = 0.3) + geom_errorbarh(aes(xmin = est - 1.645 * se, xmax = est + 1.645 * se), height = 0, linewidth = 0.5, position = position_dodge(0.5)) +
  geom_point(size = 1.6, position = position_dodge(0.5)) +
  scale_colour_manual(values = c(H_wet_z = BLUE, H_dry_z = GOLD), labels = c(H_wet_z = "Ướt", H_dry_z = "Khô–nóng"), name = NULL) +
  labs(x = "Hệ số trên 1 độ lệch chuẩn hiểm họa; KTC 90%; nền xám: miền tương đương ±0,20", y = NULL, title = "Độ vững với giá trị ngoại lai và kiểm định tương đương") + th()
savef(pC4, "RC_equivalence", 6.6, 3.0)

# =============================================================================
# D. CHI SO LUOI vs THIEN TAI XA BAO CAO (AUC)
# =============================================================================
cm <- s[!is.na(cm_drought) & year >= 2004, .(cm_drought = max(cm_drought), cm_flood = max(cm_flood), cm_storm = max(cm_storm), hot = mean(cl_hot35_12), spei = mean(cl_spei3_min),
  prdz = mean(pr_dry_z), dry_pc = mean(H_dry_pca), heavy = mean(cl_heavy50_12), rx1 = mean(cl_rx1_12), pranom = mean(cl_pr12_anom), wet_pc = mean(H_wet_pca), tc = mean(tc_di_12),
  mekong = max(mekong), prov_year = prov_year[1], province = province[1], cl_geo = cl_geo[1]), by = .(year, commune_id)]
auc <- function(score, y) { ok <- !is.na(score) & !is.na(y); sc <- score[ok]; yy <- y[ok]; n1 <- sum(yy == 1); n0 <- sum(yy == 0); if (n1 < 5 || n0 < 5) return(NA_real_)
  (sum(rank(sc)[yy == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
PAIRS <- data.table(event = c(rep("cm_drought", 4), rep("cm_flood", 4), "cm_storm"), index = c("hot", "spei", "prdz", "dry_pc", "heavy", "rx1", "pranom", "wet_pc", "tc"),
                    sign = c(1, -1, -1, 1, 1, 1, 1, 1, 1), label = c("Hạn: ngày Tmax≥35°C", "Hạn: SPEI-3 nhỏ nhất", "Hạn: mưa mùa khô", "Hạn: khô–nóng (PC1)", "Lũ: ngày mưa ≥50 mm", "Lũ: mưa ngày lớn nhất", "Lũ: dị thường mưa 12 tháng", "Lũ: chỉ số ướt (PC1)", "Bão: chỉ số phá hủy"))
cm[, py := as.integer(factor(prov_year))]; ids_py <- split(seq_len(nrow(cm)), cm$py)
D1 <- rbindlist(lapply(seq_len(nrow(PAIRS)), function(k) { sc <- PAIRS$sign[k] * cm[[PAIRS$index[k]]]; y <- cm[[PAIRS$event[k]]]
  within <- sc - ave(sc, cm$prov_year, FUN = function(z) mean(z, na.rm = TRUE)); within[cm$cl_geo != "district"] <- NA
  a0 <- auc(sc, y); aw <- auc(within, y)
  # (i) bootstrap theo cum tinh-nam (thay vi xa): bao ve phu thuoc khong gian trong cung tinh-nam
  bt <- replicate(500, { b <- sample(names(ids_py), replace = TRUE); i <- unlist(ids_py[b], use.names = FALSE); c(auc(sc[i], y[i]), auc(within[i], y[i])) })
  # (ii) hoan vi nhan su kien trong tung nam (H0: chi so khong chua thong tin) -- B = 1000, AUC gop
  yr_ix <- split(seq_len(nrow(cm)), cm$year)
  pm_ <- replicate(1000, { yp <- y; for (ix in yr_ix) yp[ix] <- y[sample(ix)]; auc(sc, yp) })
  # (iii) nguong p70/p80/p90
  th_ <- sapply(c(0.7, 0.8, 0.9), function(q) { pay <- sc >= quantile(sc, q, na.rm = TRUE); ok <- !is.na(sc) & !is.na(y); c(hit = mean(pay[ok & y == 1]), fp = mean(y[ok & pay] == 0), lift = mean(pay[ok & y == 1]) / (1 - q)) })
  # (iv) rieng ĐBSCL
  a_mk <- auc(sc[cm$mekong == 1], y[cm$mekong == 1])
  data.table(label = PAIRS$label[k], event = PAIRS$event[k], n = sum(!is.na(sc) & !is.na(y)), events = sum(y[!is.na(sc)], na.rm = TRUE), auc = a0, auc_lo_py = quantile(bt[1, ], .025, na.rm = TRUE), auc_hi_py = quantile(bt[1, ], .975, na.rm = TRUE),
             auc_within = aw, aucw_lo_py = quantile(bt[2, ], .025, na.rm = TRUE), aucw_hi_py = quantile(bt[2, ], .975, na.rm = TRUE), p_perm = mean(pm_ >= a0, na.rm = TRUE),
             auc_mekong = a_mk, hit70 = th_["hit", 1], fp70 = th_["fp", 1], lift70 = th_["lift", 1], hit80 = th_["hit", 2], fp80 = th_["fp", 2], lift80 = th_["lift", 2], hit90 = th_["hit", 3], fp90 = th_["fp", 3], lift90 = th_["lift", 3]) }))
wt(D1, "RD1_auc_robust"); lg("D1 AUC:"); print(D1[, .(label, auc = round(auc, 2), lo = round(auc_lo_py, 2), hi = round(auc_hi_py, 2), within = round(auc_within, 2), wlo = round(aucw_lo_py, 2), whi = round(aucw_hi_py, 2), p_perm, mekong = round(auc_mekong, 2), lift70 = round(lift70, 2), lift80 = round(lift80, 2), lift90 = round(lift90, 2))])
pD <- ggplot(D1, aes(auc_within, factor(label, rev(label)))) + geom_vline(xintercept = 0.5, colour = RED, linetype = 2) +
  geom_errorbarh(aes(xmin = aucw_lo_py, xmax = aucw_hi_py), height = 0, colour = BLUE, linewidth = 0.5, na.rm = TRUE) + geom_point(colour = BLUE, size = 1.8, na.rm = TRUE) +
  geom_errorbarh(aes(y = factor(label, rev(label)), xmin = auc_lo_py, xmax = auc_hi_py), height = 0, colour = GOLD, linewidth = 0.5, position = position_nudge(y = 0.28)) +
  geom_point(aes(x = auc), colour = GOLD, size = 1.6, position = position_nudge(y = 0.28)) +
  labs(x = "AUC (đỏ: 0,5 = ngẫu nhiên). Xanh: trong cùng tỉnh-năm; vàng: gộp. KTC 95% bootstrap theo cụm tỉnh-năm", y = NULL, title = "Độ vững của AUC: bootstrap theo cụm tỉnh-năm") + th()
savef(pD, "RD_auc_cluster", 6.6, 3.4)
lg("XONG"); saveRDS(list(A0 = A0[, setdiff(names(A0), "vc")], alt = alt[, !"vc"], loo = loo, ri_tab = ri_tab, B1 = B1, pb = pb, ri = unique(ri[, .(outcome, obs, p_ri, n_assign)]), B4 = B4, C0 = C0, C0s = C0s, res = res, cl_tab = cl_tab, RI_h = RI_h, C4 = C4, D1 = D1), file.path(OUT, "R28_results.rds"))
