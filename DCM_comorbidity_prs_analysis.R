# DCM co-morbidities letter analysis
#   A  OR per SD, European cases only: VCCRI vs MGRB half A, AGHA vs MGRB half B.
#   B  OR per SD, all cases vs 1000G NFE, scored on the 1000G-shared variant set.
#   C  AUC, all cases vs all MGRB controls, DCM-PRS alone and with LVESVi and the replicated co-morbidity PRS
#      that also held against 1000G (T2D excluded).
# Row order in A, B and C: DCM-PRS first (reference), then LVESVi PRS and the nine co-morbidity PRS by discovery (VCCRI) OR.
# LVESVi PRS (PGS000666): 21 of its 28 variants in A and C (allele-order-agnostic matching, no MAF filter; see README); 17 of 28 in B
#   (separate 1000G-shared run).
# Variable names: B holds the panel C (AUC) results and C holds the panel B (1000G) results.
# Usage (from a directory holding data/ and figures/, see README): Rscript DCM_comorbidity_prs_analysis.R
# Panel C: AUC 95% CI by DeLong; P (asterisks) = likelihood ratio test of each model vs DCM-PRS alone (nested models;
#   Vickers et al. BMC Med Res Methodol 2011). AUCs are apparent (fitted and evaluated in the same samples).
# The numbers table also holds panel "A_adjDCM": each PRS in A refitted with the DCM-PRS as a covariate (text only, not plotted).
# Writes figures/Figure1.pdf, figures/Figure1.png and figures/Figure1_numbers.tsv. Prints aggregates only.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(pROC) })

D <- "data"; FIG <- "figures"
AGHA_ID <- "<AGHA_ID_PATTERN>"                     # regular expression matching AGHA case IDs (not public); other cases are VCCRI
nf <- function(x) trimws(format(x, big.mark = ","))
stopifnot(dir.exists(D), dir.exists(FIG))
rd <- function(p) { d <- read.delim(p, check.names = FALSE, stringsAsFactors = FALSE); names(d)[1] <- "IID"; d$IID <- trimws(d$IID); d }

# Table rows: DCM-PRS, LVESVi PRS, then the nine co-morbidity PRS alphabetically by label (display order is ORD below).
# key = stem in the MGRB run; g1k = name in the 1000G run.
TR <- data.frame(
  key   = c("dcm", "lvesvi", "alcohol", "af", "bmi", "crf", "cad", "hypertension", "ageing", "t1d", "t2d"),
  mgrb  = c("pgs004946", "pgs000666", "alcohol", "af", "bmi", "crf", "cad", "hypertension", "ageing", "t1d", "t2d"),
  g1k   = c("DCM", "LVESVi", "alcohol", "atrial_fibrillation", "BMI", "cardiorespiratory_fitness", "coronary_artery_disease",
            "hypertension", "ageing", "T1D", "T2D"),
  label = c("DCM-PRS", "LVESVi PRS", "Alcohol use disorder", "Atrial fibrillation", "Body mass index", "Cardiorespiratory fitness",
            "Coronary heart disease", "Hypertension", "Longevity", "Type 1 diabetes", "Type 2 diabetes"),
  short = c("DCM-PRS", "LVESVi", "AUD", "AF", "BMI", "CRF", "CHD", "Hypertension", "Longevity", "T1D", "T2D"),
  stringsAsFactors = FALSE)
REF <- c("dcm", "lvesvi")                          # reference PRS, outside the BH family
FAM <- setdiff(TR$key, REF)                        # the nine co-morbidity PRS (BH family)
stopifnot(length(FAM) == 9, identical(TR$label[TR$key %in% FAM], sort(TR$label[TR$key %in% FAM])))

# ---- samples ---------------------------------------------------------------------------------------------
split <- rd(file.path(D, "samples", "mgrb_split_halves.tsv"))
stopifnot(all(split$half %in% c("A", "B")), !anyDuplicated(split$IID))
noneur <- trimws(readLines(file.path(D, "samples", "non_european_cases.txt"), warn = FALSE)); noneur <- noneur[noneur != ""]

# ---- MGRB run: cases and controls on one score scale per trait ---------------------------------------------
mg_case <- lapply(setNames(TR$mgrb, TR$key), function(s) rd(file.path(D, "mgrb_run", sprintf("prs_dcm_%s.sscore", s))))
mg_ctrl <- lapply(setNames(TR$mgrb, TR$key), function(s) rd(file.path(D, "mgrb_run", sprintf("prs_mgrb_%s.sscore", s))))

case_ids <- setdiff(Reduce(intersect, lapply(mg_case, `[[`, "IID")), noneur)
stopifnot(all(noneur %in% mg_case$af$IID))
arm <- ifelse(grepl(AGHA_ID, case_ids), "AGHA", "VCCRI")
stopifnot(all(c("VCCRI", "AGHA") %in% arm))
N <- c(VCCRI = sum(arm == "VCCRI"), AGHA = sum(arm == "AGHA"), halfA = sum(split$half == "A"), halfB = sum(split$half == "B"),
       case = length(case_ids), mgrb = nrow(split))

M <- rbind(data.frame(IID = case_ids, y = 1L, grp = arm), data.frame(IID = split$IID, y = 0L, grp = paste0("half", split$half)))
for (k in TR$key) M[[k]] <- c(mg_case[[k]]$SCORE1_AVG[match(case_ids, mg_case[[k]]$IID)],
                               mg_ctrl[[k]]$SCORE1_AVG[match(split$IID, mg_ctrl[[k]]$IID)])
stopifnot(!anyNA(M[, TR$key]), !any(case_ids %in% split$IID))

zstd <- function(d, ref, cols) { for (s in cols) d[[paste0("z_", s)]] <- (d[[s]] - mean(ref[[s]])) / sd(ref[[s]]); d }
fit_or <- function(d, s) {
  f <- glm(as.formula(paste0("y ~ z_", s)), binomial, d)
  ci <- suppressMessages(confint(f, paste0("z_", s)))         # profile likelihood
  data.frame(OR = exp(coef(f)[[2]]), lo = exp(ci[[1]]), hi = exp(ci[[2]]), p = summary(f)$coef[2, 4], beta = coef(f)[[2]])
}

# ---- Panel A ---------------------------------------------------------------------------------------------
armfit <- function(case_grp, half) {
  d <- M[M$grp %in% c(case_grp, half), ]; d <- zstd(d, d[d$y == 0, ], TR$key)
  out <- do.call(rbind, lapply(TR$key, function(s) cbind(key = s, fit_or(d, s))))
  out$n_case <- sum(d$y == 1); out$n_ctrl <- sum(d$y == 0); out
}
A_v <- armfit("VCCRI", "halfA"); A_g <- armfit("AGHA", "halfB")
A_v$q <- NA; A_v$q[A_v$key %in% FAM] <- p.adjust(A_v$p[A_v$key %in% FAM], "BH")
S <- A_v$key[A_v$key %in% FAM & A_v$q < 0.05]                                   # carried forward
A_g$q <- NA; A_g$q[A_g$key %in% S] <- p.adjust(A_g$p[A_g$key %in% S], "BH")
A_g$replicated <- A_g$key %in% S & !is.na(A_g$q) & A_g$q < 0.05 &
  sign(A_g$beta) == sign(A_v$beta[match(A_g$key, A_v$key)])
REP <- A_g$key[A_g$replicated]
cat("Discovery significant:", paste(S, collapse = ", "), "\nReplicated:", paste(REP, collapse = ", "), "\n")
stopifnot(setequal(REP, c("cad", "hypertension", "bmi", "t2d", "ageing")))   # the five named in the letter

# Same arms, each PRS adjusted for the DCM-PRS (BH families as above: nine in discovery, the carried-forward set in validation).
fit_or_adj <- function(d, s) {
  f <- glm(as.formula(paste0("y ~ z_", s, " + z_dcm")), binomial, d)
  ci <- suppressMessages(confint(f, paste0("z_", s)))
  data.frame(OR = exp(coef(f)[[2]]), lo = exp(ci[[1]]), hi = exp(ci[[2]]), p = summary(f)$coef[2, 4])
}
armfit_adj <- function(case_grp, half) {
  d <- M[M$grp %in% c(case_grp, half), ]; d <- zstd(d, d[d$y == 0, ], TR$key)
  do.call(rbind, lapply(setdiff(TR$key, "dcm"), function(s) cbind(key = s, fit_or_adj(d, s))))
}
J_v <- armfit_adj("VCCRI", "halfA"); J_g <- armfit_adj("AGHA", "halfB")
J_v$q <- NA; J_v$q[J_v$key %in% FAM] <- p.adjust(J_v$p[J_v$key %in% FAM], "BH")
J_g$q <- NA; J_g$q[J_g$key %in% S] <- p.adjust(J_g$p[J_g$key %in% S], "BH")
cat("Replicated traits with q<0.05 in validation after DCM-PRS adjustment:",
    paste(intersect(REP, J_g$key[!is.na(J_g$q) & J_g$q < 0.05]), collapse = ", "), "\n")

# ---- Panel C (AUC) -----------------------------------------------------------------------------------------
P <- zstd(M, M[M$y == 0, ], TR$key)
# Row order for every panel: DCM-PRS first, then all other PRS by discovery (VCCRI) OR, highest at the top.
REST <- setdiff(TR$key, "dcm")
ORD <- c("dcm", REST[order(-A_v$OR[match(REST, A_v$key)])])
stopifnot(setequal(ORD, TR$key))
REP <- ORD[ORD %in% REP]
fitm <- function(v) glm(as.formula(paste("y ~", paste0("z_", v, collapse = " + "))), binomial, P)
mk <- function(v) roc(P$y, fitted(fitm(v)), quiet = TRUE, direction = "<")
lrt <- function(v) anova(fitm("dcm"), fitm(v), test = "Chisq")[2, "Pr(>Chi)"]   # nested: DCM-PRS alone vs v
# T2D replicated in panel A but was not significant against 1000G (panel B), so it is left out of the AUC models.
AUC_SET <- c(setdiff(REP, "t2d"), "lvesvi")
ADD <- ORD[ORD %in% AUC_SET]                                                 # single additions, in panel A order
models <- c(list("DCM-PRS only" = "dcm"),
            setNames(lapply(ADD, function(s) c("dcm", s)), paste("+", TR$short[match(ADD, TR$key)])),
            list("+ All" = c("dcm", AUC_SET)))
B <- do.call(rbind, lapply(names(models), function(n) {
  r <- mk(models[[n]]); ci <- ci.auc(r, method = "delong")
  p <- if (n == "DCM-PRS only") NA else lrt(models[[n]])
  data.frame(model = n, AUC = as.numeric(auc(r)), lo = ci[1], hi = ci[3], p = p) }))

# ---- Panel B (1000G) -----------------------------------------------------------------------------------------
g_case <- lapply(setNames(TR$g1k, TR$key), function(s) rd(file.path(D, "1000g_run", sprintf("prs_dcm_%s.sscore", s))))
g_ref  <- lapply(setNames(TR$g1k, TR$key), function(s) rd(file.path(D, "1000g_run", sprintf("prs_1000G_%s.sscore", s))))
ref_ids <- g_ref$dcm$IID
stopifnot(all(sapply(g_ref, function(x) setequal(x$IID, ref_ids))),
          all(sapply(g_case, function(x) all(case_ids %in% x$IID))))
G <- data.frame(IID = c(case_ids, ref_ids), y = c(rep(1L, length(case_ids)), rep(0L, length(ref_ids))))
N[["g1k"]] <- length(ref_ids)
for (k in TR$key) G[[k]] <- c(g_case[[k]]$SCORE1_AVG[match(case_ids, g_case[[k]]$IID)],
                               g_ref[[k]]$SCORE1_AVG[match(ref_ids, g_ref[[k]]$IID)])
stopifnot(!anyNA(G[, TR$key]))
G <- zstd(G, G[G$y == 0, ], TR$key)
C <- do.call(rbind, lapply(TR$key, function(s) cbind(key = s, fit_or(G, s))))
C$q <- NA; C$q[C$key %in% FAM] <- p.adjust(C$p[C$key %in% FAM], "BH")

# ---- numbers table -----------------------------------------------------------------------------------------
lab <- function(k) TR$label[match(k, TR$key)]
num <- rbind(
  data.frame(panel = "A", row = lab(A_v$key), group = "VCCRI vs MGRB half A", n_case = N[["VCCRI"]], n_ctrl = N[["halfA"]],
             estimate = A_v$OR, lo = A_v$lo, hi = A_v$hi, p = A_v$p, q = A_v$q, replicated = NA),
  data.frame(panel = "A", row = lab(A_g$key), group = "AGHA vs MGRB half B", n_case = N[["AGHA"]], n_ctrl = N[["halfB"]],
             estimate = A_g$OR, lo = A_g$lo, hi = A_g$hi, p = A_g$p, q = A_g$q,
             replicated = ifelse(A_g$key %in% S, A_g$replicated, NA)),
  data.frame(panel = "C", row = B$model, group = "All cases vs all MGRB (apparent AUC; P = LRT vs DCM-PRS alone)", n_case = N[["case"]], n_ctrl = N[["mgrb"]],
             estimate = B$AUC, lo = B$lo, hi = B$hi, p = B$p, q = NA, replicated = NA),
  data.frame(panel = "B", row = lab(C$key), group = "All cases vs 1000G NFE", n_case = N[["case"]], n_ctrl = N[["g1k"]],
             estimate = C$OR, lo = C$lo, hi = C$hi, p = C$p, q = C$q, replicated = NA),
  data.frame(panel = "A_adjDCM", row = lab(J_v$key), group = "VCCRI vs MGRB half A, adjusted for DCM-PRS", n_case = N[["VCCRI"]],
             n_ctrl = N[["halfA"]], estimate = J_v$OR, lo = J_v$lo, hi = J_v$hi, p = J_v$p, q = J_v$q, replicated = NA),
  data.frame(panel = "A_adjDCM", row = lab(J_g$key), group = "AGHA vs MGRB half B, adjusted for DCM-PRS", n_case = N[["AGHA"]],
             n_ctrl = N[["halfB"]], estimate = J_g$OR, lo = J_g$lo, hi = J_g$hi, p = J_g$p, q = J_g$q, replicated = NA))
write.table(num, file.path(FIG, "Figure1_numbers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pr <- num; pr[, c("estimate", "lo", "hi")] <- round(pr[, c("estimate", "lo", "hi")], 3); pr$p <- signif(pr$p, 2); pr$q <- signif(pr$q, 2)
print(pr, row.names = FALSE)

# ---- figure ------------------------------------------------------------------------------------------------
# A and C: plot plus a text column on the same y scale; the top of each text column is its header.
# Asterisks in A and B: co-morbidity PRS use the BH-adjusted q (*q<0.05, **q<0.01, ***q<0.001); validation traits not
# carried forward from discovery have no q and no asterisk. DCM-PRS and LVESVi PRS are outside the BH family and use P.
# Asterisks in C use the likelihood ratio test P vs DCM-PRS alone.
FONT <- "Helvetica"; INK <- "grey15"; RULE <- "grey60"
COL <- c(VCCRI = "#0072B2", AGHA = "#D55E00", C = "grey15")
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
fmt <- function(e, l, h, p) sprintf("%.2f (%.2f–%.2f)%s", e, l, h, stars(p))
NR <- nrow(TR)
ypos <- function(i) NR - i + 1 + ifelse(i > 1, 0, 0.35)          # small gap below DCM-PRS, no box
HDR <- NR + 1.75                                                 # header row height
YR <- c(0.4, NR + 2.85)
base_theme <- theme_classic(base_size = 8, base_family = FONT) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.text = element_text(colour = INK, size = 7.5), axis.title.x = element_text(colour = INK, size = 8),
        axis.text.y = element_text(hjust = 0, colour = INK, size = 7.5),
        axis.line.x = element_line(colour = INK, linewidth = 0.35), axis.ticks.x = element_line(colour = INK, linewidth = 0.35),
        plot.margin = margin(2, 2, 2, 2), legend.position = "none",
        plot.tag = element_text(face = "bold", size = 11, family = FONT), plot.tag.position = c(0, 1))
OR_X <- scale_x_continuous(limits = c(0.5, 2.2), breaks = seq(0.5, 2, 0.5), expand = c(0, 0))
yrow <- scale_y_continuous(limits = YR, breaks = ypos(seq_len(NR)), labels = TR$label[match(ORD, TR$key)], expand = c(0, 0))
# h: header rows with columns x, title (bold), sub (plain, below the title)
textcol <- function(t, h, xmax) ggplot(t, aes(x, y, label = lab, colour = arm)) +
  geom_text(hjust = 0, size = 2.6, family = FONT) +
  geom_text(data = h, aes(y = HDR + 0.8, label = title), hjust = 0, size = 2.6, family = FONT, fontface = "bold") +
  geom_text(data = h, aes(y = HDR - 0.1, label = sub), hjust = 0, size = 2.5, family = FONT, lineheight = 0.95) +
  scale_colour_manual(values = COL) + scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
  scale_y_continuous(limits = YR, expand = c(0, 0)) + theme_void(base_family = FONT) +
  theme(legend.position = "none", plot.margin = margin(2, 2, 2, 6))

# Panel A
dA <- rbind(data.frame(key = A_v$key, arm = "VCCRI", OR = A_v$OR, lo = A_v$lo, hi = A_v$hi, p = ifelse(A_v$key %in% REF, A_v$p, A_v$q), dy = 0.17),
            data.frame(key = A_g$key, arm = "AGHA",  OR = A_g$OR, lo = A_g$lo, hi = A_g$hi, p = ifelse(A_g$key %in% REF, A_g$p, A_g$q), dy = -0.17))
dA$y <- ypos(match(dA$key, ORD)) + dA$dy
pA <- ggplot(dA, aes(OR, y, colour = arm)) +
  geom_vline(xintercept = 1, colour = RULE, linewidth = 0.3, linetype = "22") +
  geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 0.45) + geom_point(size = 1.6) +
  scale_colour_manual(values = COL) + OR_X + yrow +
  labs(x = "OR per 1 SD of PRS (95% CI)", tag = "A") + base_theme
tA <- data.frame(x = ifelse(dA$arm == "VCCRI", 0, 1), y = ypos(match(dA$key, ORD)), arm = dA$arm,
                 lab = fmt(dA$OR, dA$lo, dA$hi, dA$p))
hA <- data.frame(x = c(0, 1), arm = c("VCCRI", "AGHA"), title = c("Discovery", "Validation"),
                 sub = c(sprintf("VCCRI n=%s\nMGRB n=%s", nf(N[["VCCRI"]]), nf(N[["halfA"]])),
                         sprintf("AGHA n=%s\nMGRB n=%s", nf(N[["AGHA"]]), nf(N[["halfB"]]))))
pAt <- textcol(tA, hA, 2)

# Panel C (AUC; data in B): horizontal bars from 0.50, whiskers, value + asterisks at the right.
B$y <- rev(seq_len(nrow(B)))
B$fill <- ifelse(B$model == "DCM-PRS only", "base", ifelse(B$model == "+ All", "all", ifelse(B$model == "+ LVESVi", "lv", "single")))
pB <- ggplot(B, aes(y = y)) +
  geom_rect(aes(xmin = 0.5, xmax = AUC, ymin = y - 0.36, ymax = y + 0.36, fill = fill)) +
  geom_vline(xintercept = B$AUC[1], colour = "grey45", linewidth = 0.3, linetype = "22") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.3, linewidth = 0.35, colour = INK, orientation = "y") +
  geom_text(aes(x = hi + 0.003, label = sprintf("%.3f%s", AUC, stars(p))), hjust = 0, size = 2.6, family = FONT, colour = INK) +
  scale_fill_manual(values = c(base = "grey45", single = "#9BD3BF", lv = "#CC79A7", all = "#009E73")) +
  scale_x_continuous(limits = c(0.5, 0.8), breaks = seq(0.5, 0.75, 0.05), expand = c(0, 0)) +
  scale_y_continuous(breaks = B$y, labels = B$model, limits = c(0.4, nrow(B) + 0.6), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "AUC (95% CI)", tag = "C", subtitle = sprintf("All cases (n=%s) vs all MGRB controls (n=%s)", nf(N[["case"]]), nf(N[["mgrb"]]))) + base_theme +
  theme(axis.text.y = element_text(hjust = 1, colour = INK, size = 7.5),
        axis.line.y = element_line(colour = INK, linewidth = 0.35),
        plot.subtitle = element_text(size = 7.5, colour = INK, family = FONT, hjust = 0.5, margin = margin(b = 6)))

# Panel B (1000G; data in C): same rows and axis as A, one estimate per row.
dC <- data.frame(key = C$key, OR = C$OR, lo = C$lo, hi = C$hi, p = ifelse(C$key %in% REF, C$p, C$q), y = ypos(match(C$key, ORD)))
pC <- ggplot(dC, aes(OR, y)) +
  geom_vline(xintercept = 1, colour = RULE, linewidth = 0.3, linetype = "22") +
  geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 0.45, colour = COL[["C"]]) +
  geom_point(size = 1.6, colour = COL[["C"]]) + OR_X + yrow +
  labs(x = "OR per 1 SD of PRS (95% CI)", tag = "B") + base_theme
tC <- data.frame(x = 0, y = dC$y, arm = "C", lab = fmt(dC$OR, dC$lo, dC$hi, dC$p))
hC <- data.frame(x = 0, arm = "C", title = "1000G reference", sub = sprintf("All cases n=%s\n1000G NFE n=%s", nf(N[["case"]]), nf(N[["g1k"]])))
pCt <- textcol(tC, hC, 1)

top <- (pA | pAt) + plot_layout(widths = c(1.2, 1))
bottom <- (pC | pCt | plot_spacer() | pB) + plot_layout(widths = c(0.75, 0.6, 0.06, 1.15))
fig <- wrap_elements(top) / wrap_elements(bottom) + plot_layout(heights = c(1, 1))
ggsave(file.path(FIG, "Figure1.pdf"), fig, width = 180, height = 175, units = "mm", device = cairo_pdf)
ggsave(file.path(FIG, "Figure1.png"), fig, width = 180, height = 175, units = "mm", dpi = 600, bg = "white")
cat("written:", file.path(FIG, c("Figure1.pdf", "Figure1.png", "Figure1_numbers.tsv")), "\n")
