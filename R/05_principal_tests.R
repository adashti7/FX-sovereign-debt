# 05_principal_tests.R - principal repayments only: event-time coefficients, pooled windows, placebo on P_evt.
#   event-time: r on P at k = -5..+5 (k > 0 means the repayment settles k windows later) + controls
#   placebo:    reassign the real repayment amounts to random non-event windows, 300 draws, empirical p for P_evt
source("R/00_config.R"); options(mx.no_main = TRUE); source("R/03_build_panel.R"); source("R/04_run_regressions.R")
CTRL <- c("dlnVIX", "dUST2", "rBrent")
w <- window_frame("fix")

# ---- event-time block (11 coefficients, for the chart) ----
et <- w
for (k in -5:5) et[[paste0("P_k", ifelse(k < 0, "m", "p"), abs(k))]] <- if (k >= 0) lead(w$P0, k) else lag(w$P0, -k)
kcols <- grep("^P_k", names(et), value = TRUE)
f_et <- fit_hac(et, c(kcols, CTRL))
et_tab <- tibble(term = rownames(f_et$ct), coef = f_et$ct[, 1], se = f_et$ct[, 2], p = f_et$ct[, 4]) %>%
  filter(grepl("^P_k", term)) %>%
  mutate(k = as.integer(sub("P_k[mp]", "", term)) * ifelse(grepl("P_km", term), -1, 1)) %>%
  arrange(k) %>% mutate(ci_lo = coef - 1.96 * se, ci_hi = coef + 1.96 * se)
write_csv(et_tab, file.path(OUT, "mx_principal_event_time.csv"))

png(file.path(OUT, "mx_principal_event_time.png"), width = 1600, height = 900, res = 180)
par(mar = c(4.5, 5, 3, 1))
plot(et_tab$k, et_tab$coef, pch = 19, ylim = range(c(et_tab$ci_lo, et_tab$ci_hi)), xaxt = "n",
     xlab = "Fixing windows until repayment settles (positive = before, negative = after)",
     ylab = "% change in USD/MXN per MXN 1bn foreign-held principal",
     main = "Peso response around Mbono principal repayments (HAC 95% intervals)")
axis(1, at = -5:5); abline(h = 0, col = "grey50"); abline(v = 0.5, lty = 2, col = "grey60")
arrows(et_tab$k, et_tab$ci_lo, et_tab$k, et_tab$ci_hi, angle = 90, code = 3, length = 0.04)
dev.off()

# ---- pooled windows and placebo on P_evt ----
pw <- pooled_frame(w)
f_pool <- fit_hac(pw, c("P_pre", "P_evt", "P_post", CTRL))
actual <- f_pool$ct["P_evt", 1]

set.seed(7)
ev_idx <- which(w$P0 > 0); ev_amt <- w$P0[ev_idx]
pool   <- which(w$P0 == 0 & seq_len(nrow(w)) > 10 & seq_len(nrow(w)) < nrow(w) - 10)
draws <- replicate(300, {
  fake <- w; fake$P0 <- 0; fake$P0[sample(pool, length(ev_idx))] <- ev_amt
  fp <- pooled_frame(fake); fit_hac(fp, c("P_pre", "P_evt", "P_post", CTRL))$ct["P_evt", 1]
})
summ <- tibble(P_evt = actual, hac_se = f_pool$ct["P_evt", 2], hac_p = f_pool$ct["P_evt", 4],
               placebo_sd = sd(draws), placebo_p_two_sided = mean(abs(draws) >= abs(actual)),
               placebo_p_one_sided_neg = mean(draws <= actual), n_events = length(ev_idx), n_windows = nrow(w))
write_csv(summ, file.path(OUT, "mx_principal_placebo.csv"))
write_csv(tibble(draw = draws), file.path(OUT, "mx_principal_placebo_draws.csv"))

# ---- economic size: largest foreign-held maturities and their implied move ----
big <- w %>% filter(P0 > 0) %>% arrange(desc(P0)) %>% slice_head(n = 10) %>%
  transmute(date, foreign_principal_bn = round(P0, 1), implied_move_pct = round(P0 * actual, 3))
write_csv(big, file.path(OUT, "mx_principal_largest.csv"))

# ---- placebo histogram ----
png(file.path(OUT, "mx_principal_placebo_hist.png"), width = 1600, height = 900, res = 180)
par(mar = c(4.5, 5, 3, 1))
hist(draws, breaks = 40, col = "grey85", border = "white",
     main = "Payment-window coefficient: 300 placebo draws vs actual",
     xlab = "P_evt from randomly placed repayments (% per MXN 1bn)")
abline(v = actual, col = "firebrick", lwd = 3)
legend("topright", bty = "n", legend = sprintf("actual = %.4f, placebo p = %.2f", actual, summ$placebo_p_two_sided))
dev.off()

# ---- foreign share versus days to maturity (the early-exit chart) ----
hold <- read_csv(file.path(RAW, "mx_issue_holdings.csv"), show_col_types = FALSE) %>%
  mutate(date = as.Date(date), maturity = as.Date(maturity), share = foreign / outstanding,
         days_to_mat = as.numeric(maturity - date)) %>%
  filter(maturity <= max(date), days_to_mat >= 0, days_to_mat <= 3 * 365, is.finite(share))
prof <- hold %>% mutate(bucket = floor(days_to_mat / 30) * 30) %>%
  group_by(bucket) %>% summarise(median_share = median(share), q25 = quantile(share, .25),
                                  q75 = quantile(share, .75), n = n_distinct(series), .groups = "drop")
write_csv(prof, file.path(OUT, "mx_foreign_share_by_days_to_maturity.csv"))
png(file.path(OUT, "mx_foreign_share_rundown.png"), width = 1600, height = 900, res = 180)
par(mar = c(4.5, 5, 3, 1))
plot(prof$bucket, 100 * prof$median_share, type = "n", xlim = rev(range(prof$bucket)),
     ylim = c(0, 100 * max(prof$q75)), xlab = "Days until maturity",
     ylab = "Foreign share of the bond (%)",
     main = sprintf("Foreign ownership runs down before maturity (%d matured Mbonos, 2010-2026)", n_distinct(hold$series)))
polygon(c(prof$bucket, rev(prof$bucket)), 100 * c(prof$q25, rev(prof$q75)), col = "grey88", border = NA)
lines(prof$bucket, 100 * prof$median_share, lwd = 3)
legend("topright", bty = "n", legend = c("median", "interquartile range"), lwd = c(3, 10), col = c("black", "grey88"))
dev.off()

# ---- R-squared for the deck ----
r2 <- tibble(spec = c("net_flow_full_sample", "principal_pooled", "principal_event_time"),
             r2 = c(summary(fit_hac(w, c("N0", "N1", CTRL))$model)$r.squared,
                    summary(f_pool$model)$r.squared, summary(f_et$model)$r.squared))
write_csv(r2, file.path(OUT, "mx_r2.csv"))

print(et_tab %>% select(k, coef, se, p), n = 11); print(summ); print(big); print(r2)
