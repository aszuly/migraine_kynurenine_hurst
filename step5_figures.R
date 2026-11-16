# step5_figures.R
#
# Code used to draw the figures from the outputs of step3_main_analysis.m.
#
#   Figure 2   Plasma L-KYN against somatomotor-network mean Hurst exponent
#              (all 21 participants), with the leave-one-out estimates inset.
#   Figure 3   Partial Spearman correlations, 7 networks x 10 measures.
#   Figure S1  As Figure 2, participants without aura only (n = 16).
#
# Output: figures/*.pdf and figures/*.png (600 dpi), figures/captions.md

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(scales)
})

paper_theme <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(colour = "grey20", linewidth = 0.4),
      axis.ticks       = element_line(colour = "grey20", linewidth = 0.4),
      axis.text        = element_text(colour = "grey20"),
      legend.position  = "top",
      legend.title     = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", size = base_size + 1)
    )
}

# The two groups differ in colour and in marker shape.
aura_colours <- c(`Migraine without aura` = "#1f78b4", `Migraine with aura` = "#e31a1c")
aura_shapes  <- c(`Migraine without aura` = 16,        `Migraine with aura` = 17)

subjects <- read_csv("results/fig_subjects.csv", show_col_types = FALSE) |>
  mutate(group = factor(if_else(aura == 1, "Migraine with aura", "Migraine without aura"),
                        levels = names(aura_colours)))
tests <- read_csv("results/tests_70.csv", show_col_types = FALSE)
loo   <- read_csv("results/loo.csv", show_col_types = FALSE)
mwoa  <- read_csv("results/mwoa_only.csv", show_col_types = FALSE)

# plotmath label; "=" is written as text because an expression can hold only one "=="
fmt_p <- function(p) if (p < 0.001) "italic(p)*' < 0.001'" else sprintf("italic(p)*' = %.3f'", p)
rho_label <- function(rho, lo, hi, p) {
  sprintf("italic(rho)[partial]*' = %.2f (95%% CI %.2f to %.2f), '*%s", rho, lo, hi, fmt_p(p))
}

# ---------------------------------------------------------------------------
# Figure 2
# ---------------------------------------------------------------------------
main <- tests |> filter(Network == "SMN", Measure == "L-KYN")
loo_smn <- loo |> filter(Network == "SMN", Measure == "L-KYN")
stopifnot(nrow(main) == 1, nrow(loo_smn) == nrow(subjects))

scatter <- ggplot(subjects, aes(L_KYN, SMN_Hurst)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey25", fill = "grey75",
              alpha = 0.4, linewidth = 0.6) +
  geom_point(aes(colour = group, shape = group), size = 2.6, stroke = 0.4) +
  scale_colour_manual(values = aura_colours, name = NULL) +
  scale_shape_manual(values = aura_shapes, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(x = "Plasma L-kynurenine (nmol/L)",
       y = "Somatomotor network mean Hurst exponent") +
  annotate("text", x = -Inf, y = Inf, parse = TRUE, hjust = -0.04, vjust = 1.5, size = 3.3,
           label = rho_label(main$rho, main$CI_lo, main$CI_hi, main$p_perm)) +
  paper_theme()

loo_inset <- ggplot(loo_smn, aes(x = "", y = rho)) +
  geom_violin(fill = "grey90", colour = "grey60", linewidth = 0.3, width = 0.6) +
  geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1),
             size = 1.2, colour = "grey25", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
  coord_cartesian(ylim = c(0, 0.85)) +
  labs(x = NULL, y = expression("Leave-one-out " * italic(rho))) +
  paper_theme(base_size = 8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA))

fig2 <- scatter +
  inset_element(loo_inset, left = 0.68, bottom = 0.02, right = 0.995, top = 0.37,
                align_to = "panel")
ggsave("figures/Figure2_SMN_LKYN.pdf", fig2, width = 5.5, height = 4.4, device = cairo_pdf)
ggsave("figures/Figure2_SMN_LKYN.png", fig2, width = 5.5, height = 4.4, dpi = 600)

# ---------------------------------------------------------------------------
# Figure 3
# ---------------------------------------------------------------------------
net_order <- c("VN", "SMN", "DAN", "SVA", "LIM", "CON", "DMN")
met_order <- c("TRP", "L-KYN", "KYNA", "ANA", "3-HK", "XA", "3-HANA", "PICA", "QUINA", "QUINA/KYNA")
stopifnot(setequal(tests$Network, net_order), setequal(tests$Measure, met_order), nrow(tests) == 70)

hm <- tests |>
  mutate(Network = factor(Network, levels = rev(net_order)),
         Measure = factor(Measure, levels = met_order),
         sig     = q_family < 0.05,
         label   = paste0(sub("-", "−", sprintf("%.2f", rho)), if_else(sig, "*", "")))
rmax <- max(abs(hm$rho))

fig3 <- ggplot(hm, aes(Measure, Network, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = label, colour = abs(rho) > 0.5, fontface = if_else(sig, "bold", "plain")),
            size = 2.6, show.legend = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white")) +
  scale_fill_gradient2(low = "#2c7bb6", mid = "#f7f7f7", high = "#d7191c", midpoint = 0,
                       limits = c(-rmax, rmax), breaks = pretty_breaks(5),
                       name = expression(italic(rho)[partial])) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  paper_theme() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 0),
        legend.position = "right",
        legend.key.height = grid::unit(1.1, "cm"), legend.key.width = grid::unit(0.4, "cm"),
        legend.title = element_text(face = "bold", size = 10), legend.text = element_text(size = 9))
ggsave("figures/Figure3_network_by_measure.pdf", fig3, width = 6.6, height = 3.7, device = cairo_pdf)
ggsave("figures/Figure3_network_by_measure.png", fig3, width = 6.6, height = 3.7, dpi = 600)

# ---------------------------------------------------------------------------
# Figure S1
# ---------------------------------------------------------------------------
m <- mwoa |> filter(Network == "SMN")
figS1 <- ggplot(filter(subjects, aura == 0), aes(L_KYN, SMN_Hurst)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey25", fill = "grey75",
              alpha = 0.4, linewidth = 0.6) +
  geom_point(colour = aura_colours[["Migraine without aura"]], shape = 16, size = 2.6) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(x = "Plasma L-kynurenine (nmol/L)",
       y = "Somatomotor network mean Hurst exponent") +
  annotate("text", x = -Inf, y = Inf, parse = TRUE, hjust = -0.04, vjust = 1.5, size = 3.3,
           label = rho_label(m$rho, m$CI_lo, m$CI_hi, m$p_perm)) +
  paper_theme()
ggsave("figures/FigureS1_SMN_LKYN_without_aura.pdf", figS1, width = 5.0, height = 4.2, device = cairo_pdf)
ggsave("figures/FigureS1_SMN_LKYN_without_aura.png", figS1, width = 5.0, height = 4.2, dpi = 600)

# ---------------------------------------------------------------------------
# Captions
# ---------------------------------------------------------------------------
prim  <- tests |> filter(Measure == "L-KYN")
sig_n <- prim |> filter(q_family < 0.05) |> pull(Network)
cap <- c(
  "# Figure captions (numbers filled in by step5_figures.R)", "",
  sprintf(paste0(
    "**Figure 2.** Plasma L-kynurenine concentration and mean Hurst exponent of the somatomotor network ",
    "in the %d participants. Circles: migraine without aura (n = %d); triangles: migraine with aura (n = %d). ",
    "The least-squares line and its 95%% confidence band are shown as a visual guide; ",
    "the statistic is the partial Spearman correlation adjusted for age and aura status ",
    "(rho = %.2f, 95%% BCa CI %.2f to %.2f, permutation p = %.4f; q = %.3f across the seven networks, ",
    "q = %.2f across all 70 tests). Inset: the %d leave-one-out estimates (range %.2f to %.2f)."),
    nrow(subjects), sum(subjects$aura == 0), sum(subjects$aura == 1), main$rho, main$CI_lo, main$CI_hi,
    main$p_perm, main$q_family, main$q_all70, nrow(loo_smn), min(loo_smn$rho), max(loo_smn$rho)), "",
  sprintf(paste0(
    "**Figure 3.** Partial Spearman correlations (adjusted for age and aura status) between the mean Hurst ",
    "exponent of each of the seven networks and the nine plasma metabolites and the QUINA/KYNA ratio. ",
    "Asterisks and bold type mark q < 0.05 after Benjamini-Hochberg correction within the family ",
    "(L-KYN: 7 tests; all other measures: 63 tests); this was met for L-KYN with %s. ",
    "No test had q < 0.05 when all 70 tests were corrected together (smallest q = %.2f)."),
    paste(sig_n, collapse = ", "), min(tests$q_all70)), "",
  sprintf(paste0(
    "**Figure S1.** As Figure 2, restricted to the %d participants without aura; partial Spearman ",
    "correlation adjusted for age (rho = %.2f, 95%% BCa CI %.2f to %.2f, permutation p = %.4f)."),
    m$n, m$rho, m$CI_lo, m$CI_hi, m$p_perm))
writeLines(cap, "figures/captions.md")
message("figures written to figures/")
