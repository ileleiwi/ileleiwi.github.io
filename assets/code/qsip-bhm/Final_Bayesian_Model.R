
## ---------------------------
## Script: Final_Bayesian_Model.R
## ---------------------------

## Libraries ##
library(tidyverse)
library(cmdstanr)
library(brms)
library(emmeans)
library(bayesplot)
library(posterior)
library(yardstick)
library(patchwork)
     
set_cmdstan_path(cmdstanr::cmdstan_path())
options(brms.backend = "cmdstanr") 

# Output dirs 
out_dir <- "."
dir.create(file.path(out_dir, "figs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)

df <- read_tsv("data/metag_metat_eaf.tsv")

#prep data
df_prepped <- df |>
  filter(rewet %in% c("48h", "168h")) |>
  mutate(
    n = n(),
    EAF_sc = ifelse(
      mean_resampled_EAF_0 == 0, 0, #keep 0's at 0
      (mean_resampled_EAF_0 * ((n - 1) + 0.5 )) / n), #Smithson-Verkuilen formula, squeezes positive values inside (0,1)
    expr_cent = expression_ra_log10_cent,
    rewet = factor(rewet, levels = c("48h", "168h")),
    moisture = factor(moisture, levels = c("50", "100")),
    mag_id = factor(mag_id)
  ) |>
  select(mag_id, EAF_sc, expr_cent, moisture, rewet, mag_id)


#zero inflated beta model
#guide priors with this knowledge of the data
#about 70 % of MAG × sample observations are exactly zero; the positives never reach 1 
#(empirical max ≈ 0.68, and we don’t expect them to); RNA effects are small.

#initial fit
form_ziRE <- bf(
  EAF_sc ~ s(expr_cent, k = 5, by = rewet) + rewet + moisture + (1 + expr_cent | mag_id),
  phi    ~ 1 + rewet + moisture,
  zi     ~ expr_cent + rewet + moisture + (1 | mag_id)   
)

priors_ziRE <- c(
  # μ-part logit(0.2) = -1.39, qlogis(0.2)
  set_prior("normal(-1.39, 0.7)", class = "Intercept"),
  set_prior("normal(0, 0.5)",     class = "b"),
  set_prior("student_t(3, 0, 2.5)", class = "sd"),
  # φ-part (log link; intercept on log(precision))
  set_prior("normal(4, 0.5)",     class = "Intercept", dpar = "phi"),
  set_prior("normal(0, 0.3)",     class = "b",         dpar = "phi"),
  # zi-part (logit link; baseline P(Y=0)), logit(0.7) = 0.847, qlogis(0.7) 
  set_prior("normal(0.847, 0.5)", class = "Intercept", dpar = "zi"), # ~70% zeros in data 
  set_prior("normal(0, 1)",       class = "b",         dpar = "zi"),
  set_prior("student_t(3, 0, 2.5)", class = "sd",      dpar = "zi")  # RE SD in zi
)

fam <- zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit")

fit_zib_ziRE <- brm(
  formula  = form_ziRE,
  data     = df_prepped,
  family   = fam,
  prior    = priors_ziRE,
  chains   = 4, cores = 4, iter = 4000, warmup = 1000, seed = 123,
  control  = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = TRUE)
)

#prior predictions
fit_prior <- update(fit_zib_ziRE, sample_prior = "only")
y_rep_prior <- brms::posterior_predict(fit_prior, ndraws = 200)
y_rep_post  <- brms::posterior_predict(fit_zib_ziRE, ndraws = 200)

# Observed y
y_obs <- brms::get_y(fit_zib_ziRE)
if (is.matrix(y_obs)) y_obs <- y_obs[, 1]

# Pull the model data in the same row order as y_obs / yrep columns
dat <- fit_zib_ziRE$data
stopifnot(nrow(dat) == length(y_obs))

dat <- dat |>
  mutate(
    rewet = factor(rewet),
    moisture = factor(moisture)
  )

# Row indices per facet group
grp <- dat |>
  mutate(i = row_number()) |>
  group_by(rewet, moisture) |>
  summarise(idx = list(i), .groups = "drop")

zero_rate_draws <- function(yrep, idx) {
  rowMeans(yrep[, idx, drop = FALSE] == 0)
}

sample_group_mat <- function(yrep, idx, n = 30000) {
  sub <- yrep[, idx, drop = FALSE]
  N <- length(sub)
  sub[sample.int(N, min(n, N))]
}


df_zero_prior <- bind_rows(lapply(seq_len(nrow(grp)), function(g) {
  idx <- grp$idx[[g]]
  tibble(
    value = zero_rate_draws(y_rep_prior, idx),
    source = "Prior predictive",
    panel  = "Fraction of zeros  Pr(EAF_sc = 0)",
    rewet  = grp$rewet[[g]],
    moisture = grp$moisture[[g]]
  )
}))

df_zero_post <- bind_rows(lapply(seq_len(nrow(grp)), function(g) {
  idx <- grp$idx[[g]]
  tibble(
    value = zero_rate_draws(y_rep_post, idx),
    source = "Posterior predictive",
    panel  = "Fraction of zeros  Pr(EAF_sc = 0)",
    rewet  = grp$rewet[[g]],
    moisture = grp$moisture[[g]]
  )
}))

df_zero <- bind_rows(df_zero_prior, df_zero_post)

# Observed zero fraction per group for vertical reference lines
df_zero_vline <- grp |>
  transmute(
    panel = "Fraction of zeros  Pr(EAF_sc = 0)",
    rewet, moisture,
    xint = vapply(idx, function(ii) mean(y_obs[ii] == 0), numeric(1))
  )


n_per_group <- 30000

df_pos_obs <- bind_rows(lapply(seq_len(nrow(grp)), function(g) {
  idx <- grp$idx[[g]]
  v <- y_obs[idx]
  v <- v[v > 0]
  if (length(v) == 0) return(NULL)
  
  tibble(
    value = v,
    source = "Observed",
    panel  = "Non-zero values  (distribution on (0,1))",
    rewet  = grp$rewet[[g]],
    moisture = grp$moisture[[g]]
  )
}))

df_pos_prior <- bind_rows(lapply(seq_len(nrow(grp)), function(g) {
  idx <- grp$idx[[g]]
  v <- sample_group_mat(y_rep_prior, idx, n = n_per_group)
  v <- v[v > 0]
  if (length(v) == 0) return(NULL)
  
  tibble(
    value = v,
    source = "Prior predictive",
    panel  = "Non-zero values  (distribution on (0,1))",
    rewet  = grp$rewet[[g]],
    moisture = grp$moisture[[g]]
  )
}))

df_pos_post <- bind_rows(lapply(seq_len(nrow(grp)), function(g) {
  idx <- grp$idx[[g]]
  v <- sample_group_mat(y_rep_post, idx, n = n_per_group)
  v <- v[v > 0]
  if (length(v) == 0) return(NULL)
  
  tibble(
    value = v,
    source = "Posterior predictive",
    panel  = "Non-zero values  (distribution on (0,1))",
    rewet  = grp$rewet[[g]],
    moisture = grp$moisture[[g]]
  )
}))

df_pos <- bind_rows(df_pos_obs, df_pos_prior, df_pos_post)


df_plot <- bind_rows(df_zero, df_pos) |>
  mutate(
    source = factor(source, levels = c("Observed", "Prior predictive", "Posterior predictive")),
    panel  = factor(panel, levels = c("Fraction of zeros  Pr(EAF_sc = 0)",
                                      "Non-zero values  (distribution on (0,1))")),
    moisture = factor(moisture, levels = levels(dat$moisture))
  )

pdf(file.path(out_dir, "figs", "fit_zip_ziRE_observed_vs_prior_vs_post.pdf"), height = 10, width = 10)
ggplot(df_plot, aes(x = value, color = source, linetype = source)) +
  geom_density(linewidth = 1, adjust = 1) +
  geom_vline(
    data = df_zero_vline,
    aes(xintercept = xint),
    inherit.aes = FALSE,
    linetype = 2
  ) +
  facet_grid(panel + moisture ~ rewet, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Observed vs prior predictive vs posterior predictive (EAF_sc)\nFaceted by rewet and moisture bin",
    x = NULL,
    y = "Density"
  ) +
  theme_minimal()
dev.off()

## zero inflated model posterior predicitve values misaligned to observed values in each time x moisture condition pair,
## decided to add rewet*moisture interaction term and refit

priors_ziRE_int <- c(
  # μ-part logit(0.2) = -1.39, qlogis(0.2)
  set_prior("normal(-1.39, 0.7)", class = "Intercept"),
  set_prior("normal(0, 0.5)",     class = "b"),
  set_prior("student_t(3, 0, 2.5)", class = "sd"),
  # φ-part (log link; intercept on log(precision))
  set_prior("normal(4, 0.5)",     class = "Intercept", dpar = "phi"),
  set_prior("normal(0, 0.3)",     class = "b",         dpar = "phi"),
  # zi-part (logit link; baseline P(Y=0)), logit(0.7) = 0.847, qlogis(0.7) ,
  set_prior("normal(0.847, 0.5)", class = "Intercept", dpar = "zi"), # ~70% zeros in data
  set_prior("normal(0, 1)",       class = "b",         dpar = "zi"),
  set_prior("student_t(3, 0, 2.5)", class = "sd",      dpar = "zi")  # RE SD in zi
)

form_ziRE_int <- bf(
  EAF_sc ~ s(expr_cent, k = 5, by = rewet) + rewet + moisture + (1 + expr_cent | mag_id),
  phi    ~ 1 + rewet + moisture,
  zi     ~ expr_cent + rewet * moisture + (1 | mag_id)
)

fam <- zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit")

fit_zib_ziRE_int <- update(
  fit_zib_ziRE,
  formula. = form_ziRE_int,
  prior    = priors_ziRE_int, 
  newdata  = df_prepped,
  chains   = 4, cores = 4, iter = 4000, warmup = 1000, seed = 123,
  control  = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = TRUE),
  recompile = TRUE
)

# Save summary text ----
capture.output(
  { print(summary(fit_zib_ziRE_int)) },
  file = file.path(out_dir, "results", "fit_zib-summary.txt")
)

fit_prior_int <- update(fit_zib_ziRE_int, sample_prior = "only")
y_rep_prior_int <- brms::posterior_predict(fit_prior_int, ndraws = 200)
y_rep_post_int  <- brms::posterior_predict(fit_zib_ziRE_int, ndraws = 200)
# Observed y
y_obs <- brms::get_y(fit_zib_ziRE_int)
if (is.matrix(y_obs)) y_obs <- y_obs[, 1]

# Pull the model data in the same row order as y_obs / yrep columns
dat <- fit_zib_ziRE_int$data
stopifnot(nrow(dat) == length(y_obs))

dat <- dat |>
  mutate(
    rewet = factor(rewet),
    moisture = factor(moisture)
  )

# Row indices per facet group
grp <- dat |>
  mutate(i = row_number()) |>
  group_by(rewet, moisture) |>
  summarise(idx = list(i), .groups = "drop")


df_zero_prior_int <- bind_rows(lapply(seq_len(nrow(grp)), function(g) {
  idx <- grp$idx[[g]]
  tibble(
    value = zero_rate_draws(y_rep_prior_int, idx),
    source = "Prior predictive",
    panel  = "Fraction of zeros  Pr(EAF_sc = 0)",
    rewet  = grp$rewet[[g]],
    moisture = grp$moisture[[g]]
  )
}))

df_zero_post_int <- bind_rows(lapply(seq_len(nrow(grp)), function(g) {
  idx <- grp$idx[[g]]
  tibble(
    value = zero_rate_draws(y_rep_post_int, idx),
    source = "Posterior predictive",
    panel  = "Fraction of zeros  Pr(EAF_sc = 0)",
    rewet  = grp$rewet[[g]],
    moisture = grp$moisture[[g]]
  )
}))

df_zero_int <- bind_rows(df_zero_prior_int, df_zero_post_int)

# Observed zero fraction per group for vertical reference lines
df_zero_vline <- grp |>
  transmute(
    panel = "Fraction of zeros  Pr(EAF_sc = 0)",
    rewet, moisture,
    xint = vapply(idx, function(ii) mean(y_obs[ii] == 0), numeric(1))
  )

n_per_group <- 30000

df_pos_obs_int <- bind_rows(lapply(seq_len(nrow(grp)), function(g) {
  idx <- grp$idx[[g]]
  v <- y_obs[idx]
  v <- v[v > 0]
  if (length(v) == 0) return(NULL)
  
  tibble(
    value = v,
    source = "Observed",
    panel  = "Non-zero values  (distribution on (0,1))",
    rewet  = grp$rewet[[g]],
    moisture = grp$moisture[[g]]
  )
}))

df_pos_prior_int <- bind_rows(lapply(seq_len(nrow(grp)), function(g) {
  idx <- grp$idx[[g]]
  v <- sample_group_mat(y_rep_prior_int, idx, n = n_per_group)
  v <- v[v > 0]
  if (length(v) == 0) return(NULL)
  
  tibble(
    value = v,
    source = "Prior predictive",
    panel  = "Non-zero values  (distribution on (0,1))",
    rewet  = grp$rewet[[g]],
    moisture = grp$moisture[[g]]
  )
}))

df_pos_post_int <- bind_rows(lapply(seq_len(nrow(grp)), function(g) {
  idx <- grp$idx[[g]]
  v <- sample_group_mat(y_rep_post_int, idx, n = n_per_group)
  v <- v[v > 0]
  if (length(v) == 0) return(NULL)
  
  tibble(
    value = v,
    source = "Posterior predictive",
    panel  = "Non-zero values  (distribution on (0,1))",
    rewet  = grp$rewet[[g]],
    moisture = grp$moisture[[g]]
  )
}))

df_pos_int <- bind_rows(df_pos_obs_int, df_pos_prior_int, df_pos_post_int)

df_plot_int <- bind_rows(df_zero_int, df_pos_int) |>
  mutate(
    source = factor(source, levels = c("Observed", "Prior predictive", "Posterior predictive")),
    panel  = factor(panel, levels = c("Fraction of zeros  Pr(EAF_sc = 0)",
                                      "Non-zero values  (distribution on (0,1))")),
    moisture = factor(moisture, levels = levels(dat$moisture))
  )
pdf(file.path(out_dir, "figs", "fit_zip_ziRE_int_observed_vs_prior_vs_post.pdf"), height = 10, width = 10)
ggplot(df_plot_int, aes(x = value, color = source, linetype = source)) +
  geom_density(linewidth = 1, adjust = 1) +
  geom_vline(
    data = df_zero_vline,
    aes(xintercept = xint),
    inherit.aes = FALSE,
    linetype = 2
  ) +
  facet_grid(panel + moisture ~ rewet, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Observed vs prior predictive vs posterior predictive (EAF_sc)\nFaceted by rewet and moisture bin",
    x = NULL,
    y = "Density"
  ) +
  theme_minimal()
dev.off()

# Expression sequence for smooth curves (use central range to avoid extreme tails)
expr_seq <- seq(
  quantile(df_prepped$expr_cent, 0.02, na.rm = TRUE),
  quantile(df_prepped$expr_cent, 0.98, na.rm = TRUE),
  length.out = 200
)

# Grid for predictions
grid <- expand_grid(
  expr_cent = expr_seq,
  moisture  = levels(df_prepped$moisture),
  rewet     = levels(df_prepped$rewet)
) |>
  arrange(expr_cent, moisture, rewet) |>
  mutate(mag_id = df_prepped$mag_id[1])   # placeholder; ignored if re_formula = NA


# Cell weights from your data
w <- df_prepped |>
  count(moisture, rewet) |>
  arrange(moisture, rewet) |>
  mutate(w = n / sum(n))

# Posterior draws of zi probability (Pr(y=0)) on response scale
nd <- 2000  # increase for smoother ribbons; decrease if slow
p0_draws <- posterior_linpred(
  fit_zib_ziRE_int,
  dpar = "zi",
  newdata = grid,
  re_formula = NA,     # population-level (no MAG random effects)
  transform = TRUE,
  ndraws = nd
)

# Weighted average across the 4 treatment cells, within each expression value
n_expr  <- length(expr_seq)
n_combo <- nrow(w)               # should be 4
w_vec   <- w$w                   # order matches grid within each expr because of arrange()

p0_marg_draws <- sapply(seq_len(n_expr), function(i) {
  idx <- ((i - 1) * n_combo + 1):(i * n_combo)
  as.vector(p0_draws[, idx, drop = FALSE] %*% w_vec)
})
# p0_marg_draws is ndraws x n_expr

df_p0 <- tibble(
  expr_cent = expr_seq,
  estimate  = colMeans(p0_marg_draws),
  lwr       = apply(p0_marg_draws, 2, quantile, 0.025),
  upr       = apply(p0_marg_draws, 2, quantile, 0.975)
)

p_left <- ggplot(df_p0, aes(x = expr_cent, y = estimate)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25) +
  geom_line(linewidth = 1) +
  labs(
    title = "Overall transcription",
    x = "Centered expression (log10)",
    y = "Pr(EAF_sc = 0)"
  ) +
  theme_minimal(base_size = 14)


epred_draws <- posterior_epred(
  fit_zib_ziRE_int,
  newdata = grid,
  re_formula = NA,   # population-level
  ndraws = nd
)

df_epred <- bind_cols(
  grid,
  tibble(
    estimate = colMeans(epred_draws),
    lwr      = apply(epred_draws, 2, quantile, 0.025),
    upr      = apply(epred_draws, 2, quantile, 0.975)
  )
)

p_right <- ggplot(
  df_epred,
  aes(
    x = expr_cent,
    y = estimate,
    color = moisture,
    linetype = rewet,
    group = interaction(moisture, rewet)
  )
) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = moisture),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 1) +
  labs(
    x = "Centered expression (log10)",
    y = "Predicted average EAF_sc",
    color = "Moisture",
    fill  = "Moisture",
    linetype = "Rewet"
  ) +
  theme_minimal(base_size = 14)

pdf("figs/zib_ziRE_int_curves.pdf", height = 10, width = 15)
(p_left | p_right) +
  plot_layout(widths = c(1, 1.2))
dev.off()