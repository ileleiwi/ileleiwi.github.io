---
title: "Projects"
permalink: /projects/
layout: single
author_profile: true
---


Below are a few current and recent efforts. (Add images to `/assets/images/` and link out to slides/papers/data.)

### SIP/qSIP + Bayesian Hierarchical Models
I fit a robust Student-t mixed model on raw EAF with identity link, fixed effects for expression (centered), time since rewet (48h vs 168h), and moisture (50 vs 100), and random intercepts and expression slopes by MAG. Residual scale (σ) is modeled by rewet + moisture (heteroskedastic).

## Findings
  - 168h after rewet lowers mean EAF by ~0.03 (3 pp); 95% CrI ≈ [−0.04, −0.02].
  - Expression: effect ≈ 0; 95% CrI ≈ [−0.01, 0.02].
  - Interaction (expr × 168h): ≈ 0; little evidence of time-dependent change in the expression–EAF relationship.
  - Between-MAG variation is substantial: SD(intercept) ≈ 0.12; SD(expr slope) ≈ 0.06.
  - Variance shrinks at 100% moisture: sigma_moisture100 ≈ −1.72 → σ × exp(−1.72) ≈ 0.18 (≈ 82% reduction).
  - Heavy tails (ν ≈ 2.24) are appropriate given the zero mass.

<div style="display:flex; gap:1rem; align-items:flex-start; flex-wrap:wrap;">
  <div style="flex:1 1 380px;">
    <details><summary><strong>Model summary (brms)</strong></summary>
{% capture fit_summary %}{% include_relative assets/code/qsip-bhm/results/fit_t_lin-summary.txt %}{% endcapture %}
{% highlight text linenos %}{{ fit_summary }}{% endhighlight %}
</details>

<details><summary><strong>R script</strong></summary>
{% capture code_r %}{% include_relative assets/code/qsip-bhm/Final_Bayesian_Model.R %}{% endcapture %}
{% highlight r linenos %}{{ code_r }}{% endhighlight %}
</details>

  </div>
  <div style="flex:0 0 320px; max-width:320px; margin-left:auto;">
    <img src="/assets/code/qsip-bhm/figs/cond_effect_expr_by_rewet.png" alt="Conditional effects"
         style="width:100%; border-radius:12px; box-shadow:0 4px 12px rgba(0,0,0,.12); margin-bottom:12px;">
    <img src="/assets/code/qsip-bhm/figs/mcmc_intervals_fixed.png" alt="Fixed effects intervals"
         style="width:100%; border-radius:12px; box-shadow:0 4px 12px rgba(0,0,0,.12);">
  </div>
</div>


### Functional potential → MAG activity across redox treatments

We tested whether genome-encoded **functions/paths** predict which MAGs are **active** under four redox regimes (anoxic, low-frequency, high-frequency, oxic). The left panel ranks pathways by **likelihood-ratio χ²** importance; the right panels summarize **active MAG counts** and **mean gene copy per MAG** across treatments.

<div style="display:flex; gap:1rem; align-items:flex-start; flex-wrap:wrap;">
  <div style="flex:1 1 420px; min-width:320px;">

**Findings (edit to match your run):**
- Top predictors include **lactate utilization**, **arsenate reduction**, **siderophore-mediated iron acquisition**, **ETC complexes (I–IV)**, **sulfur metabolism**, and short-chain fatty-acid utilization pathways.
- Several functions show **treatment-specific enrichment** (e.g., anaerobic energy metabolism in **anoxic**; iron acquisition patterns shifting under **oxic**).
- Pathway copy number per MAG helps explain **which genomes become active** under each redox regime, beyond baseline abundance.

**Methods (1–2 sentences):**  
Binomial GLMs (active vs. inactive per MAG×treatment) with pathway indicators/copy numbers as predictors. Importance = **LR χ²** per feature; summaries are cross-tabs of active MAGs and mean gene counts.

  </div>

  <figure style="flex:0 0 420px; max-width:420px; margin-left:auto;">
    <img src="/assets/code/mag-logreg/figs/logistic_regression_barplot.pdf"
         alt="Functional predictors of MAG activity across redox treatments"
         style="width:100%; border-radius:12px; box-shadow:0 4px 12px rgba(0,0,0,.12);">
    <figcaption style="font-size:0.9em; opacity:0.85;">
      Functional predictors (LR χ²) and treatment-wise summaries of active MAGs and mean gene counts/MAG.
    </figcaption>
  </figure>
</div>

<details><summary><strong>Show analysis code</strong></summary>

{% capture code_r %}{% include_relative assets/code/mag-logreg/logistic_regression_dram_product_mags.R %}{% endcapture %}
{% highlight r linenos %}{{ code_r }}{% endhighlight %}

</details>

<p>
  <a class="btn btn--primary" href="/assets/code/mag-logreg/logistic_regression_dram_product_mags.R">Download the R script</a>
</p>


