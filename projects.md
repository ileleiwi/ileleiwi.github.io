---
title: "Projects"
permalink: /projects/
layout: single
author_profile: true
---


Below are a few current and recent efforts. 
### Using Bayesian hierarchical modeling to answer if qSIP activity correlates with transcription? 
Using R packages cmdstanr and brms, I fit a zero inflated beta model to model qSIP activity (Excess Atom Fraction, EAF) with fixed effects for expression (transcription), time since rewet (48h vs 168h), and moisture (50% vs 100%), and random intercepts and expression slopes by MAG. 

By measuring the density of DNA from microbial communities that have been treated with H<sub>2</sub><sup>18</sup>O relative to control densities we can determine if bacteria were actively dividing and incorporating heavy oxygen in their genomes. This is a technique called quantitative stable isotope probing (qSIP), and the resulting metric of activity is EAF.

## Model structure
<div style="flex:0 0 320px; max-width:320px; margin-left:auto;">
    <img src="/assets/code/qsip-bhm/figs/model_structure_slide6.png" alt="Conditional effects"
         style="width:100%; border-radius:12px; box-shadow:0 4px 12px rgba(0,0,0,.12); margin-bottom:12px;">
  </div>

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
<div markdown="1" style="flex:1 1 420px; min-width:320px;">

**Findings:**
- Top predictors include **lactate utilization**, **arsenate reduction**, **siderophore-mediated iron acquisition**, **ETC complexes (I–IV)**, **sulfur metabolism**, and short-chain fatty-acid utilization pathways.
- Iron acquisition is potentially important for active bacteria in anoxic conditions
- Activity in anoxic soils may also depend on the ability of bacteria to utilize various short-chain fatty acids as carbon sources and sulfur compounds as terminal electron acceptors for anaerobic respiration.

**Methods:**  
Binomial GLMs (active vs. inactive per MAG×treatment) with pathway indicators/copy numbers as predictors. Importance = **LR χ²** per feature; summaries are cross-tabs of active MAGs and mean gene counts.
</div>

  <figure style="flex:0 0 420px; max-width:420px; margin-left:auto;">
    <img src="/assets/code/mag-logreg/figs/logistic_regression_barplot.svg"
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


