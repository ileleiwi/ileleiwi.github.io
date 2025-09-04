---
title: "Projects"
permalink: /projects/
layout: single
author_profile: true
---


Below are a few current and recent efforts. (Add images to `/assets/images/` and link out to slides/papers/data.)

### SIP/qSIP + Bayesian Hierarchical Models
I fit a robust Student-t mixed model on raw EAF with identity link, fixed effects for expression (centered), time since rewet (48h vs 168h), and moisture (50 vs 100), and random intercepts and expression slopes by MAG. Residual scale (σ) is modeled by rewet + moisture (heteroskedastic).

Findings.

168h after rewet lowers mean EAF by ~0.03 (3 pp); 95% CrI ≈ [−0.04, −0.02].

Expression: effect ≈ 0; 95% CrI ≈ [−0.01, 0.02].

Interaction (expr × 168h): ≈ 0; little evidence of time-dependent change in the expression–EAF relationship.

Between-MAG variation is substantial: SD(intercept) ≈ 0.12; SD(expr slope) ≈ 0.06.

Variance shrinks at 100% moisture: sigma_moisture100 ≈ −1.72 → σ × exp(−1.72) ≈ 0.18 (≈ 82% reduction).

Heavy tails (ν ≈ 2.24) are appropriate given the zero mass.

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


### Genome‑resolved metagenomics + metatranscriptomics (MAG‑centric)
How you normalize (GeTMM; length/size) and how you analyze activity vs abundance.  
**Links:** [methods](), [datasets](), [workflow]().

### Forest floor nitrification under warming
Key findings and why it matters for N cycling.  
**Links:** [paper](), [data](), [press]().


