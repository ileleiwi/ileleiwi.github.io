---
title: "Research Projects"
permalink: /projects/
toc: false
---

Below are a few current and recent efforts. (Add images to `/assets/images/` and link out to slides/papers/data.)

### SIP/qSIP + Bayesian Hierarchical Models
Here I fit a robust **t** mixed model with fixed effects accounting for moisture and time and random intercepts/slopes per MAG.
- Predictors: expression (centered), time since rewet (48h vs 168h), moisture; heteroskedastic σ by rewet + moisture.
- Identity link handled with t family and informative priors.
- Key findings: 
    - Substantial vairation in EAF between MAGs. 168h post rewet lowers mean EAF by ~3. Exp[ression isn't a reliable predictor of EAF in this model.

**Findings (qSIP + BHM).** Baseline EAF ≈ 0.07 (95% CI 0.04–0.10). Rewet 168h lowers mean EAF by ~0.03 (95% CI −0.04 to −0.02). Expression effect ≈ 0 (−0.01 to 0.02), and the 168h interaction is ~0. Residual variance drops ~82% at 100% moisture (σ × exp(−1.72) ≈ 0.18). Consider a zero-one inflated Beta to capture the 0-spike explicitly.

<div style="display:flex; gap:1rem; align-items:flex-start; flex-wrap:wrap;">
  <div style="flex:1 1 380px;">
    <details><summary><strong>Model summary</strong></summary>
{% capture fit_summary %}{% include_relative ../assets/code/qsip-bhm/results/fit_t_lin-summary.txt %}{% endcapture %}
<pre>{{ fit_summary }}</pre>
    </details>
    <details><summary><strong>R script</strong></summary>
{% capture code_r %}{% include_relative ../assets/code/qsip-bhm/Final_Bayesian_Model.R %}{% endcapture %}
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


