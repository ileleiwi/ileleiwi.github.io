---
layout: splash
title: "Ikaia (Kai) Leleiwi"
author_profile: true
header:
  overlay_color: "#000"
  overlay_filter: "0.25"
  overlay_image: /assets/images/hero_soil.jpg
  caption: "Photo: your credit"
excerpt: "Soil microbiomes • metagenomics • SIP/qSIP • Bayesian hierarchical models • Machine Learning • HPC pipelines"

# data for rows lives in YAML (front matter)
intro:
  - excerpt: "Postdoctoral Researcher at LLNL specializing in machine learning methods for microbial ecology and multi-omics."

feature_row:
  - alt: "About"
    title: "About"
    excerpt: "Who I am and what I like to do."
    url: "/about/"
    btn_label: "About me"
    btn_class: "btn--primary"
  - alt: "Research Projects"
    title: "Research Projects"
    excerpt: "What I'm building and studying, from SIP/qSIP to genome-resolved metagenomics."
    url: "/projects/"
    btn_label: "Learn more"
    btn_class: "btn--primary"
  - alt: "Publications"
    title: "Publications"
    excerpt: "Peer-reviewed papers and preprints."
    url: "/publications/"
    btn_label: "See publications"
    btn_class: "btn--primary"
---

{% include author-profile.html %}

<style>
.page__hero--overlay .page__title,
.page__hero--overlay .page__lead{
  background: rgba(255,255,255,0.2);
  backdrop-filter: blur(6px);
  -webkit-backdrop-filter: blur(6px);
  border-radius: 12px;
  padding: 0.5em 1em;
  display: inline-block;
  box-shadow: 0 4px 12px rgba(0,0,0,0.3);
  color: #fff;
}
</style>

{%- comment -%} Optional: show the intro block centered {%- endcomment -%}
{% include feature_row id="intro" type="center" %}

{% include feature_row id="feature_row" %}
