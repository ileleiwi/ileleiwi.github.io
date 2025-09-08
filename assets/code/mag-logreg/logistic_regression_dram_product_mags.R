## ---------------------------
##
## Script: logistic_regression_dram_product_mags
##
## Purpose: classify activity using logistic regression
##          
## Author: Ikaia Leleiwi
##
## Date Created: November 12th, 2024
##
## Copyright (c) Ikaia Leleiwi, 2024
## Email: leleiwi1@llnl.gov
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory

setwd(paste0("/Users/leleiwi1/Desktop/LLNL_postdoc/Soil_Microbiome_SFA/",
             "Great_Redox_Experiment/MHM2_analysis/Rdirectory_MHM2"))

## ---------------------------

## Libraries ##
library(tidyverse)
library(readxl)
library(compositions)
library(tidymodels)
library(themis)
library(vip)
library(usemodels)
library(furrr)
library(Boruta)
library(recipeselectors)
library(openxlsx)



##Functions##
relabund <- function(df, columns = c(NA)) 
  #takes feature table and calculates relative abundance of columns, omitting NA's
  #considers only columns listed in columns argument (character vector of column names). 
  #columns (default) = all columns
{
  if (NA %in% columns){
    df <- sweep(df, 2, colSums(df, na.rm = TRUE), '/')
  }
  else {
    df_relabund <- df %>% select(all_of(columns)) %>%
      sweep(., 2, colSums(., na.rm = TRUE), '/')
    df[,columns] <- df_relabund
  }
  
  return(df)
}

#interpretation functions
baseline_odds <- function(intercept){
  #takes intercept from summary(model)
  #outputs baseline odds
  out <- exp(intercept)
  return(out)
}

baseline_prob <- function(intercept){
  #takes intercept from summary(model)
  #outputs baseline probability
  out <- exp(intercept)/(1+exp(intercept))
  return(out)
}

adjusted_var_odds <- function(base_odds, var_or){
  #takes baseline odds and variable odds ratio
  #outputs adjusted odds 
  out <- base_odds*var_or
  return(out)
}

prob_var <- function(intercept, var_or){
  #takes odds ratio of a variable and intercept
  #outputs interpretation statement
  #interpretation - when predictor is TRUE, the probability of outcome TRUE is %, 
  #                 this represents a % increase relative to baseline (%)
  
  base_prob <- baseline_prob(intercept = intercept)
  adj_odds <- adjusted_var_odds(base_odds = baseline_odds(intercept), var_or = var_or)
  prob_adj <- (adj_odds/(1+adj_odds))*100
  prob_change <- (prob_adj - base_prob*100)
  direction_string <- ifelse(prob_change > 0, "increase", "decrease")
  direction_string_plur <- ifelse(prob_change > 0, "increases", "decreases")
  rel_change <- ((prob_adj - (base_prob*100))/(base_prob*100))*100
  out_string <- str_glue("When predictor is TRUE, the probability of outcome being TRUE is {round(prob_adj, 3)}%, ",
                         "this represents an absolute {direction_string} of {round(prob_change,3)} percentage points",
                         " compared to the baseline probability of ({round(base_prob*100,3)}%). ",
                         "this corresponds to a {round(rel_change,3)}% relative {direction_string} in the probability of outcome being TRUE compared",
                         " to the baseline probability")
  return(out_string)
}

#colors
Bulk = "black"
C12 = "#4f94cd"
C13 = "#8B0000"
HF = "#e00641"
LF = "#ffa91f"
Oxic = "#4d45a4"
Anoxic = "#6ea458"
Field = "#aeadad"


## Data ##

drep <- read_csv("data/drep_data_tables/Wdb.csv") |>
  mutate(genome = str_remove(genome, ".fa"))

gtdb_b <- read_tsv("data/gtdb_v2.1.1_r214/gtdbtk.bac120.summary.tsv") |>
  select(genome = user_genome, classification) |>
  filter(genome %in% drep$genome)

lowest_gtdb <- gtdb_b |>
  mutate(classification = str_remove_all(classification, ".__")) |>
  separate(classification,
           into = c("domain", "phylum", "class", "order", 
                    "family", 'genus', "species"),
           sep = ";") |>
  rowid_to_column(var = "num") |> #add indexing column
  mutate(lowest_gtdb = case_when(phylum == "" ~ paste(domain, "domain"),
                                 class == "" ~ paste(phylum, "phylum"),
                                 order == "" ~ paste(class, "class"),
                                 family == "" ~ paste(order, "order"),
                                 genus == "" ~ paste(family, "family"),
                                 species == "" ~ paste(genus, "genus"),
                                 TRUE ~ species),
         lowest_gtdb = paste(lowest_gtdb, num, sep = "_")) |>
  select(-num)


metadata <- read_xlsx("data/metadata/GRE-info-Jose.xlsx") |>
  slice_head(n=-1) |>
  select(sample = Sample,
         redox = data,
         isotope = data2,
         sample_type = data3)

product <- read_tsv("data/distillate_all_mags_dram1.5/product.tsv") |> 
  right_join(lowest_gtdb,
             by = "genome") |>
  select(genome, lowest_gtdb, domain, phylum, class, order, family, genus, species,
         everything()) |>
  mutate(I = ifelse(`Complex I: NAD(P)H:quinone oxidoreductase, chloroplasts and cyanobacteria` > 0 |
                           `Complex I: NADH dehydrogenase (ubiquinone) 1 alpha subcomplex` > 0 |
                           `Complex I: NADH:quinone oxidoreductase, prokaryotes` > 0, TRUE, FALSE),
         II = ifelse(`Complex II: Fumarate reductase, prokaryotes` > 0 |
                       `Complex II: Succinate dehydrogenase (ubiquinone)` > 0 |
                       `Complex II: Succinate dehydrogenase, prokaryotes` > 0, TRUE, FALSE),
         III = ifelse(`Complex III: Cytochrome bc1 complex` > 0 |
                        `Complex III: Cytochrome bc1 complex respiratory unit` > 0 |
                        `Complex III: Cytochrome bd ubiquinol oxidase` > 0, TRUE, FALSE),
         IV = ifelse(`Complex IV High affinity: Cytochrome bd ubiquinol oxidase` > 0 |
                       `Complex IV High affinity: Cytochrome c oxidase, cbb3-type` > 0 |
                       `Complex IV Low affinity: Cytochrome aa3-600 menaquinol oxidase` > 0 |
                       `Complex IV Low affinity: Cytochrome c oxidase` > 0 |
                       `Complex IV Low affinity: Cytochrome c oxidase, prokaryotes` > 0 |
                       `Complex IV Low affinity: Cytochrome o ubiquinol oxidase` > 0, TRUE, FALSE),
         V = ifelse(`Complex V: F-type ATPase, eukaryotes` > 0 |
                      `Complex V: F-type ATPase, prokaryotes and chloroplasts` > 0 |
                      `Complex V: V-type ATPase, eukaryotes` > 0 |
                      `Complex V: V/A-type ATPase, prokaryotes` > 0, TRUE, FALSE),
         ETC = ifelse(I & II & III & IV & V, TRUE, FALSE), 
         butyrate = as.logical(rowSums(across(contains("Butyrate, pt")))),
         propionate = as.logical(rowSums(across(contains("Propionate, pt")))),
         acetate = as.logical(rowSums(across(contains("acetate, pt")))),
         lactate = as.logical(rowSums(across(contains("lactate ")))),
         acetylCoA = as.logical(rowSums(across(contains("pyruvate => acetyl")))),
         arsenate_reduction = as.logical(rowSums(across(contains("arsenate reduction, pt")))),
         sulfur = as.logical(rowSums(across(contains("Sulfur metabolism")))),
         nitrogen = as.logical(rowSums(across(contains("Nitrogen metabolism")))))

product_long <- product |>
  pivot_longer(cols = -c(genome, lowest_gtdb, domain, phylum, class, order, family, genus, species),
               names_to = "metabolism",
               values_to = "value") |>
  mutate(funct = case_when(metabolism %in% c("Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate",
                                             "Pentose phosphate pathway (Pentose phosphate cycle)",
                                             "Entner-Doudoroff pathway, glucose-6P => glyceraldehyde-3P + pyruvate",
                                             "Citrate cycle (TCA cycle, Krebs cycle)",
                                             "Glyoxylate cycle",
                                             "Reductive pentose phosphate cycle (Calvin cycle)",
                                             "Reductive citrate cycle (Arnon-Buchanan cycle)",
                                             "Dicarboxylate-hydroxybutyrate cycle",
                                             "Hydroxypropionate-hydroxybutylate cycle",
                                             "3-Hydroxypropionate bi-cycle",
                                             "Reductive acetyl-CoA pathway (Wood-Ljungdahl pathway)",
                                             "Acetyl-CoA pathway, CO2 => acetyl-CoA",
                                             "Methanogenesis, CO2 => methane") ~ "Module",
                           str_detect(metabolism, "Complex I:") ~ "I",
                           str_detect(metabolism, "Complex II:") ~ "II",
                           str_detect(metabolism, "Complex III:") ~ "III",
                           str_detect(metabolism, "Complex IV High affinity:") ~ "IV High affinity",
                           str_detect(metabolism, "Complex IV Low affinity:") ~ "IV Low affinity",
                           str_detect(metabolism, "Complex V:") ~ "V",
                           T ~ str_extract(metabolism, ".+:")),
         funct = str_remove(funct, ":"),
         ETC = ifelse(funct %in% c("Module", "I", "II", "III", 
                                   "IV High affinity", "IV Low affinity", 
                                   "V"), "ETC Complex", "Other"),
         funct = factor(funct, levels = c("Module", "I", "II", "III", 
                                          "IV High affinity", "IV Low affinity", 
                                          "V", "CAZy", "Methanogenesis and methanotrophy",
                                          "Nitrogen metabolism", "Other Reductases",
                                          "Photosynthesis", "SCFA and alcohol conversions",
                                          "Sulfur metabolism"))) 
  

funct_remove <- product_long |>
  filter(funct %in% c("Module", "Methanogenesis and methanotrophy") |
           ETC == "ETC Complex") |>
  pull(metabolism) |>
  unique()

gsf <- read_tsv("/Users/leleiwi1/Desktop/LLNL_postdoc/Side_projects/DRAM/DRAM_forms/genome_summary_form.tsv")
fhf <- read_tsv("/Users/leleiwi1/Desktop/LLNL_postdoc/Side_projects/DRAM/DRAM_forms/function_heatmap_form.tsv")


rename_cols <- function(x){
  
  x <- str_remove(x, "GRE\\.bulkMG\\.")
  x <- str_remove(x, "GRE\\.SIPMG\\.")
  x <- str_remove(x, "_Mean")
  
  return(x)
}

#Fegenie
file_list <- list.files(path = "data/FeGenie", pattern = ".csv", full.names = T)
fegenie_data <- lapply(file_list, read_csv)
names(fegenie_data) <- gsub("-" ,"_", gsub("\\.csv$", "", basename(file_list)))
list2env(fegenie_data, envir = .GlobalEnv)

FeGenie_geneSummary <- FeGenie_geneSummary[-1,]

FeGenie_heatmap_data <- FeGenie_heatmap_data |>
  rename(category = X) |>
  pivot_longer(cols = -category,
               values_to = "count",
               names_to = "fasta") |>
  mutate(fasta = str_remove(fasta, "\\.fa")) |>
  pivot_wider(names_from = category, 
              values_from = count,
              values_fill = 0)

fe_join <- FeGenie_heatmap_data |>
  drop_na() |>
  rowwise() |>
  mutate(iron_oxidation = sum(iron_oxidation, possible_iron_oxidation_and_possible_iron_reduction),
         iron_reduction = sum(iron_reduction, possible_iron_oxidation_and_possible_iron_reduction, probable_iron_reduction),
         siderophore = sum(`iron_aquisition-siderophore_synthesis`, 
                           `iron_aquisition-siderophore_transport`, 
                           `iron_aquisition-siderophore_transport_potential`)) |>
  select(fasta, iron_oxidation, iron_reduction, siderophore)



rabund <- read_tsv("data/drep_mags_relabund_10cov_95id.tsv") |>
  rename_with(rename_cols) |>
  #filter(if_any(!matches("bins"), ~ . != 0)) |>
  rename(genome = bins)

enriched <- read_tsv("data/enriched_mags.tsv") |>
  rowwise() |>
  mutate(total = LF+HF+AN+OX) 


#make ml dataframe
rabund_by_redox <- rabund |>
  column_to_rownames(var = "genome") |>
  clr() |>
  as.data.frame() |>
  rownames_to_column(var = "genome") |>
  pivot_longer(cols = -genome,
               names_to = "sample",
               values_to = "relative_abundance") |>
  left_join(metadata) |>
  group_by(genome, redox) |>
  summarise(mean_relabund = mean(relative_abundance)) |>
  pivot_wider(names_from = "redox",
              values_from = "mean_relabund") |>
  rename(High_Freq = HF, Low_Freq = LF)

df <- product |>
  left_join(enriched, by = "genome") |>
  left_join(fe_join, by = c("genome" = "fasta")) |>
  select(-c(total, genome, lowest_gtdb, domain, phylum, class, order, family,
            genus, species, any_of(funct_remove), starts_with("CAZy"), I, II, 
            III, IV, V, contains("arsenate reduction, pt"), starts_with("SCFA and alcohol"),
            starts_with("Sulfur metabolism"), starts_with("Nitrogen metabolism"))) |>
  mutate(across(c(LF, OX, AN, HF), ~ifelse(is.na(.), 0, .)),
         across(c(LF, OX, AN, HF), ~as.logical(.))) |>
  select(where(~any(. != 0))) |>
  mutate(iron_oxidation = ifelse(is.na(iron_oxidation), 0, iron_oxidation),
         iron_reduction = ifelse(is.na(iron_reduction), 0, iron_reduction),
         siderophore = ifelse(is.na(siderophore), 0, siderophore))

sum(is.na(df))
sum(is.na(fe_join))



df |>
  pivot_longer(cols = c(AN, LF, HF, OX),
               names_to = "redox",
               values_to = "mean_relabund") |>
  ggplot(aes(x = mean_relabund)) +
  geom_histogram(stat = "count") +
  facet_wrap(~redox)

##Histogram and Bar chart plot function
hist_all <- function(df, type=""){
  # List of packages for session
  .packages = c("ggplot2", "cowplot")
  # Install CRAN packages (if not already installed)
  .inst <- .packages %in% installed.packages()
  if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
  # Load packages into session 
  lapply(.packages, require, character.only=TRUE)
  
  # Filter dataframe to only numeric columns
  df_n <- df %>%
    select_if(is.numeric)
  
  # Build numeric variable plot list
  build_plots <- function(dn){
    nlist <- list()
    for(cat in colnames(dn)){
      p <- dn %>%
        ggplot(aes(x = .data[[cat]])) +
        geom_histogram(bins = 100) +
        labs(y = "count",
             title = cat) +
        theme_bw() +
        theme(axis.title.x = element_blank())
      nlist[[length(nlist)+1]] <- p
    }
    names(nlist) <- colnames(dn)
    return(nlist)
  }
  
  plist <- build_plots(df_n)
  
  # Arrange histograms in one plot
  out_plot <- do.call("plot_grid", plist)
  
  return(out_plot)
}

hist_all(df)

convert_to_factor <- function(column) {
  if (is.logical(column)) {
    return(as.factor(column))
  } else if (is.numeric(column) && length(unique(column)) <= 7) {
    return(factor(column, levels = sort(unique(column)), ordered = TRUE))
  }else{
    return(column)
  }
}

convert_logic_to_factor <- function(column) {
  if (is.logical(column)) {
    return(as.factor(column))
  } else{
    return(column)
  }
}

df <- df |>
  mutate(across(everything(), ~ convert_to_factor(.))) 

sum(is.na(df))


for(i in 1:4){
  redox_cols <- c("AN", "OX", "LF", 'HF')
  df_names <- paste0("df", "_", redox_cols)
  current_col <- redox_cols[i]
  
  df_out <- df |>
    select(-setdiff(redox_cols, current_col)) |>
    janitor::clean_names()
  
  assign(df_names[i], df_out)
}

#fit logistic regression
m_an <- glm(formula = an ~ ., data = df_AN, family = binomial)
m_ox <- glm(formula = ox ~ ., data = df_OX, family = binomial)
m_lf <- glm(formula = lf ~ ., data = df_LF, family = binomial)
m_hf <- glm(formula = hf ~ ., data = df_HF, family = binomial)

#model summaries
summary(m_an) #intercept -1.87260, 
prob_var(intercept = -1.87260, var_or = 1.71) #1.71 is odds ratio from gtsummary below for fullETC
prob_var(intercept = -1.87260, var_or = 2.05) #lactate metabolism

summary(m_ox) #intercept -1.6450248
prob_var(intercept = -1.6450248, var_or = 0.35) #sulfur
summary(m_lf) #intercept -2.21532
prob_var(intercept = -2.21532, var_or = 0.34) #sulfur
summary(m_hf) #intercept -1.941279
prob_var(intercept = -1.941279, var_or = 0.47) #sulfur

#check for nonlinearity
# resp <- "an"
# nms <- colnames(df_AN)
# nms <- nms[nms != resp]
# rhs <- paste('s(', nms, ')', sep = '', collapse = ' + ')
# fml <- paste(resp, '~', rhs, collapse = ' ')
# fml <- as.formula(fml)
# m1_an <- mgcv::gam(formula = fml, data = df_AN, family = binomial)

#check model assumptions
pdf("figures/anoxic_logreg_model_assumptions.pdf", height = 10, width = 10)
performance::check_model(m_an, check = "all", residual_type = "normal")
dev.off()
pdf("figures/oxic_logreg_model_assumptions.pdf", height = 10, width = 10)
performance::check_model(m_ox, check = "all", residual_type = "normal")
dev.off()
pdf("figures/lowfreq_logreg_model_assumptions.pdf", height = 10, width = 10)
performance::check_model(m_lf, check = "all", residual_type = "normal")
dev.off()
pdf("figures/highfreq_logreg_model_assumptions.pdf", height = 10, width = 10)
performance::check_model(m_hf, check = "all", residual_type = "normal")
dev.off()
#extract predictions
ggeffects::ggeffect(m_an) |>
  plot() |>
  sjPlot::plot_grid()

ggeffects::ggeffect(m_ox) |>
  plot() |>
  sjPlot::plot_grid()

ggeffects::ggeffect(m_lf) |>
  plot() |>
  sjPlot::plot_grid()

ggeffects::ggeffect(m_hf) |>
  plot() |>
  sjPlot::plot_grid()


# display contrasts / post-hocs /pairwise comparisons

#DFE genes encode multiheme cytochromes
an_table <- gtsummary::tbl_regression(
  m_an,
  exponentiate = TRUE,
  add_pairwise_contrasts = TRUE,
  contrasts_adjust = "bonferroni",
  pairwise_reverse = TRUE,
  pvalue_fun = ~gtsummary::style_pvalue(.x, digits = 3),
  label = list(other_reductases_mercury_reduction = "Mercury reduction (EC:1.16.1.1)",
               etc = "Full Electron Transport Chain (Complex I, II, III, IV, and V)",
               butyrate = "Butyrate utilization (K00634, K00929, K01896)",
               propionate = "Propionate utilization (K19697, K01026)",
               acetate = "Acetate utilization (K00625, K00925, K01905, K01067)",
               lactate = "Lactate utilization (K00016, K00101, K03778, K03777)",
               acetyl_co_a = "Pyruvate => Acetyl CoA (K00163, K00174, K00656)",
               arsenate_reduction = "Arsenate reduction (EC:1.20.4.1, EC:1.20.99.1)",
               sulfur = "Sulfur metabolism (K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357)",
               nitrogen = "Nitrogen metabolism (EC:1.14.99.39, K10944, K10945, K10946, K10535, K20932, K20933, K20934, K20935, K00370, K00371, K02588, K02586, K22896)",
               iron_oxidation = "Iron oxidation (Cyc1, Cyc2, MtoA, MtrB, sulfocyanin, CymA, FoxY, FoxE)",
               iron_reduction = "Iron reduction (MtrA, MtrB, MtrC, MtoA, OmcF, OmcS, DFE_0462, DFE_0463, DFE_0461, DFE_0464, DFE_0465)",
               siderophore = "Iron aquisition (Siderophore)")
) |>
  gtsummary::add_significance_stars(hide_p = F, hide_se = T, hide_ci = F) |>
  gtsummary::bold_p()

fhf |> #count(subcategory)
  filter(str_detect(subcategory, "Sulfur")) |> #count(function_name)
  filter(str_detect(function_name, "arsenate")) 


fhf |> #count(subcategory)
  filter(str_detect(subcategory, "Nitrogen") | str_detect(subcategory, "Nitrification") | str_detect(subcategory, "denitrification")) |>
  separate_rows(function_ids, sep = ", ") |>
  pull(function_ids)

cytochrome_genes <- gsf |> 
  filter(sheet == "Energy") |>
  filter(str_detect(module, "Cytochrome c oxidase") |
           str_detect(module, "Cytochrome o ubiquinol") |
           str_detect(module, "Cytochrome aa3") |
           str_detect(module, "Cytochrome bc1") |
           str_detect(module, "Cytochrome bd ubiquinol") |
           str_detect(module, "Fumarate reductase") |
           str_detect(module, "Succinate dehydrogenase") |
           str_detect(module, "NADH dehydrogenase") |
           str_detect(module, "NAD(P)") |
           str_detect(module, "NADH:quinone") |
           str_detect(module, "F-type ATPase") |
           str_detect(module, "V/A-type ATPase")) |>
  mutate(I = ifelse(str_detect(module, "NADH dehydrogenase") |
                      str_detect(module, "NAD(P)") |
                      str_detect(module, "NADH:quinone"), T, F),
         II = ifelse(str_detect(module, "Fumarate reductase") |
                       str_detect(module, "Succinate dehydrogenase"), T, F),
         III = ifelse(str_detect(module, "Cytochrome bd ubiquinol") |
                        str_detect(module, "Cytochrome bc1"), T, F),
         IV = ifelse(str_detect(module, "Cytochrome c oxidase") |
                       str_detect(module, "Cytochrome bd ubiquinol") |
                     str_detect(module, "Cytochrome o ubiquinol") |
                     str_detect(module, "Cytochrome aa3"), T, F),
         V = ifelse(str_detect(module, "F-type ATPase") |
                      str_detect(module, "V/A-type ATPase"), T, F)) |>
  pivot_longer(cols = c(I, II, III, IV, V),
               names_to = "complex",
               values_to = "pa") |>
  filter(pa) |>
  group_by(complex) |>
  summarise(genes = paste(gene_id, collapse = ", "))

FeGenie_geneSummary |> 
  filter(str_detect(category, "siderophore")) |>
  pull(HMM) |>
  unique()

ox_table <- gtsummary::tbl_regression(
  m_ox,
  exponentiate = TRUE,
  add_pairwise_contrasts = TRUE,
  contrasts_adjust = "bonferroni",
  pairwise_reverse = TRUE,
  pvalue_fun = ~gtsummary::style_pvalue(.x, digits = 3),
  label = list(other_reductases_mercury_reduction = "Mercury reduction (EC:1.16.1.1)",
               etc = "Full Electron Transport Chain (Complex I, II, III, IV, and V)",
               butyrate = "Butyrate utilization (K00634, K00929, K01896)",
               propionate = "Propionate utilization (K19697, K01026)",
               acetate = "Acetate utilization (K00625, K00925, K01905, K01067)",
               lactate = "Lactate utilization (K00016, K00101, K03778, K03777)",
               acetyl_co_a = "Pyruvate => Acetyl CoA (K00163, K00174, K00656)",
               arsenate_reduction = "Arsenate reduction (EC:1.20.4.1, EC:1.20.99.1)",
               sulfur = "Sulfur metabolism (K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357)",
               nitrogen = "Nitrogen metabolism (EC:1.14.99.39, K10944, K10945, K10946, K10535, K20932, K20933, K20934, K20935, K00370, K00371, K02588, K02586, K22896)",
               iron_oxidation = "Iron oxidation (Cyc1, Cyc2, MtoA, MtrB, sulfocyanin, CymA, FoxY, FoxE)",
               iron_reduction = "Iron reduction (MtrA, MtrB, MtrC, MtoA, OmcF, OmcS, DFE_0462, DFE_0463, DFE_0461, DFE_0464, DFE_0465)",
               siderophore = "Iron aquisition (Siderophore)")
) |>
  gtsummary::add_significance_stars(hide_p = F, hide_se = T, hide_ci = F) |>
  gtsummary::bold_p()

lf_table <- gtsummary::tbl_regression(
  m_lf,
  exponentiate = TRUE,
  add_pairwise_contrasts = TRUE,
  contrasts_adjust = "bonferroni",
  pairwise_reverse = TRUE,
  pvalue_fun = ~gtsummary::style_pvalue(.x, digits = 3),
  label = list(other_reductases_mercury_reduction = "Mercury reduction (EC:1.16.1.1)",
               etc = "Full Electron Transport Chain (Complex I, II, III, IV, and V)",
               butyrate = "Butyrate utilization (K00634, K00929, K01896)",
               propionate = "Propionate utilization (K19697, K01026)",
               acetate = "Acetate utilization (K00625, K00925, K01905, K01067)",
               lactate = "Lactate utilization (K00016, K00101, K03778, K03777)",
               acetyl_co_a = "Pyruvate => Acetyl CoA (K00163, K00174, K00656)",
               arsenate_reduction = "Arsenate reduction (EC:1.20.4.1, EC:1.20.99.1)",
               sulfur = "Sulfur metabolism (K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357)",
               nitrogen = "Nitrogen metabolism (EC:1.14.99.39, K10944, K10945, K10946, K10535, K20932, K20933, K20934, K20935, K00370, K00371, K02588, K02586, K22896)",
               iron_oxidation = "Iron oxidation (Cyc1, Cyc2, MtoA, MtrB, sulfocyanin, CymA, FoxY, FoxE)",
               iron_reduction = "Iron reduction (MtrA, MtrB, MtrC, MtoA, OmcF, OmcS, DFE_0462, DFE_0463, DFE_0461, DFE_0464, DFE_0465)",
               siderophore = "Iron aquisition (Siderophore)")
) |>
  gtsummary::add_significance_stars(hide_p = F, hide_se = T, hide_ci = F) |>
  gtsummary::bold_p()

hf_table <- gtsummary::tbl_regression(
  m_hf,
  exponentiate = TRUE,
  add_pairwise_contrasts = TRUE,
  contrasts_adjust = "bonferroni",
  pairwise_reverse = TRUE,
  pvalue_fun = ~gtsummary::style_pvalue(.x, digits = 3),
  label = list(other_reductases_mercury_reduction = "Mercury reduction (EC:1.16.1.1)",
               etc = "Full Electron Transport Chain (Complex I, II, III, IV, and V)",
               butyrate = "Butyrate utilization (K00634, K00929, K01896)",
               propionate = "Propionate utilization (K19697, K01026)",
               acetate = "Acetate utilization (K00625, K00925, K01905, K01067)",
               lactate = "Lactate utilization (K00016, K00101, K03778, K03777)",
               acetyl_co_a = "Pyruvate => Acetyl CoA (K00163, K00174, K00656)",
               arsenate_reduction = "Arsenate reduction (EC:1.20.4.1, EC:1.20.99.1)",
               sulfur = "Sulfur metabolism (K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357)",
               nitrogen = "Nitrogen metabolism (EC:1.14.99.39, K10944, K10945, K10946, K10535, K20932, K20933, K20934, K20935, K00370, K00371, K02588, K02586, K22896)",
               iron_oxidation = "Iron oxidation (Cyc1, Cyc2, MtoA, MtrB, sulfocyanin, CymA, FoxY, FoxE)",
               iron_reduction = "Iron reduction (MtrA, MtrB, MtrC, MtoA, OmcF, OmcS, DFE_0462, DFE_0463, DFE_0461, DFE_0464, DFE_0465)",
               siderophore = "Iron aquisition (Siderophore)")
) |>
  gtsummary::add_significance_stars(hide_p = F, hide_se = T, hide_ci = F) |>
  gtsummary::bold_p()


an_table |>
  gtsummary::as_flex_table() |>
  flextable::save_as_image(path = "figures/anoxic_logistic_regression_table.svg", res = 999)
#sig variables = ETC, propionate, acetate, lactate, acetyl_co_a, siderophore

ox_table |>
  gtsummary::as_flex_table() |>
  flextable::save_as_image(path = "figures/oxic_logistic_regression_table.svg", res = 999)
#sig variables = arsenate_reduction, sulfur

lf_table |>
  gtsummary::as_flex_table() |>
  flextable::save_as_image(path = "figures/lf_logistic_regression_table.svg", res = 999)
#sig variables = sulfur

hf_table |>
  gtsummary::as_flex_table() |>
  flextable::save_as_image(path = "figures/hf_logistic_regression_table.svg", res = 999)
#sig variables = ETC, butyrate, arsenate_reduction, sulfur

#which predictors are most important
anova_an <- car::Anova(m_an) |>
  rownames_to_column(var = "predictor")
anova_ox <- car::Anova(m_ox) |>
  rownames_to_column(var = "predictor")
anova_lf <- car::Anova(m_lf) |>
  rownames_to_column(var = "predictor")
anova_hf <- car::Anova(m_hf) |>
  rownames_to_column(var = "predictor")

anova_df <- bind_rows(anova_an, anova_ox, anova_lf, anova_hf, .id = "redox") |>
  mutate(redox = case_when(redox == 1 ~ "Anoxic",
                           redox == 2 ~ "Oxic", 
                           redox == 3 ~ "Low Freq.",
                           redox == 4 ~ "High Freq."))

predictor_levels <- anova_df |>
  filter(`Pr(>Chisq)` < 0.05) |>
  arrange((`LR Chisq`)) |>
  pull(predictor) |>
  unique()

anova_plot_df <- anova_df |>
  filter(`Pr(>Chisq)` < 0.05) |>
  mutate(redox = factor(redox, levels = c("Anoxic", "Low Freq.", "High Freq.", "Oxic")),
         predictor = factor(predictor, levels = predictor_levels))

short_labels <- c(
  "acetate" = "Acetate utilization",
  "butyrate" = "Butyrate utilization",
  "acetyl_co_a" = "Pyruvate => Acetyl CoA",
  "propionate" = "Propionate utilization",
  "sulfur" = "Sulfur metabolism",
  "etc" = "Full Electron Transport Chain (Complex I, II, III, IV, and V)",
  "siderophore" = "Iron aquisition (Siderophore)",
  "arsenate_reduction" = "Arsenate reduction",
  "lactate" = "Lactate utilization"
)

pdf("figures/logistic_regression_predictor_anova.pdf", height = 15, width = 10)
anova_plot_df |>
  ggplot(aes(x = `LR Chisq`, y = predictor)) +
  geom_bar(stat = "identity", 
           position = "dodge", 
           aes(group = redox, fill = redox),
           alpha = 0.5,
           width = 0.7) +
  scale_fill_manual(values = c(Anoxic, LF, HF, Oxic)) +
  scale_y_discrete(labels = short_labels) +
  labs(x = "Log Ratio ChiSq.",
       y = NULL,
       fill = NULL) +
  theme_bw()
dev.off()


#model performance evaluation
performance::performance(m_an) 
performance::performance(m_ox)
performance::performance(m_lf)
performance::performance(m_hf)

performance_metrics <- list(
  anoxic = performance::performance(m_an),
  oxic = performance::performance(m_ox),
  low_flux = performance::performance(m_lf),
  high_flux = performance::performance(m_hf)
)

metrics_df <- bind_rows(
  lapply(performance_metrics, as.data.frame),
  .id = "Model"
)

metrics_table <- metrics_df %>%
  flextable::flextable() %>%
  flextable::autofit() %>%
  flextable::set_caption(caption = "Model Performance Metrics")

flextable::save_as_image(metrics_table, path = "figures/model_performance_metrics.svg", res = 999)

effectsize::interpret_r2(0.169)
effectsize::interpret_r2(0.039)
effectsize::interpret_r2(0.032)
effectsize::interpret_r2(0.082)


roc_an <- pROC::roc(an ~ fitted.values(m_an), data = df_AN, plot = TRUE, 
          legacy.axes = TRUE, print.auc = TRUE, ci = TRUE)

roc_ox <- pROC::roc(ox ~ fitted.values(m_ox), data = df_OX, plot = TRUE, 
          legacy.axes = TRUE, print.auc = TRUE, ci = TRUE)

roc_lf <- pROC::roc(lf ~ fitted.values(m_lf), data = df_LF, plot = TRUE, 
          legacy.axes = TRUE, print.auc = TRUE, ci = TRUE)

roc_hf <- pROC::roc(hf ~ fitted.values(m_hf), data = df_HF, plot = TRUE, 
          legacy.axes = TRUE, print.auc = TRUE, ci = TRUE)

auc_an <- round(pROC::auc(roc_an), 2)
auc_ox <- round(pROC::auc(roc_ox), 2)
auc_lf <- round(pROC::auc(roc_lf), 2)
auc_hf <- round(pROC::auc(roc_hf), 2)

pdf("figures/logistic_regression_ROCAUC.pdf", height = 10, width = 10)
plot(roc_an, col = Anoxic, legacy.axes = TRUE, print.auc = FALSE, 
     main = "ROC Curves", xlab = "1 - Specificity", ylab = "Sensitivity")
lines(roc_lf, col = LF)
lines(roc_hf, col = HF)
lines(roc_ox, col = Oxic)
legend("bottomright", 
       legend = c(paste("Anoxic (AUC =", auc_an, ")"),
                  paste("Low Freq. (AUC =", auc_lf, ")"),
                  paste("High Freq. (AUC =", auc_hf, ")"),
                  paste("Oxic (AUC =", auc_ox, ")")), 
       col = c(Anoxic, LF, HF, Oxic), 
       lwd = 2)
dev.off()

# AUC values
# 0.9 - 1.0 Excelent
# 0.8 - 0.9 Very good
# 0.7 - 0.8 Good
# 0.6 - 0.7 Satisfactory
# 0.5 - 0.6 Unsatisfactory



an_tot <- df |>
  mutate_if(is.factor, as.logical) |>
  filter(AN) |> 
  pivot_longer(cols = everything(),
               names_to = "var",
               values_to = "value") |>
  mutate(value = ifelse(var %in% c("iron_oxidation",
                                   "iron_reduction",
                                   "siderophore"), 
                        value/sum(as.logical(df$AN)),
                        value)) |>
  group_by(var) |>
  summarise(total = sum(value), .groups = "drop")

ox_tot <- df |>
  mutate_if(is.factor, as.logical) |>
  filter(OX) |> 
  pivot_longer(cols = everything(),
               names_to = "var",
               values_to = "value") |>
  mutate(value = ifelse(var %in% c("iron_oxidation",
                                   "iron_reduction",
                                   "siderophore"), 
                        value/sum(as.logical(df$OX)),
                        value)) |> 
  group_by(var) |>
  summarise(total = sum(value))

lf_tot <- df |>
  mutate_if(is.factor, as.logical) |>
  filter(LF) |> 
  pivot_longer(cols = everything(),
               names_to = "var",
               values_to = "value") |>
  mutate(value = ifelse(var %in% c("iron_oxidation",
                                   "iron_reduction",
                                   "siderophore"), 
                        value/sum(as.logical(df$LF)),
                        value)) |> 
  group_by(var) |>
  summarise(total = sum(value))

hf_tot <- df |>
  mutate_if(is.factor, as.logical) |>
  filter(HF) |> 
  pivot_longer(cols = everything(),
               names_to = "var",
               values_to = "value") |>
  mutate(value = ifelse(var %in% c("iron_oxidation",
                                   "iron_reduction",
                                   "siderophore"), 
                        value/sum(as.logical(df$HF)),
                        value)) |> 
  group_by(var) |>
  summarise(total = sum(value))

df_tot <- bind_rows(an_tot, ox_tot, lf_tot, hf_tot, .id = "redox") |>
  mutate(redox = case_when(redox == 1 ~ "Anoxic",
                           redox == 2 ~ "Oxic", 
                           redox == 3 ~ "LF",
                           redox == 4 ~ "HF")) 

df_tot_p <- df_tot |>
  mutate(facet = case_when(str_detect(var, "iron") | str_detect(var, "siderophore") ~ "Mean Gene Count",
                           var %in% c("AN", "OX", "LF", 'HF') ~ "Active Redox Conditions",
                           T ~ "Functional Potential"),
         redox = case_when(redox == "Anoxic" ~ "Anoxic",
                           redox == "Oxic" ~ "Oxic",
                           redox == "LF" ~ "Low Freq.",
                           redox == "HF" ~ "High Freq."),
         redox = factor(redox, levels = c("Anoxic", "Low Freq.", "High Freq.", "Oxic")),
         var = case_when(
           var == "Other Reductases: mercury reduction" ~ "Mercury reduction (EC:1.16.1.1)",
           var == "ETC" ~ "Full Electron Transport Chain (Complex I, II, III, IV, and V)",
           var == "butyrate" ~ "Butyrate utilization (K00634, K00929, K01896)",
           var == "propionate" ~ "Propionate utilization (K19697, K01026)",
           var == "acetate" ~ "Acetate utilization (K00625, K00925, K01905, K01067)",
           var == "lactate" ~ "Lactate utilization (K00016, K00101, K03778, K03777)",
           var == "acetylCoA" ~ "Pyruvate => Acetyl CoA (K00163, K00174, K00656)",
           var == "arsenate_reduction" ~ "Arsenate reduction (EC:1.20.4.1, EC:1.20.99.1)",
           var == "sulfur" ~ "Sulfur metabolism (K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357)",
           var == "nitrogen" ~ "Nitrogen metabolism (EC:1.14.99.39, K10944, K10945, K10946, K10535, K20932, K20933, K20934, K20935, K00370, K00371, K02588, K02586, K22896)",
           var == "iron_oxidation" ~ "Iron oxidation (Cyc1, Cyc2, MtoA, MtrB, sulfocyanin, CymA, FoxY, FoxE)",
           var == "iron_reduction" ~ "Iron reduction (MtrA, MtrB, MtrC, MtoA, OmcF, OmcS, DFE_0462, DFE_0463, DFE_0461, DFE_0464, DFE_0465)",
           var == "siderophore" ~ "Iron aquisition (Siderophore)",
           var == "AN" ~ "Anoxic",
           var == "OX" ~ "Oxic",
           var == "LF" ~ "Low Freq.",
           var == "HF" ~ "High Freq.",
           TRUE ~ var 
         ),
         var = factor(var, levels = rev(c(
           "Anoxic", 
           "Low Freq.", 
           "High Freq.", 
           "Oxic", 
           "Full Electron Transport Chain (Complex I, II, III, IV, and V)", 
           "Mercury reduction (EC:1.16.1.1)", 
           "Arsenate reduction (EC:1.20.4.1, EC:1.20.99.1)", 
           "Sulfur metabolism (K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357)", 
           "Nitrogen metabolism (EC:1.14.99.39, K10944, K10945, K10946, K10535, K20932, K20933, K20934, K20935, K00370, K00371, K02588, K02586, K22896)", 
           "Butyrate utilization (K00634, K00929, K01896)", 
           "Propionate utilization (K19697, K01026)", 
           "Acetate utilization (K00625, K00925, K01905, K01067)", 
           "Lactate utilization (K00016, K00101, K03778, K03777)", 
           "Pyruvate => Acetyl CoA (K00163, K00174, K00656)", 
           "Iron aquisition (Siderophore)", 
           "Iron reduction (MtrA, MtrB, MtrC, MtoA, OmcF, OmcS, DFE_0462, DFE_0463, DFE_0461, DFE_0464, DFE_0465)", 
           "Iron oxidation (Cyc1, Cyc2, MtoA, MtrB, sulfocyanin, CymA, FoxY, FoxE)"
         )))) 


df_tot_p_subset <- df_tot_p |>
  mutate(keep = case_when(
    redox == "Anoxic" & var %in% c(
      "Full Electron Transport Chain (Complex I, II, III, IV, and V)", 
      "Propionate utilization (K19697, K01026)", 
      "Acetate utilization (K00625, K00925, K01905, K01067)", 
      "Lactate utilization (K00016, K00101, K03778, K03777)", 
      "Pyruvate => Acetyl CoA (K00163, K00174, K00656)", 
      "Iron aquisition (Siderophore)"
    ) ~ TRUE,
    redox == "Oxic" & var %in% c(
      "Arsenate reduction (EC:1.20.4.1, EC:1.20.99.1)", 
      "Sulfur metabolism (K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357)"
    ) ~ TRUE,
    redox == "Low Freq." & var %in% c(
      "Sulfur metabolism (K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357)"
    ) ~ TRUE,
    redox == "High Freq." & var %in% c(
      "Full Electron Transport Chain (Complex I, II, III, IV, and V)", 
      "Butyrate utilization (K00634, K00929, K01896)", 
      "Arsenate reduction (EC:1.20.4.1, EC:1.20.99.1)", 
      "Sulfur metabolism (K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357)"
    ) ~ TRUE,
    TRUE ~ FALSE
  )) |>
  filter(keep)

df_tot_p_subset2 <- df_tot_p |>
  mutate(keep = case_when(
    redox == "Anoxic" & var == "Anoxic" ~ TRUE,
    redox == "Oxic" & var == "Oxic" ~ TRUE,
    redox == "Low Freq." & var == "Low Freq." ~ TRUE,
    redox == "High Freq." & var == "High Freq." ~ TRUE,
    TRUE ~ FALSE
  )) |>
  filter(keep)

pdf("figures/logistic_regression_sig_table.pdf", height = 10, width = 10)
df_tot_p |>
  ggplot(aes(x = redox, y = var)) +
  geom_tile(fill = "white", color = "black") +
  geom_tile(data = df_tot_p_subset2, aes(fill = redox), alpha = 0.5, color = "black") +
  geom_tile(data = df_tot_p_subset, fill = "transparent", color = "black", linewidth = 1) +
  geom_text(aes(label = round(total, 3))) +
  scale_fill_manual(values = c(Anoxic, LF, HF, Oxic)) +
  facet_wrap(~facet, nrow = 3, scales = "free") +
  theme_minimal()
dev.off()

df_tot_p_subset <- df_tot_p |>
  mutate(keep = case_when(
    redox == "Anoxic" & var %in% c(
      "Full Electron Transport Chain (Complex I, II, III, IV, and V)", 
      "Propionate utilization (K19697, K01026)", 
      "Acetate utilization (K00625, K00925, K01905, K01067)", 
      "Lactate utilization (K00016, K00101, K03778, K03777)", 
      "Pyruvate => Acetyl CoA (K00163, K00174, K00656)", 
      "Iron aquisition (Siderophore)"
    ) ~ TRUE,
    redox == "Oxic" & var %in% c(
      "Arsenate reduction (EC:1.20.4.1, EC:1.20.99.1)", 
      "Sulfur metabolism (K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357)"
    ) ~ TRUE,
    redox == "Low Freq." & var %in% c(
      "Sulfur metabolism (K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357)"
    ) ~ TRUE,
    redox == "High Freq." & var %in% c(
      "Full Electron Transport Chain (Complex I, II, III, IV, and V)", 
      "Butyrate utilization (K00634, K00929, K01896)", 
      "Arsenate reduction (EC:1.20.4.1, EC:1.20.99.1)", 
      "Sulfur metabolism (K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357)"
    ) ~ TRUE,
    TRUE ~ FALSE
  )) |>
  filter(keep)|>
  mutate(var = str_replace_all(
    var, 
    "\\((?!(Siderophore|Complex I, II, III, IV, and V))[^)]*\\)", 
    ""
  ),
  var = str_trim(var)) |>
  mutate(var = factor(
    var, 
    levels = rev(c(
      "Anoxic", 
      "Low Freq.", 
      "High Freq.", 
      "Oxic", 
      "Full Electron Transport Chain (Complex I, II, III, IV, and V)", 
      "Mercury reduction", 
      "Arsenate reduction", 
      "Sulfur metabolism", 
      "Nitrogen metabolism", 
      "Butyrate utilization", 
      "Propionate utilization", 
      "Acetate utilization", 
      "Lactate utilization", 
      "Pyruvate => Acetyl CoA", 
      "Iron aquisition (Siderophore)", 
      "Iron reduction", 
      "Iron oxidation"
    ))
  ))

df_tot_p_subset2 <- df_tot_p |>
  mutate(keep = case_when(
    redox == "Anoxic" & var == "Anoxic" ~ TRUE,
    redox == "Oxic" & var == "Oxic" ~ TRUE,
    redox == "Low Freq." & var == "Low Freq." ~ TRUE,
    redox == "High Freq." & var == "High Freq." ~ TRUE,
    TRUE ~ FALSE
  )) |>
  filter(keep)|>
  mutate(var = str_replace_all(
    var, 
    "\\([^()]*\\)(?!\\s*Siderophore\\)|\\s*Complex I, II, III, IV, and V\\))", 
    ""
  ),
  var = str_trim(var)) |>
  mutate(var = factor(
    var, 
    levels = rev(c(
      "Anoxic", 
      "Low Freq.", 
      "High Freq.", 
      "Oxic", 
      "Full Electron Transport Chain (Complex I, II, III, IV, and V)", 
      "Mercury reduction", 
      "Arsenate reduction", 
      "Sulfur metabolism", 
      "Nitrogen metabolism", 
      "Butyrate utilization", 
      "Propionate utilization", 
      "Acetate utilization", 
      "Lactate utilization", 
      "Pyruvate => Acetyl CoA", 
      "Iron aquisition (Siderophore)", 
      "Iron reduction", 
      "Iron oxidation"
    ))
  ))


pdf("figures/logistic_regression_sig_table_no_genes.pdf", height = 10, width = 10)
df_tot_p |>
  mutate(var = str_replace_all(
    var, 
    "\\((?!(Siderophore|Complex I, II, III, IV, and V))[^)]*\\)", 
    ""
  ),
  var = str_trim(var)) |> 
  mutate(var = factor(
    var, 
    levels = rev(c(
      "Anoxic", 
      "Low Freq.", 
      "High Freq.", 
      "Oxic", 
      "Full Electron Transport Chain (Complex I, II, III, IV, and V)", 
      "Mercury reduction", 
      "Arsenate reduction", 
      "Sulfur metabolism", 
      "Nitrogen metabolism", 
      "Butyrate utilization", 
      "Propionate utilization", 
      "Acetate utilization", 
      "Lactate utilization", 
      "Pyruvate => Acetyl CoA", 
      "Iron aquisition (Siderophore)", 
      "Iron reduction", 
      "Iron oxidation"
    ))
  )) |>
  ggplot(aes(x = redox, y = var)) +
  geom_tile(fill = "white", color = "black") +
  geom_tile(data = df_tot_p_subset2, aes(fill = redox), alpha = 0.5, color = "black") +
  geom_tile(data = df_tot_p_subset, fill = "transparent", color = "black", linewidth = 1) +
  geom_text(aes(label = round(total, 3))) +
  scale_fill_manual(values = c(Anoxic, LF, HF, Oxic)) +
  facet_wrap(~facet, nrow = 3, scales = "free") +
  theme_minimal()
dev.off()

#genes table for logistic regression variables
# Create the dataframe with additional columns I, II, III, IV, and V
df_table_out <- data.frame(
  `function` = c(
    "Mercury reduction (EC:1.16.1.1)", 
    "Full Electron Transport Chain (Complex I, II, III, IV, and V)", 
    "Butyrate utilization (K00634, K00929, K01896)", 
    "Propionate utilization (K19697, K01026)", 
    "Acetate utilization (K00625, K00925, K01905, K01067)", 
    "Lactate utilization (K00016, K00101, K03778, K03777)", 
    "Pyruvate => Acetyl CoA (K00163, K00174, K00656)", 
    "Arsenate reduction (EC:1.20.4.1, EC:1.20.99.1)", 
    "Sulfur metabolism (K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357)", 
    "Nitrogen metabolism (EC:1.14.99.39, K10944, K10945, K10946, K10535, K20932, K20933, K20934, K20935, K00370, K00371, K02588, K02586, K22896)", 
    "Iron oxidation (Cyc1, Cyc2, MtoA, MtrB, sulfocyanin, CymA, FoxY, FoxE)", 
    "Iron reduction (MtrA, MtrB, MtrC, MtoA, OmcF, OmcS, DFE_0462, DFE_0463, DFE_0461, DFE_0464, DFE_0465)", 
    "Iron aquisition (Siderophore)"
  ),
  genes = c(
    "EC:1.16.1.1", 
    "Complex I, II, III, IV, and V", 
    "K00634, K00929, K01896", 
    "K19697, K01026", 
    "K00625, K00925, K01905, K01067", 
    "K00016, K00101, K03778, K03777", 
    "K00163, K00174, K00656", 
    "EC:1.20.4.1, EC:1.20.99.1", 
    "K11180, K11181, K08352, K17222, K17224, K17225, K22622, K17226, K17227, K08357", 
    "EC:1.14.99.39, K10944, K10945, K10946, K10535, K20932, K20933, K20934, K20935, K00370, K00371, K02588, K02586, K22896", 
    "Cyc1, Cyc2, MtoA, MtrB, sulfocyanin, CymA, FoxY, FoxE", 
    "MtrA, MtrB, MtrC, MtoA, OmcF, OmcS, DFE_0462, DFE_0463, DFE_0461, DFE_0464, DFE_0465", 
    "FpvE-family-permease, FpvD-family-siderophore-transport, PirA-family-siderophore-receptor, TonB-family, ExbD-family, ExbB-family, FpvC-family-siderophore-transport, PvdT-family-siderophore-export, PvdR-family-sideropore-export, HatD-family-substrate-binding-protein, VabG-family-siderophore-synthesis, PchI-family-siderophore-synthesis, PchH-family-siderophore-synthesis, EntS-family-siderophore-export, LbtU-family-siderophore receptor, PiuA-family-siderophore-receptor, FpvG-family-siderophore-transport, PvsE-family-siderophore-synthesis, PvuD-family-ATP-binding-protein, VabE-family-siderophore-synthesis, RhbA-family-siderophore-synthesis, RhbB-family-siderophore-synthesis, FptX-family-siderophore-transport, PvdM-family-siderophore-synthesis, VabF-family-siderophore-synthesis, PvsC-family-siderophore-synthesis, PvdN-family-siderophore-synthesis, PvdE-family-siderophore-synthesis, FpvF-family-siderophore-transport, FptC-family-siderophore-transport, IucD-family-siderophore-synthesis, RhbD-family-siderophore-synthesis, PvsD-family-siderophore-synthesis, PvuCD-family-permease, YqjH-family-iron-reductase, ViuB-family-siderophore-utilization, PvuE-family-ATP-binding-protein, PchC-family-siderophore-synthesis, IroC-family-ATP-binding-protein"
  ),
  source = c(
    "DRAM", 
    "DRAM", 
    "DRAM", 
    "DRAM", 
    "DRAM", 
    "DRAM", 
    "DRAM", 
    "DRAM", 
    "DRAM", 
    "DRAM", 
    "FeGenie", 
    "FeGenie", 
    "FeGenie"
  ),
  I = c(
    NA, 
    "K00330, K00331, K00332, K00333, K13378, K13380, K00334, K00335, K00336, K00337, K00338, K00339, K00340, K00341, K00342, K15863, K00343, K03934, K03935, K03936, K03937, K03938, K03939, K03940, K03941, K03942, K03943, K03944, K03945, K03946, K03947, K03948, K03949, K03950, K03951, K03952, K03953, K03954, K03955, K03956, K11352, K11353, K03957, K03958, K03959, K03960, K03961, K03962, K03963, K03964, K03965, K03966, K11351, K03967, K03968", 
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
  ),
  II = c(
    NA, 
    "K00236, K00237, K00234, K00235, K00241, K00242, K00239, K00240, K00244, K00245, K00246, K00247", 
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
  ),
  III = c(
    NA, 
    "K00412, K00413, K00410, K00411, K03886, K03887, K03888, K03890, K03891, K03889, K00412, K00413, K00410, K00411, K00414, K00415, K00416, K00417, K00418, K00419, K00420, K00425, K00426, K00424, K22501", 
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
  ),
  IV = c(
    NA, 
    "K00425, K00426, K00424, K22501, K02257, K02262, K02256, K02261, K02263, K02264, K02265, K02266, K02267, K02268, K02269, K02270, K02271, K02272, K02273, K02258, K02259, K02260, K02275, K02274, K02276, K15408, K02277, K00404, K00405, K15862, K00407, K00406, K02827, K02826, K02828, K02829, K02297, K02298, K02299, K02300", 
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
  ),
  V = c(
    NA, 
    "K02111, K02112, K02113, K02114, K02115, K02108, K02109, K02110, K02117, K02118, K02119, K02120, K02121, K02122, K02107, K02123, K02124, K02132, K02133, K02136, K02134, K02135, K02137, K02126, K02127, K02128, K02138, K02129, K01549, K02130, K02139, K02140, K02141, K02131, K02142, K02143, K02125", 
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
  ),
  stringsAsFactors = FALSE
)


df_table_out

wb <- createWorkbook()
addWorksheet(wb, "log_reg_predictor_genes")
writeData(wb, "log_reg_predictor_genes", df_table_out)
saveWorkbook(wb, "/Users/leleiwi1/Desktop/LLNL_postdoc/Soil_Microbiome_SFA/Great_Redox_Experiment/MAG_Paper/Table3.xlsx", overwrite = TRUE)

#mags with function for text
df |>
  filter()