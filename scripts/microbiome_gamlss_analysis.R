# Microbiome Composition and Covariates

# --- Libraries ---
library(compositions)
library(broom)
library(gamlss)
library(dplyr)
library(gtsummary)
library(tidyr)
library(ggplot2)

# --- Functions ---
estimate_gamlss <- function(data, taxa, covariates) {
  results_list <- list()
  
  for (taxon in taxa) {
    res <- lapply(covariates, function(cov) {
      if (length(unique(data[[cov]])) < 2) return(NULL)
      
      formula <- as.formula(paste0("`", taxon, "` ~ `", cov, "`"))
      model <- tryCatch(gamlss(formula, family = BEZI, data = data), error = function(e) NULL)
      if (is.null(model)) return(NULL)
      
      vcov_mu <- tryCatch(vcov(model, what = "mu"), error = function(e) NULL)
      if (is.null(vcov_mu)) return(NULL)
      
      coefs <- coef(model, what = "mu")
      idx <- grep(cov, names(coefs), fixed = TRUE)
      
      if (length(idx) == 1) {
        est <- coefs[idx]
        se <- sqrt(diag(vcov_mu))[idx]
        data.frame(
          taxon = taxon,
          covariate = cov,
          estimate = est,
          std.error = se,
          statistic = est / se,
          p.value = 2 * (1 - pnorm(abs(est / se)))
        )
      } else {
        NULL
      }
    })
    
    res <- bind_rows(Filter(Negate(is.null), res))
    if (nrow(res) > 0) {
      res$p_adj <- p.adjust(res$p.value, method = "BH")
      results_list[[taxon]] <- res
    }
  }
  
  bind_rows(results_list)
}

plot_heatmap <- function(results, ylab) {
  results %>%
    mutate(
      label = ifelse(p_adj < 0.05, round(estimate, 4), ""),
      taxon = factor(taxon, levels = unique(taxon)),
      covariate = factor(covariate, levels = unique(covariate))
    ) %>%
    ggplot(aes(x = covariate, y = taxon, fill = estimate)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(x = "Covariate", y = ylab, fill = "Estimate (β)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
}

# --- Load preprocessed .RData files ---
load("path/to/class&covariates.RData")   # must contain object 'class'
load("path/to/family&covariates.RData")  # must contain object 'family'
load("path/to/genus&covariates.RData")   # must contain object 'genus'
load("path/to/phyla&covariates.RData")   # must contain object 'phyla'

# --- Covariates (shared across all levels) ---
covariates <- c("Female gender", "Age", "BMI", "Antibiotics",
                "Anti-inflammatories", "Proton pump inhibitors", "Probiotics/Prebiotics",
                "Disorders", "Hours of sleep", "Tobacco",
                "Pathologies", "Physical activity", "Allergies")

# --- CLASS analysis ---
class_taxa <- c("Clostridia %", "Bacteroidia %", "Actinobacteria %", "Negativicutes %",
                "Erysipelotrichia %", "Bacilli %", "Betaproteobacteria %", "Gammaproteobacteria %",
                "Deltaproteobacteria %", "Alphaproteobacteria %")
class_results <- estimate_gamlss(class, class_taxa, covariates)
plot_heatmap(class_results, ylab = "Class")

# --- FAMILY analysis ---
family_taxa <- c("Lachnospiraceae %", "Porphyromonadaceae %", "Ruminococcaceae %", "Bacteroidaceae %",
                 "Prevotellaceae %", "Coriobacteriaceae %", "Rikenellaceae %", "Eubacteriaceae %",
                 "Bifidobacteriaceae %", "Veillonellaceae %")
family_results <- estimate_gamlss(family, family_taxa, covariates)
plot_heatmap(family_results, ylab = "Family")

# --- GENUS analysis ---
genus_taxa <- c("Eubacterium %", "Collinsella %", "Bifidobacterium %", "Dorea %",
                "Clostridium IV %", "Prevotella %", "Barnesiella %", "Rummeliibacillus %",
                "Eisenbergiella %", "Oscillibacter %")
genus_results <- estimate_gamlss(genus, genus_taxa, covariates)
plot_heatmap(genus_results, ylab = "Genus")

# --- PHYLUM analysis ---
phyla_taxa <- c("Firmicutes %", "Bacteroidetes %", "Actinobacteria %",
                "Proteobacteria %", "Verrucomicrobia %", "Acidobacteria %",
                "Armatimonadetes %", "Tenericutes %", "Synergistetes %",
                "Cyanobacteria/Chloroplast %")
phyla_results <- estimate_gamlss(phyla, phyla_taxa, covariates)
plot_heatmap(phyla_results, ylab = "Phylum")
