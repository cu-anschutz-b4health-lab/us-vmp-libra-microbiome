#NFLPA VMP metagenomic analysis
setwd('/Users/mattigo/Documents/CU_Folders/Projects_Ongoing/NFLPA_Microbiome/Test_folder')
#Loads all data required to run this code
load("NFLPA_LIBRA_US_VMP_Microbiome.RData")
new_dir_path <- "Figures" 
set.seed(64)
# Check if the directory exists
if (!dir.exists(new_dir_path)) {
  # If it doesn't exist, create it
  dir.create(new_dir_path, recursive = TRUE) 
  message(paste("Directory created:", new_dir_path))
} else {
  message(paste("Directory already exists:", new_dir_path))
}
library(reshape2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ANCOMBC)
library(qiime2R)
library(phyloseq)
library(microbiome)
library(vegan)
library(knitr)
library(forcats)
library(tidyverse)
library(DT)
library(outliers)
library(beeswarm)
library(WGCNA)
library(flextable)
library(magrittr)
library(dplyr)
library(tidyheatmaps)
library(devtools)
library(ggvegan)
theme_set(theme_bw())
#####Demographics#####
#summarySE funciton 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#Data sets
######################################################
#after all these data were combined which was done in the **metadata-wrangling.R** 
# file, this was written out to check for consistency and make the Libra variable
libra.alpha.all
libra.alpha.all$LIBRA.score.quartile <- as.factor(libra.alpha.all$LIBRA.score.quartile)
#Examining the distribution of the libra variable
#this has already been added to the metadata and will be commented out 
lib.var <- libra.alpha.all$LIBRA.score
lib.var <- data.frame(lib.var)
myCol <- lapply(lib.var, function(x) cut(x, breaks = quantile(x, na.rm = T), labels = FALSE))
myCol
myCol$lib.var <- as.factor(myCol$lib.var)
#the equation could not comute 2 values that should be in Q1 FYI
#adding in quartile data and writing out for future use
#write.csv(libra.alpha.all, "LIBRA/LIBRA_variables_noNAs_021125-alpha1.csv", row.names = T)
#libra.alpha.all$LIBRA.score.quartile <- myCol$lib.var
#table(myCol)
########################################################
#Making some metadata summary plots
library(table1)
#p-value function
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  x <- x[-length(x)]  # Remove "overall" group
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    ano <- aov(y ~ g)
    p <- summary(ano)[[1]][[5]][1]
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=2, eps=0.001)))
  
}
p <- aov(Age ~ as.factor(LIBRA.score.quartile) , data = libra.alpha.all)
summary(p)
#Demographics
libra_table_2 = table1(~ Age + Race + Ethnicity + BMI + Homeless.ever + Homeless.Times + Multiple.Mild.TBIs, data = libra.alpha.all)
libra_table_3 = table1(~ LIBRA.score + Age + Race + Ethnicity + BMI + Homeless.ever + Homeless.Times + Multiple.Mild.TBIs | LIBRA.score.quartile, 
                       data = libra.alpha.all, overall = c(right="Total"),extra.col=list(`P-value`=pvalue))

#Libra vaiables
libra_table_2a = table1(~ LIBRA.score +  Obese +  Lifetime_smoking + Diabetes + Hypertension + 
                          Hypercholesterolemia + Coronary.heart.disease + 
                          Chronic_kidney_disease +  Current_alcohol + Depressed + Healthy_diet + Physical.activity.H.L, data = libra.alpha.all)


t1flex(libra_table_2) %>% 
  save_as_docx(path="libra_table_2.docx")

t1flex(libra_table_2a) %>% 
  save_as_docx(path="libra_table_2a.docx")


t1flex(libra_table_3) %>% 
  save_as_docx(path="libra_table_3.docx")

#### 
####Supplemental Figure 1
library(ggbeeswarm)
colvec <- c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A")

p_swarm <- ggplot(data = libra.alpha.all, aes(y = LIBRA.score, x="",color = LIBRA.score.quartile)) +
  geom_beeswarm(cex = 3)+
  scale_y_continuous(breaks = seq(0, 10, by=2))+
  scale_color_manual(values = colvec)+
  guides(color = guide_legend(reverse = TRUE))+
  labs(y="LIBRA Score", x="LIBRA Score Distribution", color="LIBRA Score Quartile")+
  theme_bw()

ggsave("Figures/S1_LIBRA_swarm_figure.jpeg", p_swarm, w=6, h=6, units="in")


####Figure 1 Alpha and Beta
##############Alpha#####
#Alpha plots
libra.alpha.all$LIBRA.score.quartile <- factor(libra.alpha.all$LIBRA.score.quartile)
melt.libra.alpha <- melt(libra.alpha.all, id.vars = c("SampleID", "LIBRA.score", "LIBRA.score.quartile"), 
                         variable.name = "alpha", measure.vars = c("observed_otus", "shannon"))
new_labels <- c("observed_otus" = "Observed OTUs", "shannon" = "Shannon Index")



p1a <- ggplot(data = melt.libra.alpha %>% dplyr::filter(alpha=="shannon"), aes(x = LIBRA.score, y = value)) +
  geom_point() + 
  geom_smooth(se = F, method = "lm") +
  stat_cor(label.y.npc="top", label.x = 5.25, size = 3)+
  theme_bw() + 
  labs(tag = "A", x="LIBRA Score", y="Shannon Index") +
  theme(plot.tag = element_text(size = 12, face = "bold"))


ggsave(paste0("Figures/Figure_1a_Alpha_Points", Sys.Date(),".jpeg"),p1a, height = 10, width = 14, units="in" )


p1b <- ggplot(data = melt.libra.alpha %>% dplyr::filter(alpha=="shannon"), aes(x = LIBRA.score.quartile, y = value)) +
  geom_boxplot(aes(fill = LIBRA.score.quartile)) + 
  geom_jitter(aes(color = LIBRA.score.quartile)) + 
  scale_color_manual(values=(colvec)) + 
  scale_fill_manual(values=(colvec))+
  theme_bw() + 
  labs(tag = "B", x="LIBRA Score Quartile", y="Shannon Index") +
  stat_compare_means(method = "kruskal.test", label.y.npc="top", label.x = 2.75, size = 3) + # Add Kruskal-Wallis p-value
  theme(plot.tag = element_text(size = 12, face = "bold"),
        legend.position = "none") 
ggsave(paste0("Figures/Figure_1b_Alpha_Boxplots", Sys.Date(),".jpeg"),p1b, height = 10, width = 14, units="in" )

summary(lm(shannon~LIBRA.score+ Age+ Homeless.ever, data = libra.alpha.all))



#####Beta######
x_beta
#next line will determine which metadata is being piped in so make sure to check before proceeding 
meta.in 
#matching up metadata with the distance matrix
ind <- match(x = rownames(x_beta), table = meta.in$SampleID)
ind1 <- !is.na(ind)
x1 <- x_beta[ind1,ind1]
ind2 <- match(x = rownames(x1), table = meta.in$SampleID)
meta.in1 <- meta.in[ind2,]
#this next line should be all true - if not, then something is wrong
table(rownames(x1) == meta.in1$SampleID) 
rownames(x1)
colnames(x1)
meta.in1$SampleID
#read in taxa
taxa_species 

Taxa <- t(taxa_species)
head(Taxa)

ind1 <- match(x = colnames(x1), table = rownames(Taxa))
Taxa1 <- Taxa[ind1,]
head(Taxa1)


#**Special plot** CCA start
#this code is for coloring dots on biplot, if not wanting that then next 3 lines can be removed
meta.in$LIBRA.score.quartile <- as.factor(meta.in$LIBRA.score.quartile)
with(meta.in, levels(LIBRA.score.quartile))
scl <- 2

#CCA start
pca_w=capscale(x1~1, data=meta.in1[,-1])
colvec <- c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A")
perm_result <- with(meta, adonis2(formula = x_beta~LIBRA.score.quartile+Age+BMI+Homeless.ever+Race+Ethnicity, data = meta, permutations = 999))



ggplot_beta = fortify(pca_w)
ggplot_beta = ggplot_beta %>% rename(SampleID = label)
beta_data_set = merge(ggplot_beta, meta.in, by="SampleID")
p1c = ggplot(beta_data_set, aes(x=MDS1, y=MDS2, color = factor(LIBRA.score.quartile, levels = c("4", "3", "2", "1")))) +
  geom_point()+
  theme_bw()+
  coord_cartesian()+
  scale_color_manual(values=rev(colvec))+
  stat_ellipse(type = "norm", linetype = 2) +
  labs(color = "Libra Quartile", tag = "C") +
  theme(plot.tag = element_text(size = 12, face = "bold")) +
  annotate("text",label=paste0("PERMANOVA = ",round(perm_result$`Pr(>F)`[1],2)), x=1.55,y=1.9, size =3.5)


ggsave(paste0("Figures/Figure_1c_Beta", Sys.Date(),".jpeg"),p1c, height = 10, width = 14, units="in" )

gg_top = ggarrange(p1a, p1b, ncol=2)
p1 = ggarrange(gg_top, p1c, nrow=2)

ggsave(paste0("Figures/Figure_1_Alpha_Beta_Diversity", Sys.Date(),".jpeg"),p1, height = 8, width = 7, units="in" )

#####Figure 2 Heatmap and DEA#####
#MaAsLin2 for significant species
#MaAsLin for taxa 
library(Maaslin2)

meta 
table(meta$Diabetes)
table(meta$LIBRA.score.quartile, useNA = "ifany")
table(meta$Age.baseline, useNA = "ifany")
table(meta$Race, useNA = "ifany")
############
otu_species <- t(otu_table(physeq_mindsource_final))
df_maslin <- as.data.frame(otu_species)
colnames(df_maslin) <- as.data.frame(tax_table(physeq_mindsource_final))$Species
df_maslin = t(df_maslin)
colnames(df_maslin) <- sub("^X", "", colnames(df_maslin))
input_data = (df_maslin)
table(meta$LIBRA.score, useNA = "ifany")
names(meta)


run_maaslin <- function(df, meta, output, formula, min_prevalence, min_abundance,
                        plot_associations = FALSE, ...) {
  maaslin3(
    input_data = df, input_metadata = meta, output = output, formula = formula,
    min_prevalence = min_prevalence, min_abundance = min_abundance,
    plot_associations = plot_associations, ...
  )
}

#' Filter rows by abundance threshold across samples
#'
#' Keeps rows (e.g., taxa) whose value is >= `min_value` in at least
#' `min_prop_samples` of the samples.
#'
#' @param x Numeric matrix or data frame of relative abundances
#'          (rows = taxa, columns = samples).
#' @param min_value Minimum abundance a row must reach in a sample to count
#'                  (0.25 = 25% if x is proportions; use 25 if x is percentages).
#' @param min_prop_samples Proportion of samples that must meet `min_value`
#'                         (0.5 = 50% of samples).
#' @return `x` with only the rows that pass the filter.
filter_rows <- function(x, min_value = 0.25, min_prop_samples = 0.5) {
  m <- as.matrix(x)
  if (!is.numeric(m)) stop("`x` must contain only numeric values.")
  
  n_samples <- ncol(m)
  # Small tolerance guards against floating-point error (e.g., 0.3 * 10 = 3.0000000000000004)
  min_samples <- ceiling(min_prop_samples * n_samples - 1e-9)
  
  n_passing <- rowSums(m >= min_value, na.rm = TRUE)
  keep <- n_passing >= min_samples
  
  message("Requiring >= ", min_value, " in at least ", min_samples, " of ",
          n_samples, " samples: ", sum(keep), " of ", nrow(m), " rows kept.")
  x[keep, , drop = FALSE]
}
input_data_rel = input_data
meta = libra.alpha.all
meta <- as.data.frame(meta)
rownames(meta) <- meta$SampleID
meta$Sex = factor(meta$Sex, levels = c("Male", "Female"))

input_data_rel_maas1a <- filter_rows(input_data_rel, min_value = 0.000001, min_prop_samples = 0.1)
input_data_rel_maas1b <- filter_rows(input_data_rel, min_value = 0.00001, min_prop_samples = 0.1)
input_data_rel_maas1c <- filter_rows(input_data_rel, min_value = 0.0001, min_prop_samples = 0.1)
input_data_rel_maas1d <- filter_rows(input_data_rel, min_value = 0.001, min_prop_samples = 0.1)
input_data_rel_maas1e <- filter_rows(input_data_rel, min_value = 0.01, min_prop_samples = 0.1)
input_data_rel_maas1e <- filter_rows(input_data_rel, min_value = 0.1, min_prop_samples = 0.1)


input_data_rel_maas <- filter_rows(input_data_rel, min_value = 0.00001, min_prop_samples = 0.25)
input_data_rel_maas2 <- filter_rows(input_data_rel, min_value = 0.0001, min_prop_samples = 0.25)
input_data_rel_maas3 <- filter_rows(input_data_rel, min_value = 0.0001, min_prop_samples = 0.5)
input_data_rel_maas4 <- filter_rows(input_data_rel, min_value = 0.001, min_prop_samples = 0.5)


fit_species_filter <- run_maaslin(input_data_rel, meta, paste0("MaAsLin3_full_model_no_confoudners_", Sys.Date()),
                                  "~LIBRA.score", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_species_filter <- run_maaslin(input_data_rel_maas, meta, paste0("MaAsLin3_full_model_no_confoudners_", Sys.Date()),
                                  "~LIBRA.score", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_species_filter2 <- run_maaslin(input_data_rel_maas2, meta, paste0( "MaAsLin3_full_model_no_confoudners2_", Sys.Date()),
                                   "~LIBRA.score", min_prevalence = 0, min_abundance = 0, max_significance=0.2)

fit_species_filter3 <- run_maaslin(input_data_rel_maas3, meta, paste0( "MaAsLin3_full_model_no_confoudners3_", Sys.Date()),
                                   "~LIBRA.score", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_species_filter4 <- run_maaslin(input_data_rel_maas4, meta, paste0("MaAsLin3_full_model_no_confoudners4_", Sys.Date()),
                                   "~LIBRA.score", min_prevalence = 0, min_abundance = 0,max_significance=0.2)




fit_species_filter <- run_maaslin(input_data_rel, meta, paste0("MaAsLin3_full_model_confoudners_", Sys.Date()),
                                  "~LIBRA.score + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_species_filter <- run_maaslin(input_data_rel_maas, meta, paste0("MaAsLin3_full_model_confoudners_", Sys.Date()),
                                  "~LIBRA.score + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_species_filter2 <- run_maaslin(input_data_rel_maas2, meta, paste0( "MaAsLin3_full_model_confoudners2_", Sys.Date()),
                                   "~LIBRA.score + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0, max_significance=0.2)

fit_species_filter3 <- run_maaslin(input_data_rel_maas3, meta, paste0( "MaAsLin3_full_model_confoudners3_", Sys.Date()),
                                   "~LIBRA.score*Sex + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_species_filter4 <- run_maaslin(input_data_rel_maas4, meta, paste0("MaAsLin3_full_model_confoudners4_", Sys.Date()),
                                   "~LIBRA.score + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)



#####test
fit_species_filter1a <- run_maaslin(input_data_rel_maas1a, meta, paste0("MaAsLin3_full_model_confoudners1a_", Sys.Date()),
                                  "~LIBRA.score*Sex + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_species_filter1b <- run_maaslin(input_data_rel_maas1b, meta, paste0("MaAsLin3_full_model_confoudners1b_", Sys.Date()),
                                    "~LIBRA.score*Sex + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_species_filter1c <- run_maaslin(input_data_rel_maas1c, meta, paste0("MaAsLin3_full_model_confoudners1c_", Sys.Date()),
                                    "~LIBRA.score*Sex + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)
fit_species_filter1d <- run_maaslin(input_data_rel_maas1d, meta, paste0("MaAsLin3_noBMI_model_confoudners1d_", Sys.Date()),
                                    "~LIBRA.score*Sex + BMI+  Age +group(Homeless.ever)", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_species_filter1d_sex <- run_maaslin(input_data_rel_maas1d, meta, paste0("MaAsLin3_Sex_strat_model_confoudners1d_", Sys.Date()),
                                    "~LIBRA.score+Sex + BMI+  Age +group(Homeless.ever)", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_species_filter1e <- run_maaslin(input_data_rel_maas1e, meta, paste0("MaAsLin3_full_model_confoudners1e_", Sys.Date()),
                                    "~LIBRA.score*Sex + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)


library(dplyr)
library(flextable)
library(officer)

# ---- 1. Prep the data ----
fmt_p <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "<0.001",
    p < 0.01  ~ formatC(p, format = "f", digits = 3),
    TRUE      ~ formatC(p, format = "f", digits = 2)
  )
}

p_cut <- 0.05
abund_table <- fit_species_filter1d$fit_data_abundance$results
prev_tabel = fit_species_filter1d$fit_data_prevalence$results
full_table = dplyr::bind_rows(abund_table,prev_tabel) %>% mutate(model = case_when(model == "linear" ~"Abundnace",
                                                                                   T~"Prevalence")) %>%
  dplyr::filter(feature != "s__Other")

full_table_update <- full_table |>
  mutate(
    ci_lwr = coef - qt(.975, df = N_not_zero - 6) * stderr,
    ci_upr = coef + qt(.975, df = N_not_zero - 6) * stderr
  )


tbl_dat <- full_table_update %>%
  dplyr::select(feature, name, coef, stderr, pval_individual, qval_individual,model,ci_lwr,ci_upr) %>%
  dplyr::filter(!is.na(pval_individual), pval_individual < p_cut & name=="LIBRA.score") %>%
  arrange(desc(coef)) %>%
  mutate(
    coef_ci = sprintf("%.3f [%.3f, %.3f]", coef, ci_lwr,ci_upr),
    # tidy taxon labels (e.g., "s__Genus_species" -> "Genus species")
    feature = feature %>%
      gsub("^s__", "", .) %>%
      gsub("[_.]", " ", .),
    coef  = formatC(coef,   format = "f", digits = 3),
    stderr = formatC(stderr, format = "f", digits = 3),
    pval_individual = fmt_p(pval_individual),
    qval_individual = fmt_p(qval_individual)
  )%>%
  dplyr::select(feature, name, coef_ci,stderr, pval_individual, qval_individual,model)

# ---- 2. Build the flextable ----
sup_table_mad <- flextable(tbl_dat) %>%
  set_header_labels(
    feature         = "Species",
    name            = "Term",
    coef_ci            = "Coefficient [95% CI]",
    stderr          = "Standard Error",
    pval_individual = "p-value",
    qval_individual = "q-value (FDR)",
    model = "Model"
  ) %>%
  set_caption(
    caption = "Supplementary Table X. Species significantly associated with LIBRA Score."
  ) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  italic(j = "feature", part = "body") %>%       # italicize species names
  align(j = c("coef_ci", "stderr", "pval_individual", "qval_individual"),
        align = "center", part = "all") %>%
  align(j = c("feature", "name","model"), align = "left", part = "all") %>%
  # three-line (journal-style) borders
  border_remove() %>%
  hline_top(border = fp_border(width = 1.5), part = "all") %>%
  hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  hline(border = fp_border(width = 0.75), part = "header") %>%
  autofit()

# ---- 3. Save ----
save_as_docx(
  sup_table_mad,
  path = paste0("Supplementary_Table_libra_taxa_", Sys.Date(), ".docx")
)



####LIBRA x SEx
tbl_dat2 <- full_table_update %>%
  dplyr::select(feature, name, coef, stderr, pval_individual, qval_individual,model,ci_lwr,ci_upr) %>%
  dplyr::filter(!is.na(pval_individual), pval_individual < p_cut & name=="LIBRA.score:SexFemale") %>%
  arrange(desc(coef)) %>%
  mutate(
    coef_ci = sprintf("%.3f [%.3f, %.3f]", coef, ci_lwr,ci_upr),
    # tidy taxon labels (e.g., "s__Genus_species" -> "Genus species")
    feature = feature %>%
      gsub("^s__", "", .) %>%
      gsub("[_.]", " ", .),
    coef  = formatC(coef,   format = "f", digits = 3),
    stderr = formatC(stderr, format = "f", digits = 3),
    pval_individual = fmt_p(pval_individual),
    qval_individual = fmt_p(qval_individual)
  )%>%
  dplyr::select(feature, name, coef_ci,stderr, pval_individual, qval_individual,model)

# ---- 2. Build the flextable ----
sup_table_mad2 <- flextable(tbl_dat2) %>%
  set_header_labels(
    feature         = "Species",
    name            = "Term",
    coef_ci            = "Coefficient",
    stderr          = "Standard Error",
    pval_individual = "p-value",
    qval_individual = "q-value (FDR)",
    model = "Model"
  ) %>%
  set_caption(
    caption = "Supplementary Table X. Species significantly associated with LIBRA Score x Sex Interaction."
  ) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  italic(j = "feature", part = "body") %>%       # italicize species names
  align(j = c("coef_ci", "stderr", "pval_individual", "qval_individual"),
        align = "center", part = "all") %>%
  align(j = c("feature", "name","model"), align = "left", part = "all") %>%
  # three-line (journal-style) borders
  border_remove() %>%
  hline_top(border = fp_border(width = 1.5), part = "all") %>%
  hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  hline(border = fp_border(width = 0.75), part = "header") %>%
  autofit()

# ---- 3. Save ----
save_as_docx(
  sup_table_mad2,
  path = paste0("Supplementary_Table_libra_sex_taxa_", Sys.Date(), ".docx")
)


######filter for SCFA bacteria
##########
clean_taxa <- function(x) {
  x <- sub("^.*s__", "", as.character(x))  # drop "s__" (and any lineage before it)
  x <- gsub("_", " ", x)                   # underscores -> spaces
  trimws(x)
}

filter_scfa <- function(df, scfa_list) {
  # If the organism names are in a first column instead of rownames, move them to rownames
  if (!is.numeric(df[[1]])) {
    rownames(df) <- make.unique(as.character(df[[1]]))
    df <- df[, -1, drop = FALSE]
  }
  
  cleaned <- clean_taxa(rownames(df))
  keep <- tolower(cleaned) %in% tolower(clean_taxa(scfa_list))
  
  out <- df[keep, , drop = FALSE]
  rownames(out) <- cleaned[keep]
  
  # Report anything in the list that wasn't found
  missing <- setdiff(tolower(clean_taxa(scfa_list)), tolower(cleaned))
  if (length(missing) > 0) message(length(missing), " scfa_list entries not found in data")
  
  out
}

scfa_df <- filter_scfa(input_data_rel, scfa_list)




fit_scfa_species_filter <- run_maaslin(scfa_df, meta, paste0("MaAsLin3_scfa_confoudners_", Sys.Date()),
                                  "~LIBRA.score + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_species_filter <- run_maaslin(input_data_rel_maas, meta, paste0("MaAsLin3_full_model_confoudners_", Sys.Date()),
                                  "~LIBRA.score + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_species_filter2 <- run_maaslin(input_data_rel_maas2, meta, paste0( "MaAsLin3_full_model_confoudners2_", Sys.Date()),
                                   "~LIBRA.score + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0, max_significance=0.2)

fit_species_filter3 <- run_maaslin(input_data_rel_maas3, meta, paste0( "MaAsLin3_full_model_confoudners3_", Sys.Date()),
                                   "~LIBRA.score + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_species_filter4 <- run_maaslin(input_data_rel_maas4, meta, paste0("MaAsLin3_full_model_confoudners4_", Sys.Date()),
                                   "~LIBRA.score + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)





fit_species_filter <- run_maaslin(input_data_rel, meta, paste0("MaAsLin3_full_model_confoudners_", Sys.Date()),
                                  "~LIBRA.score + Age +group(Homeless.ever)+ BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)





gen.fit_LIBRA <- Maaslin2(input_data = input_data, input_metadata = meta, min_prevalence = 0.5, 
                          min_abundance = 0.1, normalization = "TSS", 
                          output = "MaAsLin2_Out_final", 
                          fixed_effects = c("LIBRA.score","Age","BMI"),
                          plot_scatter = F,
                          plot_heatmap = F,
                          correction="BH",
                          analysis_method = "LM",
                          transform = "LOG",
                          max_significance = 0.25)

library(maaslin3)


gen.fit_LIBRA = maaslin3(
  input_data = input_data_rel_maas,
  input_metadata = meta,
  #output = paste0(new_dirc, "/", name, "/", name, "_MaAsLin3_output"),
  output= paste0("MaAsLin3_full_model_no_confoudners_",Sys.Date()),
  #formula = "~group_cat+age+group(race_cat)+group(Picogreen_cat)+BMI",
  formula = "~LIBRA.score",
  min_prevalence = 0,
  min_abundance = 0,
  plot_associations = F,
  max_significance = 0.2)

gen.fit_LIBRA = maaslin3(
  input_data = input_data_rel_maas,
  input_metadata = meta,
  #output = paste0(new_dirc, "/", name, "/", name, "_MaAsLin3_output"),
  output= paste0("MaAsLin3_full_model_confoudners_",Sys.Date()),
  #formula = "~group_cat+age+group(race_cat)+group(Picogreen_cat)+BMI",
  formula = "~LIBRA.score + Age +group(Homeless.ever)",
  min_prevalence = 0,
  min_abundance = 0,
  plot_associations = F,
  max_significance = 0.2)


gen.fit_LIBRA = maaslin3(
  input_data = input_data_rel_maas,
  input_metadata = meta,
  #output = paste0(new_dirc, "/", name, "/", name, "_MaAsLin3_output"),
  output= paste0("MaAsLin3_full_model_confoudners_",Sys.Date()),
  #formula = "~group_cat+age+group(race_cat)+group(Picogreen_cat)+BMI",
  formula = "~LIBRA.score.quartile + Age +group(Homeless.ever)",
  min_prevalence = 0,
  min_abundance = 0,
  plot_associations = F,
  max_significance = 0.2)


gen.fit_LIBRA <- Maaslin2(input_data = input_data_rel_maas, input_metadata = meta, min_prevalence = 0, 
                          min_abundance = 0, normalization = "TSS", 
                          output = "MaAsLin2_Out_homless-final-10-5-26", 
                          fixed_effects = c("LIBRA.score","Age","BMI","Homeless.ever"),
                          plot_scatter = F,
                          plot_heatmap = F,
                          correction="BH",
                          analysis_method = "LM",
                          transform = "LOG",
                          max_significance = 0.3)


#####Heatmap for relative abundances
mal3_input = read.table('/Users/mattigo/Documents/CU_Folders/Projects_Ongoing/NFLPA_Microbiome/Test_folder/MaAsLin3_noBMI_model_confoudners1d_2026-10-30/all_results.tsv', header = T)
libra_mas = data.frame(mal3_input)
mas.sig = libra_mas %>% dplyr::filter(pval_individual < 0.05 & name == "LIBRA.score:SexFemale")

input_data2 = data.frame(input_data)
colnames(input_data2) <- sub("^X", "", colnames(input_data2))

relative_taxa_adbund_filtered = input_data2 %>% dplyr::filter(row.names(.) %in% mas.sig$feature)

relative_taxa_adbund_filtered_update <- relative_taxa_adbund_filtered %>%
  rownames_to_column(var = "Species")

relative_taxa_adbund_filtered_update2 <- data.frame(t(relative_taxa_adbund_filtered_update[-1]))
colnames(relative_taxa_adbund_filtered_update2) <- relative_taxa_adbund_filtered_update[, 1]

relative_taxa_adbund_filtered_update3 <- relative_taxa_adbund_filtered_update2 %>%
  rownames_to_column(var = "SampleID")


relative_taxa_adbund_filtered_update3_melt = melt(relative_taxa_adbund_filtered_update3, id.var="SampleID", variable.name="Species", value.name = "Abundance") 


all_taxa_meta_data = plyr::join(relative_taxa_adbund_filtered_update3_melt, meta.in, by="SampleID") %>% 
  mutate(Species = str_remove(Species, "s__")) %>% 
  mutate(Species = str_replace(Species, "_", " "))

all_taxa_meta_data2 = all_taxa_meta_data %>% mutate(Abundance_log = log10(Abundance + 1E-10))
all_taxa_meta_data2$LIBRA.score.quartile = factor(all_taxa_meta_data2$LIBRA.score.quartile, levels = c("4", "3", "2", "1"))
all_taxa_meta_data3 = all_taxa_meta_data2 %>% arrange(Species, (LIBRA.score)) %>%
  mutate(Libra.score.half = case_when(LIBRA.score.quartile == "1" | LIBRA.score.quartile == "2" ~ "Lower",
                                      T~"Upper"))

gen.fit_LIBRA
mas.sig2 = mas.sig %>% mutate(Up_Down = case_when(coef > 0 ~ "Up",
                                                  T ~ "Down")) %>%
  arrange((pval_individual))%>% 
  mutate(feature = str_remove(feature, "s__")) %>% 
  mutate(feature = str_replace(feature, "_", " "))

mas.sig2$Up_Down = factor(mas.sig2$Up_Down, levels = c("Up", "Down")) 
#Creates figure for horizontal LFC
mas_sig_list = mas.sig2 %>% arrange(desc(coef))


all_taxa_meta_data_test = all_taxa_meta_data2 %>% 
  arrange(factor(Species, levels =  mas_sig_list$feature), LIBRA.score)

colvec <- c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A")


cols = list( `LIBRA Score Quartile`  = c( "4" = "#4A4A4A", "3" = "#CD3333","2" = "#FFC125","1" = "#4F94CD" ))



#####Final heatmap
all_taxa_meta_data_test = all_taxa_meta_data_test %>%      
  rename(
    `LIBRA Score` = LIBRA.score,
    `LIBRA Score Quartile` = LIBRA.score.quartile)

tdy_hm_rel2=tidyheatmap(df = all_taxa_meta_data_test,
                        rows = Species,
                        columns = SampleID,
                        values = Abundance_log,
                        #cluster_cols = T,
                        #clustering_distance_cols = "euclidean",
                        #clustering_method = "complete",
                        annotation_col = c(`LIBRA Score Quartile`,  `LIBRA Score`, `Sex`) ,
                        annotation_colors = cols,
                        colors = c("#ffffff","#145afc"),
                        scale= "none",
                        display_numbers = F,
                        gaps_col = `LIBRA Score Quartile`,
                        show_colnames = F)
dev.off()
#adjusts the font in the y.axis labels
library(grid)
tdy_hm_rel2$gtable$grobs[[2]]$gp = gpar(fontsize = 14, fontface="italic")



library(ggplotify)
#Turns the heatmap to a ggplot item to make it easier to work with
test_hm = as.ggplot(tdy_hm_rel2) 
test_hm = test_hm +# labs(tag = "B") +
  annotate("text", x = 0.325, y = 0.00001, label = "Samples") +  
  coord_cartesian(ylim = c(0, 1), clip = "off")+
  theme(plot.margin = margin(t = 0, r = -2.4, b = 0.3, l = 0, unit = "cm"),
        plot.tag = element_text(size = 12, face = "bold"))
test_hm
ggsave(paste0("Figures/Figure_2_LIBRA*SEX_INTERACTION_Quartile_Heatmap_",Sys.Date(),".jpeg"), test_hm, jpeg, height = 8, width = 16, units = "in" )


###Sex heatmap
all_taxa_meta_data_test = all_taxa_meta_data2 %>% 
  arrange(factor(Species, levels =  mas_sig_list$feature), Sex, LIBRA.score)
all_taxa_meta_data_test = all_taxa_meta_data_test %>%      
  rename(
    `LIBRA Score` = LIBRA.score,
    `LIBRA Score Quartile` = LIBRA.score.quartile)


tdy_hm_rel2_sex=tidyheatmap(df = all_taxa_meta_data_test,
                        rows = Species,
                        columns = SampleID,
                        values = Abundance_log,
                        #cluster_cols = T,
                        #clustering_distance_cols = "euclidean",
                        #clustering_method = "complete",
                        annotation_col = c(`Sex`,`LIBRA Score Quartile`,  `LIBRA Score`) ,
                        annotation_colors = cols,
                        colors = c("#ffffff","#145afc"),
                        scale= "none",
                        display_numbers = F,
                        gaps_col = `LIBRA Score Quartile`,
                        show_colnames = F)
dev.off()
#adjusts the font in the y.axis labels
library(grid)
tdy_hm_rel2_sex$gtable$grobs[[2]]$gp = gpar(fontsize = 8, fontface="italic")



library(ggplotify)
#Turns the heatmap to a ggplot item to make it easier to work with
test_hm_sex = as.ggplot(tdy_hm_rel2_sex) 
test_hm_sex = test_hm_sex +# labs(tag = "B") +
  annotate("text", x = 0.325, y = 0.00001, label = "Samples") +  
  coord_cartesian(ylim = c(0, 1), clip = "off")+
  theme(plot.margin = margin(t = 0, r = -2.4, b = 0.3, l = 0, unit = "cm"),
        plot.tag = element_text(size = 12, face = "bold"))
test_hm_sex

set_text_size <- function(g, size, face = NULL) {
  if (inherits(g, "text")) {
    g$gp$fontsize <- size
    if (!is.null(face)) g$gp$fontface <- face
    return(g)
  }
  if (inherits(g, "gTree")) {
    for (i in seq_along(g$children)) {
      g$children[[i]] <- set_text_size(g$children[[i]], size, face)
    }
  }
  g
}

stretch_y <- function(g, k) {
  if (inherits(g, "gTree")) {
    for (i in seq_along(g$children)) g$children[[i]] <- stretch_y(g$children[[i]], k)
    return(g)
  }
  if (!is.null(g$y)) {
    grid.newpage()
    pushViewport(viewport(height = unit(1000, "pt")))
    yn <- convertY(g$y, "npc", valueOnly = TRUE)   # fractions of a 1000-pt viewport
    popViewport()
    g$y <- unit(1, "npc") - unit(k * (1 - yn) * 1000, "pt")
  }
  g
}


gt <- tdy_hm_rel2_sex$gtable

# Color-scale legend and annotation legend
for (nm in c("legend", "annotation_legend")) {
  idx <- which(gt$layout$name == nm)
  if (length(idx)) {
    gt$grobs[[idx]] <- set_text_size(gt$grobs[[idx]], size = 11)
    # widen the legend column so larger text isn't clipped
    col <- gt$layout$l[idx]
    gt$widths[col] <- gt$widths[col] * (11 / 10)   # 10 = pheatmap's default fontsize
  }
}

# See what's in there (run once)
# gt$layout$name

# --- y-axis (row) labels: bigger, bold italic ---
rn <- which(gt$layout$name == "row_names")
gt$grobs[[rn]]$gp <- gpar(fontsize = 14, fontface = "bold.italic")

# Widen the column holding the labels so they don't hit the legend
gt$widths[gt$layout$l[rn]] <- grobWidth(gt$grobs[[rn]]) + unit(10, "pt")
gt$layout$name
# --- annotation bar names at the top: bigger ---
an <- which(gt$layout$name == "col_annotation_names")
if (length(an)) {
  gt$grobs[[an]]$gp <- gpar(fontsize = 11, fontface="bold")
  # annotation names share the label column, so make sure it fits the longest one
  gt$widths[gt$layout$l[an]] <- max(
    gt$widths[gt$layout$l[rn]],
    grobWidth(gt$grobs[[an]]) + unit(10, "pt")
  )
}

# --- annotation legend text (optional, keeps it proportional) ---
lg <- which(gt$layout$name == "annotation_legend")
if (length(lg)) gt$grobs[[lg]]$gp <- gpar(fontsize = 12)

#k_ann <- 1.1
#for (nm in c("col_annotation", "annotation_names_col")) {
#  i <- which(gt$layout$name == nm)
#  if (length(i)) gt$grobs[[i]] <- stretch_y(gt$grobs[[i]], k_ann)
#}
#ia <- which(gt$layout$name == "col_annotation")
#if (length(ia)) {
#  r <- gt$layout$t[ia]
#  gt$heights[r] <- gt$heights[r] * k_ann      # give the stretched bars room
#}

# ---- legends (color scale + annotation legend) ----
k_leg <- 1.1
for (nm in c("legend", "annotation_legend")) {
  i <- which(gt$layout$name == nm)
  if (length(i)) gt$grobs[[i]] <- stretch_y(gt$grobs[[i]], k_leg)
}

tdy_hm_rel2_sex$gtable <- gt

library(ggplotify)
test_hm_sex <-  as.ggplot(tdy_hm_rel2_sex) +
  annotate("text", x = 0.325, y = 0.00000004, label = "Samples", fontface="bold", size = 8,vjust = 1) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  theme(plot.margin = margin(t = 0, r = 1, b = 1.0, l = 3.5, unit = "cm"))
test_hm_sex

ggsave(paste0("Figures/Figure_2_LIBRA*SEX_INTERACTION_SEX_REORDERED_Quartile_Heatmap_",Sys.Date(),".jpeg"), test_hm_sex, jpeg, height = 8, width = 16, units = "in" )


#####figure 2b
######LFC of significant species####3
gen.fit_LIBRA
mas.sig2 = mas.sig %>% mutate(Up_Down = case_when(coef > 0 ~ "Up",
                                                  T ~ "Down")) %>%
  arrange((qval))%>% 
  mutate(feature = str_remove(feature, "s__")) %>% 
  mutate(feature = str_replace(feature, "_", " "))

mas.sig2$Up_Down = factor(mas.sig2$Up_Down, levels = c("Up", "Down")) 
#Creates figure for horizontal LFC
mas_sig_list = mas.sig2 %>% arrange(desc(coef))

p3 = ggplot(mas.sig2, aes(x=reorder(feature,coef), y=coef, fill=Up_Down)) +
  geom_bar(stat='identity') +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))+
  scale_x_discrete(position = "top")+
  #theme(axis.text.y = element_text(face = "italic")) + # Assign specific colors to groups
  coord_flip() +
  labs(#title="Species LFC by LIBRA Score", 
    x="Species", 
    y="Log Fold Change") +
  scale_fill_manual(values = c("Up" = "darkred", "Down" = "skyblue"))

ggsave(paste0("Figures/Figure_2a_LFC_Significant_Species_",Sys.Date(),".jpeg"), p3, jpeg, height = 8, width = 12, units = "in" )


p3 = p3 + labs(tag = "A") +
  theme(plot.tag = element_text(size = 12, face = "bold"))
#Merge the heatmap and LFC figures together
p2 = ggarrange(p3, test_hm, ncol=1)
ggsave(paste0("Figures/Figure_2_Significant_Species_LFC_HM_",Sys.Date(),".jpeg"), p2, jpeg, height = 10, width = 12, units = "in" )

mas.sig2_outtable = mas.sig2 %>% dplyr::select(feature, coef, stderr, pval, qval) %>% arrange(desc(coef))

sup_table_2 = flextable(mas.sig2_outtable)
sup_table_2 = set_caption(sup_table_2, caption = "Supplementary Table 2. Species associated with LIBRA Scores") %>% autofit()
save_as_docx(sup_table_2, path = "Figures/Supplementary_Table_2_MaAsLin_Results.docx")






###########LIBRA only heatmap
libra_mas = data.frame(mal3_input)
mas.sig = libra_mas %>% dplyr::filter(pval_individual < 0.05 & name == "LIBRA.score")


input_data2 = data.frame(input_data)
colnames(input_data2) <- sub("^X", "", colnames(input_data2))

relative_taxa_adbund_filtered = input_data2 %>% dplyr::filter(row.names(.) %in% mas.sig$feature)

relative_taxa_adbund_filtered_update <- relative_taxa_adbund_filtered %>%
  rownames_to_column(var = "Species")

relative_taxa_adbund_filtered_update2 <- data.frame(t(relative_taxa_adbund_filtered_update[-1]))
colnames(relative_taxa_adbund_filtered_update2) <- relative_taxa_adbund_filtered_update[, 1]

relative_taxa_adbund_filtered_update3 <- relative_taxa_adbund_filtered_update2 %>%
  rownames_to_column(var = "SampleID")


relative_taxa_adbund_filtered_update3_melt = melt(relative_taxa_adbund_filtered_update3, id.var="SampleID", variable.name="Species", value.name = "Abundance") %>%
  dplyr::filter(Species !="s__Other")



all_taxa_meta_data = plyr::join(relative_taxa_adbund_filtered_update3_melt, meta.in, by="SampleID") %>% 
  mutate(Species = str_remove(Species, "s__")) %>% 
  mutate(Species = str_replace(Species, "_", " "))

all_taxa_meta_data2 = all_taxa_meta_data %>% mutate(Abundance_log = log10(Abundance + 1E-10))
all_taxa_meta_data2$LIBRA.score.quartile = factor(all_taxa_meta_data2$LIBRA.score.quartile, levels = c("4", "3", "2", "1"))
all_taxa_meta_data3 = all_taxa_meta_data2 %>% arrange(Species, (LIBRA.score)) %>%
  mutate(Libra.score.half = case_when(LIBRA.score.quartile == "1" | LIBRA.score.quartile == "2" ~ "Lower",
                                      T~"Upper"))

gen.fit_LIBRA
mas.sig2 = mas.sig %>% mutate(Up_Down = case_when(coef > 0 ~ "Up",
                                                  T ~ "Down")) %>%
  arrange((pval_individual))%>% 
  mutate(feature = str_remove(feature, "s__")) %>% 
  mutate(feature = str_replace(feature, "_", " "))

mas.sig2$Up_Down = factor(mas.sig2$Up_Down, levels = c("Up", "Down")) 
#Creates figure for horizontal LFC
mas_sig_list = mas.sig2 %>% arrange(desc(coef))


all_taxa_meta_data_test = all_taxa_meta_data2 %>% 
  arrange(factor(Species, levels =  mas_sig_list$feature), LIBRA.score)

colvec <- c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A")


cols = list( `LIBRA Score Quartile`  = c( "4" = "#4A4A4A", "3" = "#CD3333","2" = "#FFC125","1" = "#4F94CD" ))



#####Final heatmap
all_taxa_meta_data_test = all_taxa_meta_data_test %>%      
  rename(
    `LIBRA Score` = LIBRA.score,
    `LIBRA Score Quartile` = LIBRA.score.quartile)

tdy_hm_rel2=tidyheatmap(df = all_taxa_meta_data_test,
                        rows = Species,
                        columns = SampleID,
                        values = Abundance_log,
                        #cluster_cols = T,
                        #clustering_distance_cols = "euclidean",
                        #clustering_method = "complete",
                        annotation_col = c(`LIBRA Score Quartile`,  `LIBRA Score`, `Sex`) ,
                        annotation_colors = cols,
                        colors = c("#ffffff","#145afc"),
                        scale= "none",
                        display_numbers = F,
                        gaps_col = `LIBRA Score Quartile`,
                        show_colnames = F)
dev.off()
#adjusts the font in the y.axis labels
library(grid)
tdy_hm_rel2$gtable$grobs[[2]]$gp = gpar(fontsize = 14, fontface="italic")



library(ggplotify)
#Turns the heatmap to a ggplot item to make it easier to work with
test_hm = as.ggplot(tdy_hm_rel2) 
test_hm = test_hm +# labs(tag = "B") +
  annotate("text", x = 0.325, y = 0.00001, label = "Samples") +  
  coord_cartesian(ylim = c(0, 1), clip = "off")+
  theme(plot.margin = margin(t = 0, r = -2.4, b = 0.3, l = 0, unit = "cm"),
        plot.tag = element_text(size = 12, face = "bold"))
test_hm



# Set font size on every text element inside a grob (works through nested gTrees)
set_text_size <- function(g, size, face = NULL) {
  if (inherits(g, "text")) {
    g$gp$fontsize <- size
    if (!is.null(face)) g$gp$fontface <- face
    return(g)
  }
  if (inherits(g, "gTree")) {
    for (i in seq_along(g$children)) {
      g$children[[i]] <- set_text_size(g$children[[i]], size, face)
    }
  }
  g
}

stretch_y <- function(g, k) {
  if (inherits(g, "gTree")) {
    for (i in seq_along(g$children)) g$children[[i]] <- stretch_y(g$children[[i]], k)
    return(g)
  }
  if (!is.null(g$y)) {
    grid.newpage()
    pushViewport(viewport(height = unit(1000, "pt")))
    yn <- convertY(g$y, "npc", valueOnly = TRUE)   # fractions of a 1000-pt viewport
    popViewport()
    g$y <- unit(1, "npc") - unit(k * (1 - yn) * 1000, "pt")
  }
  g
}


gt <- tdy_hm_rel2$gtable

# Color-scale legend and annotation legend
for (nm in c("legend", "annotation_legend")) {
  idx <- which(gt$layout$name == nm)
  if (length(idx)) {
    gt$grobs[[idx]] <- set_text_size(gt$grobs[[idx]], size = 11)
    # widen the legend column so larger text isn't clipped
    col <- gt$layout$l[idx]
    gt$widths[col] <- gt$widths[col] * (11 / 10)   # 10 = pheatmap's default fontsize
  }
}

# See what's in there (run once)
# gt$layout$name

# --- y-axis (row) labels: bigger, bold italic ---
rn <- which(gt$layout$name == "row_names")
gt$grobs[[rn]]$gp <- gpar(fontsize = 14, fontface = "bold.italic")

# Widen the column holding the labels so they don't hit the legend
gt$widths[gt$layout$l[rn]] <- grobWidth(gt$grobs[[rn]]) + unit(10, "pt")
gt$layout$name
# --- annotation bar names at the top: bigger ---
an <- which(gt$layout$name == "col_annotation_names")
if (length(an)) {
  gt$grobs[[an]]$gp <- gpar(fontsize = 11, fontface="bold")
  # annotation names share the label column, so make sure it fits the longest one
  gt$widths[gt$layout$l[an]] <- max(
    gt$widths[gt$layout$l[rn]],
    grobWidth(gt$grobs[[an]]) + unit(10, "pt")
  )
}

# --- annotation legend text (optional, keeps it proportional) ---
lg <- which(gt$layout$name == "annotation_legend")
if (length(lg)) gt$grobs[[lg]]$gp <- gpar(fontsize = 12)

#k_ann <- 1.1
#for (nm in c("col_annotation", "annotation_names_col")) {
#  i <- which(gt$layout$name == nm)
#  if (length(i)) gt$grobs[[i]] <- stretch_y(gt$grobs[[i]], k_ann)
#}
#ia <- which(gt$layout$name == "col_annotation")
#if (length(ia)) {
#  r <- gt$layout$t[ia]
#  gt$heights[r] <- gt$heights[r] * k_ann      # give the stretched bars room
#}

# ---- legends (color scale + annotation legend) ----
k_leg <- 1.1
for (nm in c("legend", "annotation_legend")) {
  i <- which(gt$layout$name == nm)
  if (length(i)) gt$grobs[[i]] <- stretch_y(gt$grobs[[i]], k_leg)
}

tdy_hm_rel2$gtable <- gt

library(ggplotify)
test_hm <- as.ggplot(tdy_hm_rel2) +
  annotate("text", x = 0.325, y = 0.00000004, label = "Samples", fontface="bold", size = 8,vjust = 1) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  theme(plot.margin = margin(t = 0, r = 1, b = 1.0, l = 3.5, unit = "cm"))
test_hm
ggsave(paste0("Figures/Figure_2_LIBRA_ONLY_Quartile_Heatmap_",Sys.Date(),".jpeg"), test_hm, jpeg, height = 8, width = 16, units = "in" )

######Figure 3#####Biplots
#next line will determine which metadata is being piped in so make sure to check before proceeding 
meta.in
meta_env_fit = libra.alpha.all %>% dplyr::select(SampleID,LIBRA.score.quartile,Coronary.heart.disease, Physical.inactivity ,Diabetes,Hypercholesterolemia, BMI,Hypertension, Depression, current_alcohol_emr, lifetime_chronic_kindey_emr, smoking_status, healthy_diet, LIBRA.score)

meta_env_fit_update = meta_env_fit %>% mutate(Coronary.heart.disease = if_else(Coronary.heart.disease == "Yes", 1 , 0),
                                              Physical.activity.H.L = if_else(Physical.inactivity == "High", 1, 0), 
                                              Diabetes = if_else(Diabetes == "Yes", 1,0),
                                              Hypercholesterolemia = if_else(Hypercholesterolemia=="Yes",1,0),
                                              Hypertension = if_else(Hypertension=="Yes", 1,0),
                                              Current_alcohol = if_else(current_alcohol_emr =="Yes",1,0),
                                              Chronic_kidney_disease = if_else(lifetime_chronic_kindey_emr=="Yes",1,0),
                                              Lifetime_smoking = if_else(smoking_status=="Yes",1,0),
                                              Healthy_diet = if_else(healthy_diet=="Yes",1,0))
#matching up metadata with the distance matrix
x_beta
ind <- match(x = rownames(x_beta), table = meta_env_fit_update$SampleID)
ind1 <- !is.na(ind)
x1 <- x_beta[ind1,ind1]
ind2 <- match(x = rownames(x1), table = meta_env_fit_update$SampleID)
meta.in1 <- meta_env_fit_update[ind2,]
#this next line should be all true - if not, then something is wrong
table(rownames(x1) == meta.in1$SampleID) 
rownames(x1)
colnames(x1)
meta.in1$SampleID

#read in taxa
#This needs to get renamed
taxa_species
input_data_rel = data.frame(input_data)
colnames(input_data_rel) <- sub("^X", "", colnames(input_data_rel))

input_data_rel$names = rownames(input_data_rel) 

input_data_rel = input_data_rel %>% mutate(names = str_remove(names, "s__")) %>% 
  mutate(names = str_replace(names, "_", " "))
rownames(input_data_rel) = input_data_rel$names
input_data_rel = input_data_rel %>% dplyr::select(!names)
Taxa <- t(input_data_rel)
Taxa[1:10,1:10]

ind1 <- match(x = colnames(x1), table = rownames(Taxa))
Taxa1 <- Taxa[ind1,]
Taxa1[1:10,1:10]

#**Special plot** CCA start
#this code is for coloring dots on biplot, if not wanting that then next 3 lines can be removed
meta.in1$LIBRA.score.quartile <- factor(meta.in1$LIBRA.score.quartile, levels=c("4", "3", "2","1"))
with(meta.in1, levels(LIBRA.score.quartile))
scl <- 2
#colvec <- c("red", "black", "lightgreen", "royalblue")
colvec <- rev(c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A"))

#CCA start
pca_w=capscale(x1~1, data=meta.in1[,-c(1,2)])
#####envfit for libra variables
(fit1=envfit(pca_w, meta.in1[,-c(1,2)], na.rm = T))
#p.max in the code below will determine which variables end up on the graph - adjust for more or less
#plot(fit1, col="skyblue", p.max=0.3)
fit1
###envfit for species
#p.max in the code below will determine which taxa end up on the graph - adjust for more or less
Taxa1 = as.data.frame(Taxa1)
fit2=envfit(pca_w, Taxa1, na.rm = T)
#plot(fit2, col="darkred",p.max=0.002)
fit2

####Pulls together data from pca and envfit
pca_data_mds = data.frame(vegan::scores(pca_w)$sites)
#pca_data_mds = pca_data %>% select(MDS1,MDS2)
pca_data_mds$SampleID = rownames(pca_data_mds)


libra_fit = data.frame(vegan::scores(fit1, "vectors")) *ordiArrowMul(fit1)
libra_pca_pval = data.frame(fit1$vectors$pvals) %>% dplyr::filter(fit1.vectors.pvals < 0.3)
libra_fit_filter <- libra_fit %>%
  dplyr::filter(row.names(.) %in% rownames(libra_pca_pval)) %>%
  mutate(indicator= "libra")


species_fit = data.frame(vegan::scores(fit2, "vectors")) *ordiArrowMul(fit2)
species_pca_pval = data.frame(fit2$vectors$pvals) %>% dplyr::filter(fit2.vectors.pvals <= 0.002)
species_fit_filter <- species_fit %>%
  dplyr::filter(row.names(.) %in% rownames(species_pca_pval)) %>%
  mutate(indicator= "species")


pca_data_mds = data.frame(vegan::scores(pca_w)$sites)
#pca_data_mds = pca_data %>% select(MDS1,MDS2)
pca_data_mds$SampleID = rownames(pca_data_mds)
pca_data_full = full_join(meta_env_fit_update,pca_data_mds, by="SampleID" )
pca_data_full$LIBRA.score.quartile = factor(pca_data_full$LIBRA.score.quartile, levels=c("4","3","2","1"))

colvec <- rev(c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A"))


####ggplot for species biplot
species_biplot= ggplot(pca_data_full, aes(x=MDS1, y=MDS2))+
  geom_point(aes(color= LIBRA.score.quartile))+
  scale_colour_manual(values = colvec)+
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), 
               data = libra_fit_filter, size =1, alpha = 0.5, colour = "skyblue")+
  ggrepel::geom_text_repel(max.overlaps=nrow(libra_fit_filter),data = libra_fit_filter, aes(x = MDS1, y = MDS2), 
                           label = row.names(libra_fit_filter), colour = "skyblue",nudge_x=-0.05,nudge_y=0.05)+
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), 
               data = species_fit_filter, size =1, alpha = 0.5, colour = "darkred")+
  ggrepel::geom_text_repel(max.overlaps=15,data = species_fit_filter, aes(x = MDS1, y = MDS2), 
                           label = row.names(species_fit_filter), colour = "darkred",nudge_y=-0.03,nudge_x=0.01) 

species_biplot = species_biplot + labs(tag = "B") +
  theme(plot.tag = element_text(size = 12, face = "bold"))

ggsave(paste0("Figures/Figure_3b_Biplot_Species_",Sys.Date(),".jpeg"), species_biplot, jpeg, height = 8, width = 12, units = "in" )



#######Takes bray-curtis distance for phylum level biplot
ps_phylum <- tax_glom(physeq_mindsource_final, taxrank = "Phylum", NArm = TRUE)
otu_phylum <- t(otu_table(ps_phylum))
df_phylum <- as.data.frame(otu_phylum)
colnames(df_phylum) <- as.data.frame(tax_table(ps_phylum))$Phylum
df_phylum = t(df_phylum)
colnames(df_phylum) <- sub("^X", "", colnames(df_phylum))

input_data_phylum = data.frame(df_phylum)
colnames(input_data_phylum) <- sub("^X", "", colnames(input_data_phylum))

#Removing virus and archae samples from analysis
input_data_phylum = input_data_phylum[!(row.names(input_data_phylum) %in% c("p__Viruses_noname","p__Eukaryota_noname", "p__Euryarchaeota")), ]

trans_taxa_table_1 = t(data.frame(input_data_phylum))
trans_taxa_table_1 = data.frame(trans_taxa_table_1)
trans_taxa_table_1$names = rownames(trans_taxa_table_1)
trans_taxa_table_1 = trans_taxa_table_1 %>% mutate(names = str_replace(names, "X", "")) 
rownames(trans_taxa_table_1) = trans_taxa_table_1$names
trans_taxa_table_1= trans_taxa_table_1 %>% dplyr::select(!names)

#matching up metadata with the distance matrix
ind <- match(x = rownames(x_beta), table = meta_env_fit_update$SampleID)
ind1 <- !is.na(ind)
x1 <- x_beta[ind1,ind1]
ind2 <- match(x = rownames(x1), table = meta_env_fit_update$SampleID)
meta.in1 <- meta_env_fit_update[ind2,]
#this next line should be all true - if not, then something is wrong
table(rownames(x1) == meta.in1$SampleID) 
rownames(x1)
colnames(x1)
meta.in1$SampleID
#read in taxa
input_data_phylum
input_data_phylum$names = rownames(input_data_phylum) 
input_data_phylum = input_data_phylum %>% mutate(names = str_remove(names, "p__")) %>% 
  mutate(names = str_replace(names, "_", " "))
rownames(input_data_phylum) = input_data_phylum$names
input_data_phylum = input_data_phylum %>% dplyr::select(!names)

Taxa <- t(input_data_phylum)
head(Taxa)

ind1 <- match(x = colnames(x1), table = rownames(Taxa))
Taxa1 <- Taxa[ind1,]
head(Taxa1)

#**Special plot** CCA start
#this code is for coloring dots on biplot, if not wanting that then next 3 lines can be removed
meta.in1$LIBRA.score.quartile <- factor(meta.in1$LIBRA.score.quartile, levels=c("4", "3", "2","1"))
with(meta.in1, levels(LIBRA.score.quartile))
scl <- 2
#colvec <- c("red", "black", "lightgreen", "royalblue")
colvec <- rev(c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A"))

#CCA start
pca_w=capscale(x1~1, data=meta.in1[,-c(1,2)])
#envfit for LIBRA factors
(fit1=envfit(pca_w, meta.in1[,-c(1,2)], na.rm = T))
#p.max in the code below will determine which variables end up on the graph - adjust for more or less
fit1
###
#envfit for phylum level
#p.max in the code below will determine which taxa end up on the graph - adjust for more or less
fit2=envfit(pca_w, Taxa1, na.rm = T)
fit2


pca_data_mds = data.frame(vegan::scores(pca_w)$sites)
#pca_data_mds = pca_data %>% select(MDS1,MDS2)
pca_data_mds$SampleID = rownames(pca_data_mds)

libra_fit = data.frame(vegan::scores(fit1, "vectors")) *ordiArrowMul(fit1)
libra_pca_pval = data.frame(fit1$vectors$pvals) %>% dplyr::filter(fit1.vectors.pvals <= 0.9)
libra_fit_filter <- libra_fit %>%
  dplyr::filter(row.names(.) %in% rownames(libra_pca_pval)) %>%
  mutate(indicator= "libra")


phylum_fit = data.frame(vegan::scores(fit2, "vectors")) *ordiArrowMul(fit2)
phylum_pca_pval = data.frame(fit2$vectors$pvals) %>% dplyr::filter(fit2.vectors.pvals <= 0.9)
phylum_fit_filter <- phylum_fit %>%
  dplyr::filter(row.names(.) %in% rownames(phylum_pca_pval)) %>%
  mutate(indicator= "species")

all_fit_lines = bind_rows(libra_fit_filter,phylum_fit_filter)
pca_data_full = full_join(meta_env_fit_update,pca_data_mds, by="SampleID" )
pca_data_full$LIBRA.score.quartile = factor(pca_data_full$LIBRA.score.quartile, levels=c("4","3","2","1"))

colvec <- rev(c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A"))

phylum_biplot = ggplot(pca_data_full, aes(x=MDS1, y=MDS2))+
  geom_point(aes(color= LIBRA.score.quartile))+
  scale_colour_manual(values = colvec)+
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), 
               data = libra_fit_filter, size =1, alpha = 0.5, colour = "skyblue")+
  ggrepel::geom_text_repel(max.overlaps=nrow(libra_fit_filter),data = libra_fit_filter, aes(x = MDS1, y = MDS2), 
                           label = row.names(libra_fit_filter), colour = "skyblue",nudge_x=-0.05,nudge_y=0.05)+
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), 
               data = phylum_fit_filter, size =1, alpha = 0.5, colour = "darkred")+
  ggrepel::geom_text_repel(max.overlaps=15,data = phylum_fit_filter, aes(x = MDS1, y = MDS2), 
                           label = row.names(phylum_fit_filter), colour = "darkred",nudge_y=0.01,nudge_x=0.01) 

phylum_biplot = phylum_biplot + labs(tag = "A") +
  theme(plot.tag = element_text(size = 12, face = "bold"))

ggsave(paste0("Figures/Figure_3a_Biplot_Species_",Sys.Date(),".jpeg"), phylum_biplot, jpeg, height = 8, width = 12, units = "in" )
biplot_organisms = ggarrange(phylum_biplot + theme(text = element_text(size = 12, face = "bold")), species_biplot+ theme(text = element_text(size = 12, face = "bold")), nrow=2)
ggsave(paste0("Figures/Figure_3_Biplot_Organisms_",Sys.Date(),".jpeg"), biplot_organisms, jpeg, height = 10, width = 12, units = "in" )


#####Pathway Analysis
#Figure 3
####Volcano
#DeSeq for genes and pathways
#Differential Abundance with DESeq2 - done at genus level
detach("package:microbiome", unload = TRUE)
detach("package:phyloseq", unload = TRUE)
library(DESeq2)
library(pheatmap)
library(org.Mm.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
library(apeglm)


path_table_1 <- read.csv("/Volumes/B4_Backup/Metagemomic-sequences/MQ3201_0001_DeepSeq_Jul_2023_Diversigen/filtered functional tables/uniref-pathway-filtered-absolute.csv", header = TRUE, row.names = 1)
path_table_1 <- path_table_1[-c(1, 2), ]
path_table_2 <- read.csv("/Volumes/B4_Backup/Metagemomic-sequences/MQ3201_0003_DeepSeq_May_2024_Diversigen/filtered functional tables/uniref-pathway-filtered-absolute.csv", header = TRUE, row.names = 1)
path_table_2 <- path_table_2[-c(1, 2), ]

path_table_all <- merge(path_table_1, path_table_2, by = 0, all = TRUE)
path_table_all[is.na(path_table_all)] <- 0
rownames(path_table_all) <- path_table_all$Row.names
path_table_all <- path_table_all %>% dplyr::select(-Row.names)
path_table_all_df <- data.frame(path_table_all)
names(path_table_all_df) <- gsub("X", "", names(path_table_all_df))
path_table_all_df_filter <- path_table_all_df[, names(path_table_all_df) %in% libra.alpha.all$SampleID]



dat_pathway = data.frame(path_table_all_df_filter)
colnames(dat_pathway) <- sub("^X", "", colnames(dat_pathway))

#dat[-1] <-sweep(dat[-1], 2, colSums(dat[,-1]), `/`) * 10000
#Round data so that all the data are integers
dat_pathway <- dat_pathway %>% 
  dplyr::mutate_if(is.numeric, round)
#keep_genes <- rowSums( dat_pathway > 100 ) >= 42
#dat_pathway2 <- dat_pathway[keep_genes,]


#dat_pathway2 <- as.matrix(dat_pathway2)
#meta

#dat_pathway
dat_pathway_filter_1 <- filter_rows(dat_pathway, min_value = 100, min_prop_samples = 0.1)

dat_pathway_filter <- filter_rows(dat_pathway, min_value = 100, min_prop_samples = 0.25)
dat_pathway_filter2 <- filter_rows(dat_pathway, min_value = 100, min_prop_samples = 0.5)
dat_pathway_filter3 <- filter_rows(dat_pathway, min_value = 1000, min_prop_samples = 0.25)
dat_pathway_filter4 <- filter_rows(dat_pathway, min_value = 1000, min_prop_samples = 0.5)
dat_pathway_filter1a <- filter_rows(dat_pathway, min_value = 1000, min_prop_samples = 0.1)
dat_pathway_filter1b <- filter_rows(dat_pathway, min_value = 5000, min_prop_samples = 0.1)
dat_pathway_filter1c <- filter_rows(dat_pathway, min_value = 10000, min_prop_samples = 0.1)

fit_path_con_filter1 <- run_maaslin(dat_pathway_filter_1, meta, paste0("MaAsLin3_path_filter_confoudners1_", Sys.Date()),
                                    "~LIBRA.score*Sex+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_path_con_filter1a <- run_maaslin(dat_pathway_filter1a, meta, paste0("MaAsLin3_path_filter_confoudners1a_", Sys.Date()),
                                   "~LIBRA.score*Sex+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)
fit_path_con_filter1b <- run_maaslin(dat_pathway_filter1b, meta, paste0("MaAsLin3_path_filter_confoudners1b_", Sys.Date()),
                                   "~LIBRA.score*Sex+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)
fit_path_con_filter1c <- run_maaslin(dat_pathway_filter1c, meta, paste0("MaAsLin3_path_filter_confoudners1c_", Sys.Date()),
                                   "~LIBRA.score*Sex+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)






fit_path <- run_maaslin(dat_pathway, meta, paste0("MaAsLin3_path_no_confoudners_", Sys.Date()),
                                  "~LIBRA.score", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_path_filter <- run_maaslin(dat_pathway_filter, meta, paste0("MaAsLin3_path_filter_no_confoudners_", Sys.Date()),
                                  "~LIBRA.score", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_path_filter2 <- run_maaslin(dat_pathway_filter2, meta, paste0( "MaAsLin3_path_filter_no_confoudners2_", Sys.Date()),
                                   "~LIBRA.score", min_prevalence = 0, min_abundance = 0, max_significance=0.2)

fit_path_filter3 <- run_maaslin(dat_pathway_filter3, meta, paste0( "MaAsLin3_path_filter_no_confoudners3_", Sys.Date()),
                                   "~LIBRA.score", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_path_filter4 <- run_maaslin(dat_pathway_filter4, meta, paste0("MaAsLin3_path_filter_no_confoudners4_", Sys.Date()),
                                   "~LIBRA.score", min_prevalence = 0, min_abundance = 0,max_significance=0.2)


fit_path_filter1a <- run_maaslin(dat_pathway_filter1a, meta, paste0("MaAsLin3_path_filter_no_confoudners1a_", Sys.Date()),
                                "~LIBRA.score", min_prevalence = 0, min_abundance = 0,max_significance=0.2)


fit_con_path <- run_maaslin(dat_pathway, meta, paste0("MaAsLin3_path_confoudners_", Sys.Date()),
                        "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_path_con_filter <- run_maaslin(dat_pathway_filter, meta, paste0("MaAsLin3_path_filter_confoudners_", Sys.Date()),
                               "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_path_con_filter2 <- run_maaslin(dat_pathway_filter2, meta, paste0( "MaAsLin3_path_filter_confoudners2_", Sys.Date()),
                                "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0, max_significance=0.2)

fit_path_con_filter3 <- run_maaslin(dat_pathway_filter3, meta, paste0( "MaAsLin3_path_filter_confoudners3_", Sys.Date()),
                                "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_path_con_filter4 <- run_maaslin(dat_pathway_filter4, meta, paste0("MaAsLin3_path_filter_confoudners4_", Sys.Date()),
                                "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_path_con_filter1 <- run_maaslin(dat_pathway_filter_1, meta, paste0("MaAsLin3_path_filter_confoudners1_", Sys.Date()),
                                    "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)



library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(DescTools)

log2(1.0354)
log2(1.032703)
log2(1.5)
0.05*12.6
0.5849625/12.6
2^0.0464256
0.0464254*12.6
2^0.58496
2^0.55
  
volcplot_maaslin3 <- function(fit,
                              variable,                      # metadata variable to plot
                              level = NULL,                  # for categorical variables: the `value` level
                              padj_threshold = 0.2,
                              fc = 1.032703,
                              qval_col = "qval_individual",  # or "qval_joint"
                              coef_is_natural_log = F,    # convert coef to log2 scale
                              plot_title = "Volcano Plot",
                              plot_subtitle = NULL,
                              feature_vector = NULL,         # optional: only label these features
                              feature_filter = FALSE,        # TRUE = drop all features not in feature_vector
                              label_features = TRUE,
                              n_labels = 6) {
  
  # Accept the full maaslin3() object, fit$fit_data_abundance, or the results data frame itself
  if (is.data.frame(fit)) {
    res <- fit
  } else if (!is.null(fit$fit_data_abundance)) {
    res <- fit$fit_data_abundance$results
  } else if (!is.null(fit$results)) {
    res <- fit$results
  } else {
    stop("Couldn't find abundance results. Pass `fit`, `fit$fit_data_abundance`, or the results data frame.")
  }
  
  neg_log2fc <- -log2(fc)
  pos_log2fc <- log2(fc)
  
  plot_ready_data <- res %>%
    dplyr::filter(name == variable) %>%
    { if (!is.null(level)) dplyr::filter(., value == level) else . } %>%
    mutate(
      logFC = if (coef_is_natural_log) coef / log(2) else coef,
      adj.P.Val = .data[[qval_col]]
    ) %>%
    mutate(
      adj.P.Val = replace_na(adj.P.Val, 1),
      logFC = replace_na(logFC, 0),
      log2fc_threshold = case_when(
        logFC >= pos_log2fc & adj.P.Val <= padj_threshold ~ "up",
        logFC <= neg_log2fc & adj.P.Val <= padj_threshold ~ "down",
        TRUE ~ "ns"
      ),
      label = replace_na(feature, "none")
    )
  
  if (nrow(plot_ready_data) == 0) {
    stop("No rows found for variable '", variable, "'",
         if (!is.null(level)) paste0(" and level '", level, "'"), ".")
  }
  
  if (feature_filter && !is.null(feature_vector)) {
    plot_ready_data <- plot_ready_data %>% dplyr::filter(feature %in% feature_vector)
  }
  
  # Only label significant features (restricted to feature_vector if provided)
  # Label only the n_labels features with the lowest q-values
  # (restricted to feature_vector if provided)
  eligible <- plot_ready_data$log2fc_threshold != "ns"
  if (!is.null(feature_vector)) {
    eligible <- eligible & plot_ready_data$feature %in% feature_vector
  }
  q_rank <- rank(ifelse(eligible, plot_ready_data$adj.P.Val, Inf), ties.method = "first")
  
  plot_ready_data <- plot_ready_data %>%
    mutate(label = ifelse(label_features & eligible & q_rank <= n_labels, feature, ""))  
  x_axis_limits <- DescTools::RoundTo(max(abs(plot_ready_data$logFC)), 0.1, ceiling)
  y_max <- max(1.5, ceiling(max(-log10(plot_ready_data$adj.P.Val)) * 10) / 10)
  #plot_colors <- c("up" = "#1B3558", "ns" = "gray", "down" = "#A77903")
  # Set the plot colors
  plot_colors <- c(
    'up' = 'darkred',
    'ns' = 'gray',
    'down' = 'dodgerblue1'
  )
  
  plot <- ggplot(plot_ready_data) +
    geom_point(alpha = 0.25, size = 1.5) +
    aes(x = logFC, y = -log10(adj.P.Val), color = log2fc_threshold, label = label) +
    geom_vline(xintercept = c(neg_log2fc, pos_log2fc), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed") +
    scale_x_continuous("log2(FC)", limits = c(-x_axis_limits, x_axis_limits)) +
    scale_y_continuous("-log10(q)", limits = c(-0.25, y_max)) +
    scale_color_manual(values = plot_colors) +
    labs(title = plot_title, subtitle = plot_subtitle) +
    theme_bw(base_size = 24) +
    theme(text = element_text(size = 14, face = "bold"), legend.position = "none",
          aspect.ratio = 1, axis.text = element_text(color = "black"))
  
  if (label_features) {
    plot <- plot + geom_label_repel(size = 2.5, force = 1, max.overlaps = 100000,
                                    nudge_x = 0, nudge_y = 0, segment.color = "black",
                                    min.segment.length = 0, show.legend = FALSE)
  }
  plot
}


q_cut <- 0.2
abund <- fit_path_con_filter1$fit_data_abundance$results
# Number of unique features with qval < cutoff, by model and metadata variable
abund %>%
  dplyr::filter(!is.na(qval_individual), qval_individual < q_cut) %>%
  group_by(metadata, value) %>%
  summarise(n_features = n_distinct(feature), .groups = "drop")

volc = volcplot_maaslin3(fit_path_con_filter1$fit_data_abundance, variable = "LIBRA.score:SexFemale",
                         plot_title = "",n_labels = 3)
ggsave(paste0("Figures/Figure_4_LIBRA*SEX_Pathway_Volcano_Plot_",Sys.Date(),".jpeg"), volc, jpeg, height = 8, width = 12, units = "in" )


volc2 = volcplot_maaslin3(fit_path_con_filter1$fit_data_abundance, variable = "LIBRA.score",
                         plot_title = "",n_labels = 10)
ggsave(paste0("Figures/Figure_4_LIBRA_ONLY_Pathway_Volcano_Plot_",Sys.Date(),".jpeg"), volc2, jpeg, height = 8, width = 12, units = "in" )

q_cut <- 0.2
abund <- fit_path_con_filter1a$fit_data_abundance$results
# Number of unique features with qval < cutoff, by model and metadata variable
abund %>%
  dplyr::filter(!is.na(qval_individual), qval_individual < q_cut) %>%
  group_by(metadata, value) %>%
  summarise(n_features = n_distinct(feature), .groups = "drop")

volc2 = volcplot_maaslin3(fit_path_con_filter1a$fit_data_abundance, variable = "LIBRA.score",
                         plot_title = "",n_labels = 3)



# ---- 1. Prep the data ----
fmt_p <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "<0.001",
    p < 0.01  ~ formatC(p, format = "f", digits = 3),
    TRUE      ~ formatC(p, format = "f", digits = 2)
  )
}

q_cut <- 0.2
abund_table <- fit_path_con_filter1$fit_data_abundance$results
#prev_tabel = fit_species_filter1d$fit_data_prevalence$results
full_table = abund_table%>% mutate(model = case_when(model == "linear" ~"Abundnace",
                                                                                   T~"Prevalence")) %>%
  dplyr::filter(feature != "s__Other")

full_table_update <- full_table |>
  mutate(
    ci_lwr = coef - qt(.975, df = N_not_zero - 6) * stderr,
    ci_upr = coef + qt(.975, df = N_not_zero - 6) * stderr
  )


tbl_dat <- full_table_update %>%
  dplyr::select(feature, name, coef, stderr, pval_individual, qval_individual,model,ci_lwr,ci_upr) %>%
  dplyr::filter(!is.na(qval_individual), qval_individual < q_cut & name=="LIBRA.score") %>%
  arrange(desc(coef)) %>%
  mutate(
    coef_ci = sprintf("%.3f [%.3f, %.3f]", coef, ci_lwr,ci_upr),
    # tidy taxon labels (e.g., "s__Genus_species" -> "Genus species")
    feature = feature %>%
      gsub("^s__", "", .) %>%
      gsub("[_.]", " ", .),
    coef  = formatC(coef,   format = "f", digits = 3),
    stderr = formatC(stderr, format = "f", digits = 3),
    pval_individual = fmt_p(pval_individual),
    qval_individual = fmt_p(qval_individual)
  )%>%
  dplyr::select(feature, name, coef_ci,stderr, pval_individual, qval_individual,model)

# ---- 2. Build the flextable ----
sup_table_mad <- flextable(tbl_dat) %>%
  set_header_labels(
    feature         = "Pathway",
    name            = "Term",
    coef_ci            = "Coefficient [95% CI]",
    stderr          = "Standard Error",
    pval_individual = "p-value",
    qval_individual = "q-value (FDR)",
    model = "Model"
  ) %>%
  set_caption(
    caption = "Supplementary Table X. Pathways significantly associated with LIBRA Score."
  ) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  italic(j = "feature", part = "body") %>%       # italicize species names
  align(j = c("coef_ci", "stderr", "pval_individual", "qval_individual"),
        align = "center", part = "all") %>%
  align(j = c("feature", "name","model"), align = "left", part = "all") %>%
  # three-line (journal-style) borders
  border_remove() %>%
  hline_top(border = fp_border(width = 1.5), part = "all") %>%
  hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  hline(border = fp_border(width = 0.75), part = "header") %>%
  autofit()

# ---- 3. Save ----
save_as_docx(
  sup_table_mad,
  path = paste0("Supplementary_Table_LIBRA_Only_pathways_", Sys.Date(), ".docx")
)



####LIBRA x SEx
tbl_dat2 <- full_table_update %>%
  dplyr::select(feature, name, coef, stderr, pval_individual, qval_individual,model,ci_lwr,ci_upr) %>%
  dplyr::filter(!is.na(qval_individual), qval_individual < q_cut & name=="LIBRA.score:SexFemale") %>%
  arrange(desc(coef)) %>%
  mutate(
    coef_ci = sprintf("%.3f [%.3f, %.3f]", coef, ci_lwr,ci_upr),
    # tidy taxon labels (e.g., "s__Genus_species" -> "Genus species")
    feature = feature %>%
      gsub("^s__", "", .) %>%
      gsub("[_.]", " ", .),
    coef  = formatC(coef,   format = "f", digits = 3),
    stderr = formatC(stderr, format = "f", digits = 3),
    pval_individual = fmt_p(pval_individual),
    qval_individual = fmt_p(qval_individual)
  )%>%
  dplyr::select(feature, name, coef_ci,stderr, pval_individual, qval_individual,model)

# ---- 2. Build the flextable ----
sup_table_mad2 <- flextable(tbl_dat2) %>%
  set_header_labels(
    feature         = "Pathways",
    name            = "Term",
    coef_ci            = "Coefficient",
    stderr          = "Standard Error",
    pval_individual = "p-value",
    qval_individual = "q-value (FDR)",
    model = "Model"
  ) %>%
  set_caption(
    caption = "Supplementary Table X. Pathways significantly associated with LIBRA Score x Sex Interaction."
  ) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  italic(j = "feature", part = "body") %>%       # italicize species names
  align(j = c("coef_ci", "stderr", "pval_individual", "qval_individual"),
        align = "center", part = "all") %>%
  align(j = c("feature", "name","model"), align = "left", part = "all") %>%
  # three-line (journal-style) borders
  border_remove() %>%
  hline_top(border = fp_border(width = 1.5), part = "all") %>%
  hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  hline(border = fp_border(width = 0.75), part = "header") %>%
  autofit()

# ---- 3. Save ----
save_as_docx(
  sup_table_mad2,
  path = paste0("Supplementary_Table_LIBRA*SEX_pathways_", Sys.Date(), ".docx")
)







-log2(1.0354)
.05*12
ggsave(paste0("Figures/Figure_4_Pathway_Volcano_Plot_",Sys.Date(),".jpeg"), volc, jpeg, height = 8, width = 12, units = "in" )




genes = sigtab2 %>% dplyr::filter(log2FoldChange > log2(1.0354) | log2FoldChange< -log2(1.0354) ) %>%
  arrange(desc(abs(log2FoldChange))) %>% slice_head(n=6)




gen.fit_LIBRA <- Maaslin2(input_data = dat_pathway, input_metadata = meta, min_prevalence = 0, 
                          min_abundance = 0, normalization = "TSS", 
                          output = "MaAsLin2_Out_Pathways_final-10-5-26", 
                          fixed_effects = c("LIBRA.score","Age","BMI"),
                          plot_scatter = F,
                          plot_heatmap = F,
                          correction="BH",
                          analysis_method = "LM",
                          transform = "LOG",
                          max_significance = 0.25)


path.fit_LIBRA = maaslin3(
  input_data = dat_pathway,
  input_metadata = meta,
  #output = paste0(new_dirc, "/", name, "/", name, "_MaAsLin3_output"),
  output= paste0("MaAsLin3_full_path_model_no_confoudners_",Sys.Date()),
  #formula = "~group_cat+age+group(race_cat)+group(Picogreen_cat)+BMI",
  formula = "~LIBRA.score",
  min_prevalence = 0,
  min_abundance = 0,
  plot_associations = F,
  max_significance = 0.2,
  median_comparison_prevalence = FALSE)


path.fit_LIBRA = maaslin3(
  input_data = dat_pathway2,
  input_metadata = meta,
  #output = paste0(new_dirc, "/", name, "/", name, "_MaAsLin3_output"),
  output= paste0("MaAsLin3_full_path2_model_",Sys.Date()),
  #formula = "~group_cat+age+group(race_cat)+group(Picogreen_cat)+BMI",
  formula = "~LIBRA.score + Age + group(Homeless.ever)",
  min_prevalence = 0,
  min_abundance = 0,
  plot_associations = F,
  max_significance = 0.2)


#align data
ind <- match(x = colnames(dat_pathway), table = rownames(meta))
meta1 <- meta[ind,] 
meta1[1:10, 1:10]
#Check that data are aligned
all(colnames(dat_pathway) == rownames(meta1))
#make deseq object
meta1
library(caret)
#scale and center data
meta1_scale_factors = subset(meta1, select = c(Age, BMI, LIBRA.score))
meta1_scaled  = scale(meta1_scale_factors, scale = T, center = T)

#Create deseq object
dds <- DESeqDataSetFromMatrix(countData = dat_pathway, colData = meta1_scaled, design = ~ Age + BMI + LIBRA.score)
#generate normalized counts 
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

#QC for DE analysis using DESeq2
#If input file is very large > 800 MB then, use the below command
rld <- varianceStabilizingTransformation(dds, blind=T)

#PCA 
plotPCA(rld, intgroup="LIBRA.score")
#Hierarchical Clustering
### Extract the rlog matrix from the object
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(rld_cor)
#differential expression analysis
dds <- DESeq(dds)
plotDispEsts(dds)
res <- results(dds)
#shrinkage factor to normalize dataa
dds_shrink = lfcShrink(dds = dds, res = res, coef="LIBRA.score")
summary(res)
summary(dds_shrink)
summary(dds)

plotMA(res, ylim=c(-2,2))
# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(dds_shrink))

# Extract significant results
signif_res <- dds_shrink[dds_shrink$padj < 0.3 & dds_shrink$pvalue < 0.05 & !is.na(dds_shrink$padj), ]

plotMA(signif_res, ylim=c(-4,4))
#This will be empty if no genes were significant
signif_genes <- as.character(rownames(signif_res))


#Below is old code to export sig-taxa
sigtab = cbind(as(dds_shrink, "data.frame"), all_genes)


##Function to create volcano plot
volcplot <- function(data, padj_threshold = 0.05, fc = 1.0354, plot_title = 'Volcano Plot', plot_subtitle = NULL, genelist_vector = genes$all_genes, genelist_filter = FALSE) {
  
  # Set the fold-change thresholds
  neg_log2fc <- -log2(fc)
  pos_log2fc <- log2(fc)
  
  # Make a dataset for plotting, add the status as a new column
  plot_ready_data <- data %>%
    mutate_at('padj', ~replace(.x, is.na(.x), 1)) %>%
    mutate_at('log2FoldChange', ~replace(.x, is.na(.x), 0)) %>%
    mutate(
      log2fc_threshold = ifelse(log2FoldChange >= pos_log2fc & pvalue <= padj_threshold, 'up',
                                ifelse(log2FoldChange <= neg_log2fc & pvalue <= padj_threshold, 'down', 'ns')
      )
    ) %>%
    mutate(all_genes = replace_na(all_genes, 'none'))
  
  if (genelist_filter) {
    plot_ready_data <- plot_ready_data %>% dplyr::filter(all_genes %in% genelist_vector)
  }
  
  if(!is.null(genelist_vector)) {
    plot_ready_data <- plot_ready_data %>% mutate(all_genes = ifelse(all_genes %in% genelist_vector & pvalue < padj_threshold & log2fc_threshold != 'ns', all_genes, ''))
  }
  
  # Get the number of up, down, and unchanged genes
  up_genes <- plot_ready_data %>% dplyr::filter(log2fc_threshold == 'up') %>% nrow()
  down_genes <- plot_ready_data %>% dplyr::filter(log2fc_threshold == 'down') %>% nrow()
  unchanged_genes <- plot_ready_data %>% dplyr::filter(log2fc_threshold == 'ns') %>% nrow()
  
  # Make the labels for the legend
  legend_labels <- c(
    str_c('Up: ', up_genes),
    str_c('NS: ', unchanged_genes),
    str_c('Down: ', down_genes)
  )
  
  # Set the x axis limits, rounded to the next even number
  x_axis_limits <- DescTools::RoundTo(
    max(abs(plot_ready_data$log2FoldChange)),
    0.5,
    ceiling
  )
  
  # Set the plot colors
  plot_colors <- c(
    'up' = 'darkred',
    'ns' = 'gray',
    'down' = 'dodgerblue1'
  )
  
  
  # Make the plot, these options are a reasonable strting point
  plot <- ggplot(plot_ready_data) +
    geom_point(
      alpha = 0.25,
      size = 1.5
    ) +
    aes(
      x = log2FoldChange,
      y = -log10(pvalue),
      color = log2fc_threshold,
      label = all_genes
    ) +
    geom_vline(
      xintercept = c(neg_log2fc, pos_log2fc),
      linetype = 'dashed'
    ) +
    geom_hline(
      yintercept = -log10(padj_threshold),
      linetype = 'dashed'
    ) +
    scale_x_continuous(
      'log2(FC)',
      limits = c(-x_axis_limits, x_axis_limits)
    ) +
    scale_color_manual(
      values = plot_colors#,
      #labels = legend_labels
    ) +
    #labs(
    #  color = str_c(fc, '-fold, padj ≤', padj_threshold),
    #  title = plot_title,
    #  subtitle = plot_subtitle
    #) +
    theme_bw(base_size = 24) +
    theme(legend.position = "none",
          aspect.ratio = 1,
          axis.text = element_text(color = 'black'),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),  # Reduces dead area around legend
          legend.spacing.x = unit(0.2, 'cm'),
          legend.title = element_blank()
    )
  
  # Add gene labels if needed
  if (!is.null(genelist_vector)) {
    plot <- plot +
      geom_label_repel(
        size = 2.5,
        force = 1,
        max.overlaps = 100000,
        nudge_x = -.05,
        nudge_y = -.20,
        segment.color = 'black',
        min.segment.length = 0,
        show.legend = FALSE
      )
  }
  plot
}


sigtab2 = na.omit(sigtab)
rownames(sigtab2) <- NULL

#pulls out genes for labeling
genes = sigtab2 %>% dplyr::filter(log2FoldChange > log2(1.0354) | log2FoldChange< -log2(1.0354) ) %>%
  arrange(desc(abs(log2FoldChange)))# %>% slice_head(n=6)

out_table = sigtab2 %>% dplyr::filter(padj < 0.3 & pvalue < 0.05) %>%
  arrange(desc(log2FoldChange)) 

#creates volcano plot and prints
library(ggrepel)
volc = volcplot(sigtab2)

ggsave(paste0("Figures/Figure_4_Pathway_Volcano_Plot_",Sys.Date(),".jpeg"), volc, jpeg, height = 8, width = 12, units = "in" )


#Table of signifincat pathways
out_table_pathways = out_table %>% dplyr::select(all_genes, log2FoldChange, lfcSE, pvalue, padj) %>% arrange(desc(log2FoldChange))

sup_table_3 = flextable(out_table_pathways)
sup_table_3 = set_caption(sup_table_3, caption = "Supplementary Table 3. Metabolic Pathways associated with LIBRA Scores") %>% autofit()
save_as_docx(sup_table_3, path = "Figures/Supplementary_Table_3_DESeq_Results.docx")




######Extra for MaAsLin
#MaAsLin analysis for each libra factor
mal_fun= function(effect,ref){
  fit = Maaslin2(input_data = input_data, input_metadata = meta, min_prevalence = .5, 
                 min_abundance = .01, normalization = "TSS", 
                 output = paste0("MaAsLin_Out_species-level_",effect,"_",Sys.Date()), 
                 fixed_effects = effect,
                 reference = ref,
                 plot_scatter = F,
                 plot_heatmap = F,                          
                 correction="BH",
                 analysis_method = "LM",
                 transform = "LOG",
                 max_significance = 0.3)#,
  #reference = "Current_alcohol,No")}
  return(fit)
}






#EC Analysis
######
ec_table_1 = read.table('/Volumes/B4_Backup/Metagemomic-sequences/MQ3201_0001_DeepSeq_Jul_2023_Diversigen/filtered functional tables/sample_ec4.tsv', header =T,sep="\t", check.names = F, row.names = 1)

ec_table_2 = read.table('/Volumes/B4_Backup/Metagemomic-sequences/MQ3201_0003_DeepSeq_May_2024_Diversigen/filtered functional tables/sample_ec4.tsv', header =T,sep="\t", check.names = F, row.names = 1)
meta
ec_table_1 <- ec_table_1[-c(1, 2), ]
ec_table_2 <- ec_table_2[-c(1, 2), ]
ec_table_all = merge(ec_table_1,ec_table_2, by=0, all=T)
ec_table_all[is.na(ec_table_all)] <- 0
rownames(ec_table_all) = ec_table_all$Row.names
ec_table_all = ec_table_all %>% dplyr::select(-Row.names)
# Remove the first two rows
#taxa_table_all <- taxa_table_all[-c(1, 2), ]
ec_table_all_df = data.frame(ec_table_all)
names(ec_table_all_df) <- gsub("X", "", names(ec_table_all_df))
meta$Sample_ID = rownames(meta)
ec_table_all_df_filter <- ec_table_all_df[, names(ec_table_all_df) %in% meta$Sample_ID]
#meta_and_alpha %>% filter(SampleID == "19-1307.P24.8m-g")
#dat[-1] <-sweep(dat[-1], 2, colSums(dat[,-1]), `/`) * 10000
#Round data so that all the data are integers

ec_df = data.frame(ec_table_all_df_filter)
dat_ec <- ec_df %>%
  dplyr::mutate_if(is.numeric, round)
names(dat_ec) <- gsub("X", "", names(dat_ec))
keep_genes <- colSums( dat_ec > 100 ) >= 70
dat_ec <- dat_ec[keep_genes,]
dat_ec_2 <- as.matrix(dat_ec)


dat_ec_2



trans_dat_ec = t(data.frame(dat_ec))


#trans_dat_pathway = trans_dat_pathway %>% row_to_names(row_number = 1)
trans_dat_ec = data.frame(trans_dat_ec)
trans_dat_ec$Sample_ID = rownames(trans_dat_ec)


trans_dat_ec_ID = trans_dat_ec %>% mutate(Sample_ID = str_replace(Sample_ID, "X", ""))
rownames(trans_dat_ec_ID) = trans_dat_ec_ID$Sample_ID
trans_ec_update = trans_dat_ec_ID %>% 
  dplyr::select(-Sample_ID) %>% 
  dplyr::select(where(~ any(. != 0)))
names(trans_ec_update) <- gsub("X", "", names(trans_ec_update))



filter_ec_table = read.csv('/Volumes/BrainHealth/B4 Health Lab/Matt_Folder/Projects/2026/Mindsource/Data/ec_number_list_final-09-29-26.csv')
trans_ec_update_filter <- trans_ec_update %>% 
  dplyr::select(any_of(filter_ec_table$EC_Number))

trans_ec_update_filter = t(trans_ec_update_filter)
trans_ec_update2 = t(trans_ec_update)

trans_ec_update_filter1a <- filter_rows(trans_ec_update_filter, min_value = 10, min_prop_samples = 0.1)
trans_ec_update_filter1b <- filter_rows(trans_ec_update_filter, min_value = 100, min_prop_samples = 0.1)
trans_ec_update_filter1c <- filter_rows(trans_ec_update_filter, min_value = 100, min_prop_samples = 0.1)
trans_ec_update_filter1d <- filter_rows(trans_ec_update_filter, min_value = 1000, min_prop_samples = 0.1)


trans_ec_update_filter2 <- filter_rows(trans_ec_update_filter, min_value = 10, min_prop_samples = 0.25)
trans_ec_update_filter3 <- filter_rows(trans_ec_update_filter, min_value = 10, min_prop_samples = 0.5)
trans_ec_update_filter4 <- filter_rows(trans_ec_update_filter, min_value = 100, min_prop_samples = 0.25)

fit_ec_con_filter1a <- run_maaslin(trans_ec_update_filter1a, meta, paste0("MaAsLin3_ec_filter_confoudners1a_", Sys.Date()),
                                 "~LIBRA.score*Sex+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)
fit_ec_con_filter1b <- run_maaslin(trans_ec_update_filter1b, meta, paste0("MaAsLin3_ec_filter_confoudners1b_", Sys.Date()),
                                 "~LIBRA.score*Sex+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)
fit_ec_con_filter1c <- run_maaslin(trans_ec_update_filter1c, meta, paste0("MaAsLin3_ec_filter_confoudners1c_", Sys.Date()),
                                 "~LIBRA.score*Sex+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)
fit_ec_con_filter1d <- run_maaslin(trans_ec_update_filter1d, meta, paste0("MaAsLin3_ec_filter_confoudners1d_", Sys.Date()),
                                 "~LIBRA.score*Sex+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)





fit_ec_con_filter <- run_maaslin(trans_ec_update_filter1, meta, paste0("MaAsLin3_ec_filter_confoudners_", Sys.Date()),
                                   "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_ec_con_filter2 <- run_maaslin(trans_ec_update_filter2, meta, paste0( "MaAsLin3_ec_filter_confoudners2_", Sys.Date()),
                                    "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0, max_significance=0.2)

fit_ec_con_filter3 <- run_maaslin(trans_ec_update_filter3, meta, paste0( "MaAsLin3_ec_filter_confoudners3_", Sys.Date()),
                                    "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_ec_con_filter4 <- run_maaslin(trans_ec_update_filter4, meta, paste0("MaAsLin3_ec_filter_confoudners4_", Sys.Date()),
                                    "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)




trans_ec_update_all_1 <- filter_rows(trans_ec_update2, min_value = 10, min_prop_samples = 0.1)
trans_ec_update_all_2 <- filter_rows(trans_ec_update2, min_value = 10, min_prop_samples = 0.25)
trans_ec_update_all_3 <- filter_rows(trans_ec_update2, min_value = 10, min_prop_samples = 0.5)
trans_ec_update_all_4 <- filter_rows(trans_ec_update2, min_value = 100, min_prop_samples = 0.25)
trans_ec_update_all_1a <- filter_rows(trans_ec_update2, min_value = 10, min_prop_samples = 0.1)
trans_ec_update_all_1b <- filter_rows(trans_ec_update2, min_value = 100, min_prop_samples = 0.1)
trans_ec_update_all_1c <- filter_rows(trans_ec_update2, min_value = 1000, min_prop_samples = 0.1)
trans_ec_update_all_1d <- filter_rows(trans_ec_update2, min_value = 10000, min_prop_samples = 0.1)


fit_ec_con1a <- run_maaslin(trans_ec_update_all_1a, meta, paste0("MaAsLin3_ec_confoudners1a_", Sys.Date()),
                           "~LIBRA.score*Sex+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_ec_con1b <- run_maaslin(trans_ec_update_all_1b, meta, paste0("MaAsLin3_ec_confoudners1b_", Sys.Date()),
                           "~LIBRA.score*Sex+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_ec_con1c <- run_maaslin(trans_ec_update_all_1c, meta, paste0("MaAsLin3_ec_confoudners1c_", Sys.Date()),
                           "~LIBRA.score*Sex+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_ec_con1d <- run_maaslin(trans_ec_update_all_1d, meta, paste0("MaAsLin3_ec_confoudners1d_", Sys.Date()),
                           "~LIBRA.score*Sex+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)



fit_ec_con_no_filter <- run_maaslin(trans_ec_update2, meta, paste0("MaAsLin3_ec_filter_all_confoudners_", Sys.Date()),
                                 "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)


fit_ec_con1 <- run_maaslin(trans_ec_update_all_1, meta, paste0("MaAsLin3_ec_confoudners_", Sys.Date()),
                                 "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_ec_con2 <- run_maaslin(trans_ec_update_all_2, meta, paste0( "MaAsLin3_ec_confoudners2_", Sys.Date()),
                                  "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0, max_significance=0.2)

fit_ec_con3 <- run_maaslin(trans_ec_update_all_3, meta, paste0( "MaAsLin3_ec_confoudners3_", Sys.Date()),
                                  "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_ec_con4 <- run_maaslin(trans_ec_update_all_4, meta, paste0("MaAsLin3_ec_confoudners4_", Sys.Date()),
                                  "~LIBRA.score+ Age + group(Homeless.ever) + BMI", min_prevalence = 0, min_abundance = 0,max_significance=0.2)




gen.fit_LIBRA <- Maaslin2(input_data = trans_ec_update_filter, 
                          input_metadata = meta, min_prevalence = 0, 
                          min_abundance = 0, normalization = "TSS", 
                          output = "MaAsLin2_Out_EC_final-10-5-26", 
                          fixed_effects = c("LIBRA.score","Age","BMI"),
                          plot_scatter = F,
                          plot_heatmap = F,
                          correction="BH",
                          analysis_method = "LM",
                          transform = "LOG",
                          max_significance = 0.25)

path.fit_LIBRA = maaslin3(
  input_data = trans_ec_update_filter,
  input_metadata = meta,
  #output = paste0(new_dirc, "/", name, "/", name, "_MaAsLin3_output"),
  output= paste0("MaAsLin3_full_ec_model_",Sys.Date()),
  #formula = "~group_cat+age+group(race_cat)+group(Picogreen_cat)+BMI",
  formula = "~LIBRA.score + Age + group(Homeless.ever)",
  min_prevalence = 0,
  min_abundance = 0,
  plot_associations = F,
  max_significance = 0.2)

path.fit_LIBRA = maaslin3(
  input_data = trans_ec_update_filter,
  input_metadata = meta,
  #output = paste0(new_dirc, "/", name, "/", name, "_MaAsLin3_output"),
  output= paste0("MaAsLin3_no_confounders_ec_model_",Sys.Date()),
  #formula = "~group_cat+age+group(race_cat)+group(Picogreen_cat)+BMI",
  formula = "~LIBRA.score ",
  min_prevalence = 0,
  min_abundance = 0,
  plot_associations = F,
  max_significance = 0.2)




trans_ec_update_filter2 <- trans_ec_update_filter %>%
  rownames_to_column(var = "Sample_ID")
trans_ec_update_filter2_melt = reshape2::melt(trans_ec_update_filter2, id.var="Sample_ID", variable.name="EC", value.name = "Abundance")
trans_ec_update_filter2_melt2 = trans_ec_update_filter2_melt %>%
  mutate(Sample_ID = str_replace(Sample_ID, "X", "")) %>%
  dplyr::filter(Sample_ID %in% meta_and_libra_final$Sample_ID)

trans_ec_update_filter2_melt3 = trans_ec_update_filter2_melt2 %>%
  group_by(Sample_ID, EC) %>%
  mutate(total_abd = sum(as.numeric(Abundance))) %>%
  ungroup() %>%
  select(-Abundance) %>%
  distinct()
all_taxa_meta_data_ec = merge(trans_ec_update_filter2_melt3, meta_and_libra_final, by="Sample_ID")
all_taxa_meta_data_ec_update = all_taxa_meta_data_ec %>% mutate(Prev = case_when(total_abd >0 ~ 1,
                                                                                 T~0))




ec_abund_bin = ggplot(all_taxa_meta_data_ec_update, aes(x=reorder(Participant, Total_LIBRA_Score), y=total_abd, color=LIBRA_Half))+
  geom_boxplot()+
  geom_point()+
  labs(tag = "B", x="Participants \n(Ordered by LIBRA Score)", y="Summed Relative Abundance of \nSCFA Associated Bacteria", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~EC)+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )




ec_fig_supp = ggplot(all_taxa_meta_data_ec_update %>% filter(EC == "2.7.2.7"), aes(y=total_abd,x=Time_num,group=LIBRA_Half,color=LIBRA_Half))+
  #geom_line(aes(group=Participant,color=LIBRA_Half))+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F)+
  labs( x="Time (Months)", y="Relative Abundance", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  facet_wrap(~Participant, scales="free") +
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )+
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)


fit_data_ec2 = fit_data_ec
filter_ec_table$Pathway...Reaction
fit_data_ec2$fit_data_abundance$results$feature <- filter_ec_table$Enzyme[match(fit_data_ec2$fit_data_abundance$results$feature, filter_ec_table$EC_Number)]
fit_data_ec2$fit_data_abundance$results$feature <- paste0(fit_data_ec2$fit_data_abundance$results$feature, ": ", filter_ec_table$EC_Number[match(fit_data_ec2$fit_data_abundance$results$feature, filter_ec_table$Enzyme)])


de_abund_full_libra = de_graph_fun(fit_data_ec2$fit_data_abundance$results, "Differential Expression of EC Numbers", "Time_num:Total_LIBRA_Score")

de_abund_full_libra2 = de_abund_full_libra +
  labs(x="Enzyme Commission Number")







meta$Physical.inactivity_lib.weight = as.factor(meta$Physical.inactivity_lib.weight)
meta$Healthy_diet
depressed_mal = mal_fun("Depressed", "Depressed,No")

obese_mal = mal_fun("Obese", "Obese,No")

alc_mal = mal_fun("Current_alcohol","Current_alcohol,No")

chd_mal = mal_fun("Coronary.heart.disease", "Coronary.heart.disease,No")

cdd_mal = mal_fun("Chronic_kidney_disease", "Chronic_kidney_disease,No")

diabetes_mal = mal_fun("Diabetes","Diabetes,No")

Hypercholesterolemia_mal = mal_fun("Hypercholesterolemia", "Hypercholesterolemia,No")

Lifetime_smoking_mal = mal_fun("Lifetime_smoking", "Lifetime_smoking,No")

Hypertension_mal = mal_fun("Hypertension","Hypertension,No")

Physical.inactivity_lib.weight_mal = mal_fun("Physical.inactivity_lib.weight", "Physical.inactivity_lib.weight,0")

Healthy_diet_mal = mal_fun("Healthy_diet", "Healthy_diet,Yes")




#####Supplemental figure 2
####################################################################
#CCA for pathways to make biplot
#import distance matrix
x_pathway 

####
#*****************************IMPORTANT******************************************
#next line will determine which metadata is being piped in so make sure to check before proceeding 
#meta.in <- read.csv("LIBRA/LIBRA_variables_CCA_noNAs_021225.csv", stringsAsFactors = T)
#matching up metadata with the distance matrix

ind <- match(x = rownames(x_pathway), table = meta_env_fit_update$SampleID)
ind1 <- !is.na(ind)
x1 <- x_pathway[ind1,ind1]
ind2 <- match(x = rownames(x1), table = meta_env_fit_update$SampleID)
meta.in1 <- meta_env_fit_update[ind2,]
#this next line should be all true - if not, then something is wrong
table(rownames(x1) == meta.in1$SampleID) 
rownames(x1)
colnames(x1)
meta.in1$SampleID
#read in pathways
pathway_abundance <- read.csv(file = '/Volumes/B4_Backup/Stamper_Mac_VA_files_060325/NFL-grant/Metagenomic_data/filtered functional tables/uniref-pathway-filtered-absolute_NFL1.csv', 
                              row.names = 1, header = T, check.names = F)

Taxa <- t(pathway_abundance)
Taxa[1:10,1:10]

ind1 <- match(x = colnames(x1), table = rownames(Taxa))
Taxa1 <- Taxa[ind1,]
Taxa1[1:10,1:10]

#Uses only pathways that are significant for DESeq
colnames(Taxa1)
ind <- match(colnames(Taxa1), rownames(signif_res))
ind
ind1 <- !is.na(ind)
Taxa2 <- Taxa1[,ind1]
colnames(Taxa2)


#**Special plot** CCA start
#this code is for coloring dots on biplot by quartile
meta.in$LIBRA.score.quartile <- as.factor(meta.in$LIBRA.score.quartile)
with(meta.in, levels(LIBRA.score.quartile))
scl <- 2

#CCA start
pca_w=capscale(x1~1, data=meta.in1[,-1])
fit1=envfit(pca_w, meta.in1[,-1], na.rm = T)
#p.max in the code below will determine which variables end up on the graph - adjust for more or less
#p.max in the code below will determine which taxa end up on the graph - adjust for more or less
fit2=envfit(pca_w, Taxa2, na.rm = T)



pca_data_mds = data.frame(vegan::scores(pca_w)$sites)
#pca_data_mds = pca_data %>% select(MDS1,MDS2)
pca_data_mds$SampleID = rownames(pca_data_mds)

libra_fit = data.frame(vegan::scores(fit1, "vectors")) *ordiArrowMul(fit1)
libra_pca_pval = data.frame(fit1$vectors$pvals) %>% dplyr::filter(fit1.vectors.pvals <= 0.5)
libra_fit_filter <- libra_fit %>%
  dplyr::filter(row.names(.) %in% rownames(libra_pca_pval)) %>%
  mutate(indicator= "libra")


phylum_fit = data.frame(vegan::scores(fit2, "vectors")) *ordiArrowMul(fit2)
phylum_pca_pval = data.frame(fit2$vectors$pvals) %>% dplyr::filter(fit2.vectors.pvals <= 0.001)
phylum_fit_filter <- phylum_fit %>%
  dplyr::filter(row.names(.) %in% rownames(phylum_pca_pval)) %>%
  mutate(indicator= "species")

all_fit_lines = bind_rows(libra_fit_filter,phylum_fit_filter)
pca_data_full = full_join(meta_env_fit_update,pca_data_mds, by="SampleID" )
pca_data_full$LIBRA.score.quartile = factor(pca_data_full$LIBRA.score.quartile, levels=c("4","3","2","1"))

colvec <- rev(c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A"))

pathway_biplot = ggplot(pca_data_full, aes(x=MDS1, y=MDS2))+
  geom_point(aes(color= LIBRA.score.quartile))+
  scale_colour_manual(values = colvec)+
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), 
               data = libra_fit_filter, size =1, alpha = 0.5, colour = "skyblue")+
  ggrepel::geom_text_repel(max.overlaps=nrow(libra_fit_filter),data = libra_fit_filter, aes(x = MDS1, y = MDS2), 
                           label = row.names(libra_fit_filter), colour = "skyblue",nudge_x=-0.01,nudge_y=0.01)+
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), 
               data = phylum_fit_filter, size =1, alpha = 0.5, colour = "darkred")+
  ggrepel::geom_text_repel(max.overlaps=15,data = phylum_fit_filter, aes(x = MDS1, y = MDS2), 
                           label = row.names(phylum_fit_filter), colour = "darkred",nudge_y=0.01,nudge_x=0.01) 

ggsave(paste0("Figures/Figure_supp_2_Biplot_Pathway_",Sys.Date(),".jpeg"), pathway_biplot, jpeg, height = 8, width = 12, units = "in" )




###Extra
########Figure 2#######
####relative abundance of significant species######
taxa_species

mas.sig 

mas.sig_filter = mas.sig %>% dplyr::filter(qval < 0.1)


relative_taxa_adbund_filtered = taxa_species_relative %>% dplyr::filter(row.names(.) %in% mas.sig$feature)

relative_taxa_adbund_filtered_update <- relative_taxa_adbund_filtered %>%
  rownames_to_column(var = "Species")

relative_taxa_adbund_filtered_update2 <- data.frame(t(relative_taxa_adbund_filtered_update[-1]))
colnames(relative_taxa_adbund_filtered_update2) <- relative_taxa_adbund_filtered_update[, 1]

relative_taxa_adbund_filtered_update3 <- relative_taxa_adbund_filtered_update2 %>%
  rownames_to_column(var = "SampleID")

relative_taxa_adbund_filtered_update3_melt = melt(relative_taxa_adbund_filtered_update3, id.var="SampleID", variable.name="Species", value.name = "Abundance") 


all_taxa_meta_data = merge(relative_taxa_adbund_filtered_update3_melt, meta.in, by="SampleID")

abundance_plot = ggplot(all_taxa_meta_data, aes(y=Abundance, x=LIBRA.score)) + 
  geom_point() +
  facet_wrap(~Species, scales="free") +
  geom_smooth(method="lm") +
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)
ggsave(paste0("Figures/Figure_extra_2_Abundance_Change_Libra_Significant_Species_",Sys.Date(),".jpeg"), abundance_plot, jpeg, height = 20, width = 25, units = "in" )
