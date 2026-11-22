########################## set the work directory ##########################
rm(list=ls())
setwd("Path/to/your_workdir")
getwd()

############################ Load R package and data ###########################
cran_packages <- c(
  "readxl", "dplyr", "writexl", "lme4", "lmerTest", "broom.mixed",
  "tidyverse", "ggplot2", "patchwork", "ggrepel", "cowplot"
  )

# Only install packages that are not already installed
new_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(readxl)
library(dplyr)
library(writexl)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(cowplot)


tlens <- read_xlsx("tlens_chr_filtered.xlsx")
clinical_info <- read_xlsx("31family_info.xlsx", range = "A1:G125")    #clinical phenotype data: family_ID, sample_ID, family_role, sex, age


################## Combine clinical information with ChArmTLs ##################
######## Calculate the mTL of individuals and the average length of chromosome arms
tlens_dip <- tlens %>%
  group_by(sample, chr) %>%
  summarise(TL_p75 = mean(TL_p75, na.rm = TRUE), .groups = "drop") %>%
  mutate(level = "diploid") %>%
  select(sample, level, chr, TL_p75) %>%
  bind_rows(tlens %>%
              group_by(sample) %>%
              summarise(TL_p75 = mean(TL_p75, na.rm = TRUE), .groups = "drop") %>%
              mutate(level = "diploid", chr = "mTL") %>%
              select(sample, level, chr, TL_p75))


######## Calculate the haplotype-specific mTL of individuals and extract the ChArmTLs
tlens_hap <- tlens %>%
  filter(family_role %in% c("ASD","Sibling"), parent_of_origin %in% c("Pat","Mat")) %>%
  select(sample, parent_of_origin, chr, TL_p75) %>%
  rename(level = parent_of_origin) %>%
  bind_rows(tlens %>%
              filter(family_role %in% c("ASD","Sibling"), parent_of_origin %in% c("Pat","Mat")) %>%
              group_by(sample, parent_of_origin) %>%
              summarise(TL_p75 = mean(TL_p75, na.rm = TRUE), .groups = "drop") %>%
              rename(level = parent_of_origin) %>%
              mutate(chr = "mTL") %>%
              select(sample, level, chr, TL_p75))


######## merge data
tlens_clin <- bind_rows(tlens_dip, tlens_hap) %>%
  mutate(level = factor(level, levels = c("diploid","Pat","Mat")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%
  arrange(sample, level, chr)

tlens_clin <- merge(tlens_clin, clinical_info, by = "sample")
write_xlsx(tlens_clin, "tlens_clin.xlsx")


############################ Sex and age differences ###########################
######## data preparation
tlens_clin <- tlens_clin %>%
  mutate(diagnosis = case_when(family_role == "ASD" ~ "ASD",
                               family_role == "Sibling" ~ "normal",
                               TRUE ~ NA_character_),
         sex = factor(sex, levels = c("Male","Female")),
         diagnosis = factor(diagnosis, levels = c("normal","ASD"))) %>%
  relocate(diagnosis, .after = family_role)


######## Linear mixed model, sex, age and family ID
fit_model_1 <- function(data) {
  data <- data %>% filter(!is.na(TL_p75))
  
  model <- tryCatch({
    lmer(TL_p75 ~ sex + age + (1 | family), data = data)
  }, warning = function(w) {
    suppressWarnings(lmer(TL_p75 ~ sex + age + (1 | family), data = data))
  }, error = function(e) return(NULL))
  
  if (is.null(model)) return(NULL)
  
  model_output <- tidy(model, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>%
    mutate(n_obs = nrow(data), singular_flag = isSingular(model)) %>%
    filter(term != "(Intercept)")
  
  return(model_output)
}

#### all samples
diff_all <- tlens_clin %>%
  filter(level == "diploid") %>%
  group_by(chr) %>%
  nest() %>%
  mutate(model_summary = map(data, fit_model_1)) %>%
  select(-data) %>%
  unnest(cols = c(model_summary))

diff_all <- diff_all %>%
  group_by(term) %>%
  mutate(p_adj_fdr = case_when(chr == "mTL" ~ p.value,
                               TRUE ~ p.adjust(p.value, method = "fdr"))) %>%
  ungroup() %>%
  mutate(level = "diploid_all",
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%
  arrange(chr) %>%
  select(level,chr,n_obs,term,estimate,std.error,conf.low,conf.high,statistic,p.value,p_adj_fdr,singular_flag,everything())


######## Linear mixed model, sex, age, disease status and family ID
fit_model_2 <- function(data) {
  data <- data %>% filter(!is.na(TL_p75))
  
  model <- tryCatch({
    lmer(TL_p75 ~ sex + age + diagnosis + (1 | family), data = data)
  }, warning = function(w) {
    suppressWarnings(lmer(TL_p75 ~ sex + age + diagnosis + (1 | family), data = data))
  }, error = function(e) return(NULL))
  
  if (is.null(model)) return(NULL)
  
  model_output <- tidy(model, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>%
    mutate(n_obs = nrow(data), singular_flag = isSingular(model)) %>%
    filter(term != "(Intercept)")
  
  return(model_output)
}

#### all offspring
diff_offspring <- tlens_clin %>%
  filter(family_role %in% c("ASD","Sibling")) %>% 
  group_by(level,chr) %>%
  nest() %>%
  mutate(model_summary = map(data, fit_model_2)) %>%
  select(-data) %>%
  unnest(cols = c(model_summary))

diff_offspring <- diff_offspring %>%
  group_by(level,term) %>%
  mutate(p_adj_fdr = case_when(chr == "mTL" ~ p.value,
                               TRUE ~ p.adjust(p.value, method = "fdr"))) %>%
  ungroup() %>%
  mutate(level = as.character(level),
         level = ifelse(level == "diploid", "diploid_offspring", level),
         level = factor(level, levels = c("diploid_offspring","Pat","Mat")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%
  arrange(level, chr) %>%
  select(level,chr,n_obs,term,estimate,std.error,conf.low,conf.high,statistic,p.value,p_adj_fdr,singular_flag,everything())


######## merge data
diff <- rbind(diff_all, diff_offspring)
write_xlsx(diff, "diff.xlsx")


################## Visualization of gender and age differences #################
######## all samples
diff_all_plot <- diff_all %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         sig = case_when(p_adj_fdr < 0.05 ~ "FDR < 0.05",
                         p.value < 0.05 ~ "Nominal P < 0.05",
                         TRUE ~ "Not significant"),
         sig = factor(sig, levels = c("FDR < 0.05","Nominal P < 0.05","Not significant")))

#### sex
p1 <- ggplot(data = subset(diff_all_plot, term == "sexFemale"), aes(x = estimate, y = chr_number, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "black") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0.5, linewidth = 0.35) +
  geom_point(size = 0.75) +
  scale_color_manual(values = c("FDR < 0.05" = "#EB6368", "Nominal P < 0.05" = "#FCD28c", "Not significant"  = "grey80"), name = NULL) +
  scale_x_continuous(breaks = c(-500,0,500,1000), labels = c("-500","0","500","1,000")) +
  scale_y_discrete(limits = rev(c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                                  "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q", "mTL")))+
  labs(x = "Female - Male ChArmTL (bp)", y = NULL) +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        legend.position = "bottom",
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.25, unit = "cm"))

#### age
p2 <- ggplot(data = subset(diff_all_plot, term == "age"), aes(x = estimate, y = chr_number, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "black") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0.5, linewidth = 0.35) +
  geom_point(size = 0.75) +
  scale_color_manual(values = c("FDR < 0.05" = "#EB6368", "Nominal P < 0.05" = "#FCD28c", "Not significant"  = "grey80"), name = NULL) +
  scale_x_reverse() +
  scale_y_discrete(limits = rev(c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                                  "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q", "mTL")))+
  labs(x = "Annual change in ChArmTL (bp/year)", y = NULL) +
  theme_bw() +
  guides(color = "none") + 
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        # legend.text = element_text(size = 6),
        # legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        # legend.box.spacing = unit(0.01, "cm"),
        # legend.position = "bottom",
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.25, r = 0.25, unit = "cm"))


######## offspring Pat + Mat
diff_offspring_hap_plot <- diff_offspring %>%
  filter(level %in% c("Pat","Mat")) %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         sig = case_when(p_adj_fdr < 0.05 ~ "FDR < 0.05",
                         p.value < 0.05 ~ "Nominal P < 0.05",
                         TRUE ~ "Not significant"),
         sig = factor(sig, levels = c("FDR < 0.05","Nominal P < 0.05","Not significant")))

#### sex
p3 <- ggplot(data = subset(diff_offspring_hap_plot, term == "sexFemale"), aes(x = estimate, y = chr_number, shape = level, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "black") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0.5, linewidth = 0.35, position = position_dodge(width = 0.75)) +
  geom_point(size = 0.75, position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c("Pat" = 16, "Mat" = 17), labels = c("Pat", "Mat") , name = NULL) +
  scale_color_manual(values = c("FDR < 0.05" = "#EB6368", "Nominal P < 0.05" = "#FCD28c", "Not significant"  = "grey80"), name = NULL) +
  scale_x_continuous(breaks = c(-4000,-2000,0,2000,4000), labels = c("-4,000","-2,000","0","2,000","4,000")) +
  scale_y_discrete(limits = rev(c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                                  "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q", "mTL")))+
  labs(x = "Female - Male ChArmTL (bp)", y = NULL) +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        legend.position = "bottom",
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.25, unit = "cm"))

#### age
p4 <- ggplot(data = subset(diff_offspring_hap_plot, term == "age"), aes(x = estimate, y = chr_number, shape = level, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "black") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0.5, linewidth = 0.35, position = position_dodge(width = 0.75)) +
  geom_point(size = 0.75, position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c("Pat" = 16, "Mat" = 17), labels = c("Pat", "Mat") , name = NULL) +
  scale_color_manual(values = c("FDR < 0.05" = "#EB6368", "Nominal P < 0.05" = "#FCD28c", "Not significant"  = "grey80"), name = NULL) +
  scale_y_discrete(limits = rev(c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                                  "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q", "mTL")))+
  labs(x = "Annual change in ChArmTL (bp/year)", y = NULL) +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        legend.position = "bottom",
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.25, r = 0.5, unit = "cm"))

p3_4 <- (p3 + p4) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", 
        plot.tag = element_text(face = "bold", size = 8))

ggsave(filename = "Supplementary Figure 7.pdf", plot = p3_4, width = 10, height = 12, units = "cm", dpi = 600)


############## Consistency in direction of sex and age differences #############
######## Group by "level"
level <- unique(diff$level)
term_alt <- c(sexFemale = "greater",
              age       = "less")

diff_direction <- list()
for (tm in names(term_alt)) {
  diff_direction[[tm]] <- list()
  for (g in level) {
    sub <- diff %>% filter(level == g, chr != "mTL", term == tm)
    pos <- sum(sub$estimate > 0)
    neg <- sum(sub$estimate < 0)
    total <- pos + neg
    test <- binom.test(pos, total, p = 0.5, alternative = term_alt[[tm]])
    diff_direction[[tm]][[g]] <- data.frame(
      group = g,
      term = tm,
      positive = pos,
      negative = neg,
      total = total,
      prop_positive = pos / total,
      p_value = test$p.value,
      lower_ci = test$conf.int[1],
      upper_ci = test$conf.int[2])
  }
}

rm(level, term_alt, tm, g, sub, pos, neg, total, test)

#### merge data
diff_direction <- do.call(rbind, lapply(diff_direction, function(x) do.call(rbind, x)))
rownames(diff_direction) <- NULL

write_xlsx(diff_direction, "diff_direction.xlsx")


###### The correlation of sex and age among different groups of offspring ######
######## sex
diff_corr_sex_wide <- diff_offspring %>%
  filter(term == "sexFemale") %>%
  select(level, chr, estimate) %>%
  pivot_wider(names_from = level, values_from = estimate) %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]"))) %>%
  filter( chr != "mTL")

#### Descriptive check: Scatter plots are more crucial than Shapiro
pairs(diff_corr_sex_wide[, c("Pat","Mat")],
      main = "Pairwise scatter matrix")

## Shapiro, as a distribution description
shapiro.test(diff_corr_sex_wide$Pat)
shapiro.test(diff_corr_sex_wide$Mat)

#### correlation test
cor.test(diff_corr_sex_wide$Pat, diff_corr_sex_wide$Mat, method = "spearman", 
         alternative = "two.sided", exact = FALSE, conf.level = 0.95)

p5 <- ggplot(diff_corr_sex_wide, aes(x = Pat, y = Mat, label = chr_number))+
  geom_point(size = 0.5, color = "#82ADD0")+
  geom_text_repel(size = 1.5, max.overlaps = Inf, box.padding = 0.05, point.padding = 0.05, segment.color = "grey40")+
  geom_smooth(method = 'lm', formula = y ~ x, se = T, linewidth = 0.35, color = "grey40")+
  geom_rug(outside = F, length = unit(0.02, "npc"), linewidth = 0.2, color = "#82ADD0")+
  scale_x_continuous(expand = c(0,0), limits = c(-500,1800), breaks = c(-500,0,500,1000,1500), labels = c("-500","0","500","1,000","1,500"))+
  scale_y_continuous(expand = c(0,0), limits = c(-2200,2000), breaks = c(-2000,-1000,0,1000,2000), labels = c("-2,000","-1,000","0","1,000","2,000"))+
  labs(x = "Paternally inherited\n(Female - Male ChArmTL, bp)", y = "Maternally inherited\n(Female - Male ChArmTL, bp)")+
  annotate(geom = "text", x = -400, y = 1800, label = "rho = 0.028, P = 0.87", size = 2, hjust = 0)+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.5, b = 0.25, l = 0.25, r = 0.5, unit = "cm"))


######## age
diff_corr_age_wide <- diff_offspring %>%
  filter(term == "age") %>%
  select(level, chr, estimate) %>%
  pivot_wider(names_from = level, values_from = estimate) %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]"))) %>%
  filter( chr != "mTL")

#### Descriptive check: Scatter plots are more crucial than Shapiro
pairs(diff_corr_age_wide[, c("Pat","Mat")],
      main = "Pairwise scatter matrix")

## Shapiro, as a distribution description
shapiro.test(diff_corr_age_wide$diploid_offspring)
shapiro.test(diff_corr_age_wide$Pat)
shapiro.test(diff_corr_age_wide$Mat)

#### correlation test
cor.test(diff_corr_age_wide$Pat, diff_corr_age_wide$Mat, method = "spearman", 
         alternative = "two.sided", exact = FALSE, conf.level = 0.95)

p6 <- ggplot(diff_corr_age_wide, aes(x = Pat, y = Mat, label = chr_number))+
  geom_point(size = 0.5, color = "#82ADD0")+
  geom_text_repel(size = 1.5, max.overlaps = Inf, box.padding = 0.05, point.padding = 0.05, segment.color = "grey40")+
  geom_smooth(method = 'lm', formula = y ~ x, se = T, linewidth = 0.35, color = "grey40")+
  geom_rug(outside = F, length = unit(0.02, "npc"), linewidth = 0.2, color = "#82ADD0")+
  scale_x_continuous(expand = c(0,0), limits = c(-210,70), breaks = c(-200,-150,-100,-50,0,50))+
  scale_y_continuous(expand = c(0,0), limits = c(-210,40), breaks = c(-200,-150,-100,-50,0))+
  labs(x = "Paternally inherited\nannual change in ChArmTL (bp/year)", y = "Maternally inherited\nannual change in ChArmTL (bp/year)")+
  annotate(geom = "text", x = -180, y = 20, label = "rho = -0.15, P = 0.36", size = 2, hjust = 0)+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.25, r = 0.5, unit = "cm"))

p1_2 <- plot_grid(p1, p2, nrow = 1, labels = c('a','b'), label_size = 8, rel_widths = c(1,1), align = 'h')
p5_6 <- plot_grid(p5, p6, NULL, ncol = 1, labels = c('c', 'd', ''), label_size = 8, rel_heights = c(1,1,0.1), align = 'v')

p1_2_5_6 <- plot_grid(p1_2, p5_6, nrow = 1, labels = c('', ''), rel_widths = c(1,0.75), align = 'hv')
ggsave(filename = "Figure 2.pdf", plot = p1_2_5_6, width = 15, height = 12, units = "cm", dpi = 600)

