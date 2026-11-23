############################ set the work directory ############################
rm(list=ls())
setwd("Path/to/your_workdir")
getwd()


############################ Load R package and data ###########################
# Define the list of required packages
required_packages <- c("readxl", "dplyr","writexl","lme4","lmerTest","emmeans","ggeffects",
                       "ggplot2","tidyverse","broom","ggrepel","cowplot","ggh4x")

# Identify packages from the list that are not currently installed
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

# If there are missing packages, download and install them
if(length(missing_packages) > 0) {
  message("Installing missing packages...")
  install.packages(missing_packages, dependencies = TRUE)
} else {
  message("All required packages are already installed.")
}

library(readxl)
library(dplyr)
library(writexl)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggeffects)
library(ggplot2)
library(tidyverse)
library(broom)
library(ggrepel)
library(cowplot)
library(ggh4x)

tlens <- read_xlsx("tlens_chr_filtered.xlsx")
clinical_info <- read_xlsx("31family_info.xlsx", range = "A1:G125")    #clinical phenotype data: family_ID, sample_ID, family_role, sex, age


################## Combine clinical information with ChArmTLs ##################
######## Calculate the haplotype-specific mTL and extract the ChArmTLs
tlens_hap <- tlens %>%
  filter(family_role %in% c("Father","Mother"), trans_status_ASD %in% c("T-ASD","U-ASD")) %>%
  select(sample, trans_status_ASD, chr, TL_p75) %>%
  rename(haplotype = trans_status_ASD) %>%
  bind_rows(tlens %>%
              filter(family_role %in% c("Father","Mother"), trans_status_ASD %in% c("T-ASD","U-ASD")) %>%
              group_by(sample, trans_status_ASD) %>%
              summarise(TL_p75 = mean(TL_p75, na.rm = TRUE), .groups = "drop") %>%
              mutate(chr = "mTL") %>%
              select(sample, trans_status_ASD, chr, TL_p75) %>%
              rename(haplotype = trans_status_ASD)) %>%
  bind_rows(tlens %>%
              filter(family_role %in% c("Father","Mother"), trans_status_Sib %in% c("T-Sib","U-Sib")) %>%
              select(sample, trans_status_Sib, chr, TL_p75) %>%
              rename(haplotype = trans_status_Sib)) %>%
  bind_rows(tlens %>%
              filter(family_role %in% c("Father","Mother"), trans_status_Sib %in% c("T-Sib","U-Sib")) %>%
              group_by(sample, trans_status_Sib) %>%
              summarise(TL_p75 = mean(TL_p75, na.rm = TRUE), .groups = "drop") %>%
              mutate(chr = "mTL") %>%
              select(sample, trans_status_Sib, chr, TL_p75) %>%
              rename(haplotype = trans_status_Sib)) %>%
  bind_rows(tlens %>%
              filter(family_role %in% c("ASD","Sibling"), parent_of_origin %in% c("Pat","Mat")) %>%
              select(sample, parent_of_origin, chr, TL_p75) %>%
              rename(haplotype = parent_of_origin)) %>%
  bind_rows(tlens %>%
              filter(family_role %in% c("ASD","Sibling"), parent_of_origin %in% c("Pat","Mat")) %>%
              group_by(sample, parent_of_origin) %>%
              summarise(TL_p75 = mean(TL_p75, na.rm = TRUE), .groups = "drop") %>%
              mutate(chr = "mTL") %>%
              select(sample, parent_of_origin, chr, TL_p75) %>%
              rename(haplotype = parent_of_origin)) %>%
  mutate(haplotype  = factor(haplotype , levels = c("T-Sib","U-Sib","T-ASD","U-ASD","Pat","Mat")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%
  arrange(sample, haplotype, chr)


######## merge data
tlens_hap <- merge(tlens_hap, clinical_info, by = "sample")
tlens_hap$family_role <- factor(tlens_hap$family_role, levels = c("Father","Mother","Sibling","ASD"))
tlens_hap$sex <- factor(tlens_hap$sex, levels = c("Male","Female"))
write_xlsx(tlens_hap, "tlens_hap.xlsx")


####################### The distribution of Δ(Pat - Mat) #######################
######## T-Sib - U-Sib
tlens_T_U_Sib_diff <- tlens_hap %>%
  filter(family_role %in% c("Father","Mother"), haplotype %in% c("T-Sib","U-Sib")) %>%
  pivot_wider(names_from = haplotype, values_from = TL_p75) %>%
  drop_na(`T-Sib`,`U-Sib`) %>%
  mutate(diff = `T-Sib`-`U-Sib`)

tlens_T_U_Sib_diff <- tlens_T_U_Sib_diff %>%
  mutate(family_role = "Parents") %>%
  bind_rows(tlens_T_U_Sib_diff) %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         family_role = factor(family_role, levels = c("Parents","Father","Mother")))

T_U_Sib_diff_stats <- tlens_T_U_Sib_diff %>%
  group_by(family_role, chr) %>%
  summarise(mean_diff = mean(diff, na.rm = TRUE),
            sd_diff = sd(diff, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")))

p1 <- ggplot()+
  geom_errorbar(data = T_U_Sib_diff_stats, aes(x = chr_number, ymin = mean_diff - sd_diff, ymax = mean_diff + sd_diff, color = family_role),
                width = 0.4, linewidth = 0.25, position = position_dodge(width = 0.75))+
  stat_summary(data = tlens_T_U_Sib_diff, aes(x = chr_number, y = diff, color = family_role), 
               geom = "point", fun = "mean", size = 0.35, position = position_dodge(width = 0.75))+
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "black")+
  scale_color_manual(values = c("Parents" = "grey80", "Father" = "#FCD28c", "Mother" = "#A1B5A3"), name = NULL) +
  scale_y_continuous(breaks = c(-4000,-2000,0,2000,4000), labels = c("-4,000","-2,000","0","2,000","4,000")) +
  scale_x_discrete(limits = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                              "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q", "mTL"))+
  labs(x = NULL, y = "T-Sib - U-Sib\nChArmTL (bp)") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        legend.position = "top",
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.25, l = 0.5, r = 0.5, unit = "cm"))


######## T-ASD - U-ASD
tlens_T_U_ASD_diff <- tlens_hap %>%
  filter(family_role %in% c("Father","Mother"), haplotype %in% c("T-ASD","U-ASD")) %>%
  pivot_wider(names_from = haplotype, values_from = TL_p75) %>%
  drop_na(`T-ASD`, `U-ASD`) %>%
  mutate(diff = `T-ASD`-`U-ASD`)

tlens_T_U_ASD_diff <- tlens_T_U_ASD_diff %>%
  mutate(family_role = "Parents") %>%
  bind_rows(tlens_T_U_ASD_diff) %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         family_role = factor(family_role, levels = c("Parents","Father","Mother")))

T_U_ASD_diff_stats <- tlens_T_U_ASD_diff %>%
  group_by(family_role, chr) %>%
  summarise(mean_diff = mean(diff, na.rm = TRUE),
            sd_diff = sd(diff, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")))

p2 <- ggplot()+
  geom_errorbar(data = T_U_ASD_diff_stats, aes(x = chr_number, ymin = mean_diff - sd_diff, ymax = mean_diff + sd_diff, color = family_role),
                width = 0.4, linewidth = 0.25, position = position_dodge(width = 0.75))+
  stat_summary(data = tlens_T_U_ASD_diff, aes(x = chr_number, y = diff, color = family_role), 
               geom = "point", fun = "mean", size = 0.35, position = position_dodge(width = 0.75))+
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "black")+
  scale_color_manual(values = c("Parents" = "grey80", "Father" = "#F9D5B2", "Mother" = "#87BEA7"), name = NULL) +
  scale_y_continuous(breaks = c(-6000,-3000,0,3000,6000), labels = c("-6,000","-3,000","0","3,000","6,000")) +
  scale_x_discrete(limits = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                              "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q", "mTL"))+
  labs(x = NULL, y = "T-ASD - U-ASD\nChArmTL (bp)") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        legend.position = "top",
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))

p1_2 <- plot_grid(p1, p2, ncol = 1, labels = c('a','b'), label_size = 8, rel_heights = c(1,1), align = 'v', axis = 'lr')
ggsave(filename = "Supplementary Figure 10.pdf", plot = p1_2, width = 15, height = 10, units = "cm", dpi = 600)


########### Analysis of differences between transmit and untransmit ############
######## A custom function compares the telomere lengths of the transmitted and untransmitted status
T_U_difference <- function(data, t_col, u_col) {
  data %>%
    group_by(family_role, chr) %>%
    summarise(
      n_pairs = sum(!is.na(.data[[t_col]]) & !is.na(.data[[u_col]])),
      mean_T  = mean(.data[[t_col]], na.rm = TRUE),
      sd_T    = sd(.data[[t_col]], na.rm = TRUE),
      mean_U  = mean(.data[[u_col]], na.rm = TRUE),
      sd_U    = sd(.data[[u_col]], na.rm = TRUE),
      tidy(t.test(.data[[t_col]], .data[[u_col]], paired = TRUE)),
      .groups = 'drop'
    )
}


######## Comparison of differences between transmitted and untransmitted status
diff_T_U_Sib <- T_U_difference(tlens_T_U_Sib_diff, "T-Sib","U-Sib")
diff_T_U_ASD <- T_U_difference(tlens_T_U_ASD_diff, "T-ASD","U-ASD")

#### Merge data and add corrected P-values
diff_T_U_Sib <- diff_T_U_Sib %>%
  group_by(family_role) %>%
  mutate(p_adj_fdr = case_when(chr == "mTL" ~ p.value,
                               TRUE ~ p.adjust(p.value, method = "fdr")),
         family_role= factor(family_role, levels = c("Mother","Father","Parents")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%
  select(family_role,chr,n_pairs,mean_T,sd_T,mean_U,sd_U,estimate,conf.low,conf.high,statistic,p.value,p_adj_fdr,method,alternative) %>%
  rename(mean_diff = estimate) %>%
  arrange(family_role, chr)

write_xlsx(diff_T_U_Sib, "diff_T_U_Sib.xlsx")


diff_T_U_ASD <- diff_T_U_ASD %>%
  group_by(family_role) %>%
  mutate(p_adj_fdr = case_when(chr == "mTL" ~ p.value,
                               TRUE ~ p.adjust(p.value, method = "fdr")),
         family_role= factor(family_role, levels = c("Mother","Father","Parents")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%
  select(family_role,chr,n_pairs,mean_T,sd_T,mean_U,sd_U,estimate,conf.low,conf.high,statistic,p.value,p_adj_fdr,method,alternative) %>%
  rename(mean_diff = estimate) %>%
  arrange(family_role, chr)

write_xlsx(diff_T_U_ASD, "diff_T_U_ASD.xlsx")


########## Visualization of transmitted and untransmitted differences ##########
diff_T_U_Sib_plot <- diff_T_U_Sib %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         mean_diff_sign = ifelse(mean_diff > 0, "T-Sib > U-Sib", "T-Sib < U-Sib"),
         mean_diff_sign = factor(mean_diff_sign, levels =c("T-Sib > U-Sib", "T-Sib < U-Sib")),
         significance = case_when(p.value < 0.05 ~ "Nominal P < 0.05",
                                  TRUE ~ "Not significant"),
         significance = factor(significance, levels = c("Nominal P < 0.05", "Not significant")))

p3 <- ggplot(diff_T_U_Sib_plot, aes(x = chr_number, y = family_role, color = significance, shape = mean_diff_sign, size = abs(mean_diff)))+
  geom_point(stroke = 0.35)+
  scale_color_manual(values = c("Nominal P < 0.05" = "#FCD28c", "Not significant"  = "grey80"), name = NULL)+
  scale_shape_manual(values = c("T-Sib > U-Sib" = 16, "T-Sib < U-Sib" = 1), name = NULL)+
  scale_size_continuous(name = "|Mean difference|", range = c(0.5,3), breaks = c(500,1000,1500), labels = c("500","1,000","1,500"))+
  scale_x_discrete(limits = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                              "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q", "mTL"))+
  labs(x = NULL, y = NULL)+
  guides(shape = guide_legend(override.aes = list(size = 1.5), position = "top", order = 1),
         color = guide_legend(override.aes = list(size = 1.5), position = "top", order = 2), 
         size = guide_legend(override.aes = list(shape = 1, stroke = 0.35), position = "top", order = 3))+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        panel.grid = element_line(linewidth = 0.2),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.25, l = 0.5, r = 0.5, unit = "cm"))


diff_T_U_ASD_plot <- diff_T_U_ASD %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         mean_diff_sign = ifelse(mean_diff > 0, "T-ASD > U-ASD", "T-ASD < U-ASD"),
         mean_diff_sign = factor(mean_diff_sign, levels =c("T-ASD > U-ASD", "T-ASD < U-ASD")),
         significance = case_when(p.value < 0.05 ~ "Nominal P < 0.05",
                                  TRUE ~ "Not significant"),
         significance = factor(significance, levels = c("Nominal P < 0.05", "Not significant")))

p4 <- ggplot(diff_T_U_ASD_plot, aes(x = chr_number, y = family_role, color = significance, shape = mean_diff_sign, size = abs(mean_diff)))+
  geom_point(stroke = 0.35)+
  scale_color_manual(values = c("Nominal P < 0.05" = "#FCD28c", "Not significant"  = "grey80"), name = NULL)+
  scale_shape_manual(values = c("T-ASD > U-ASD" = 16, "T-ASD < U-ASD" = 1), name = NULL)+
  scale_size_continuous(name = "|Mean difference|", range = c(0.5,3), breaks = c(500,1000,1500,2000), labels = c("500","1,000","1,500","2,000"))+
  scale_x_discrete(limits = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                              "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q", "mTL"))+
  labs(x = NULL, y = NULL)+
  guides(shape = guide_legend(override.aes = list(size = 1.5), position = "top", order = 1),
         color = guide_legend(override.aes = list(size = 1.5), position = "top", order = 2), 
         size = guide_legend(override.aes = list(shape = 1, stroke = 0.35), position = "top", order = 3))+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        panel.grid = element_line(linewidth = 0.2),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.5, r = 0.5, unit = "cm"))


########## The correlation between father and mother in transmission status differences ##########
diff_corr_T_U_Sib <- diff_T_U_Sib %>%
  select(chr,family_role,mean_diff) %>%
  filter(!family_role == "Parents") %>%
  pivot_wider(names_from = family_role, values_from = mean_diff) %>% 
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]"))) %>%
  filter( chr != "mTL") 

shapiro.test(diff_corr_T_U_Sib$Mother)
shapiro.test(diff_corr_T_U_Sib$Father)

cor.test(diff_corr_T_U_Sib$Mother, diff_corr_T_U_Sib$Father, alternative = "two.sided", method = "spearman", exact = TRUE, conf.level = 0.95)

p5 <- ggplot(diff_corr_T_U_Sib, aes(x = Father, y = Mother, label = chr_number))+
  geom_point(size = 0.5, color = "#FCD28c")+
  geom_text_repel(size = 1.5, max.overlaps = Inf, box.padding = 0.05, point.padding = 0.05, segment.color = "grey40")+
  geom_smooth(method = 'lm', formula = y ~ x, se = T, linewidth = 0.35, color = "grey40")+
  geom_rug(outside = F, length = unit(0.02, "npc"), linewidth = 0.2, color = "#FCD28c")+
  scale_x_continuous(expand = c(0,0), limits = c(-1000,2000), breaks = c(-1000,0,1000,2000), labels = c("-1,000","0","1,000","2,000"))+
  scale_y_continuous(expand = c(0,0), limits = c(-600,2000), breaks = c(0,500,1000,1500,2000), labels = c("0","500","1,000","1,500","2,000"))+
  labs(x = "T-Sib - U-Sib ChArmTLs in father (bp)", y = "T-Sib - U-Sib ChArmTLs in mother (bp)")+
  annotate(geom = "text", x = -800, y = 1800, label = "rho = -0.011, P = 0.95", size = 2, hjust = 0)+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.5, r = 0.25, unit = "cm"))


diff_corr_T_U_ASD <- diff_T_U_ASD %>%
  select(chr,family_role,mean_diff) %>%
  filter(!family_role == "Parents") %>%
  pivot_wider(names_from = family_role, values_from = mean_diff) %>% 
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]"))) %>%
  filter( chr != "mTL") 

shapiro.test(diff_corr_T_U_ASD$Mother)
shapiro.test(diff_corr_T_U_ASD$Father)

cor.test(diff_corr_T_U_ASD$Mother, diff_corr_T_U_ASD$Father, alternative = "two.sided", method = "spearman", exact = TRUE, conf.level = 0.95)

p6 <- ggplot(diff_corr_T_U_ASD, aes(x = Father, y = Mother, label = chr_number))+
  geom_point(size = 0.5, color = "#A1B5A3")+
  geom_text_repel(size = 1.5, max.overlaps = Inf, box.padding = 0.05, point.padding = 0.05, segment.color = "grey40")+
  geom_smooth(method = 'lm', formula = y ~ x, se = T, linewidth = 0.35, color = "grey40")+
  geom_rug(outside = F, length = unit(0.02, "npc"), linewidth = 0.2, color = "#A1B5A3")+
  scale_x_continuous(expand = c(0,0), limits = c(-1100,2200), breaks = c(-1000,0,1000,2000), labels = c("-1,000","0","1,000","2,000"))+
  scale_y_continuous(expand = c(0,0), limits = c(-2600,1200), breaks = c(-2000,-1000,0,1000), labels = c("-2,000","-1,000","0","1,000"))+
  labs(x = "T-ASD - U-ASD ChArmTLs in father (bp)", y = "T-ASD - U-ASD ChArmTLs in mother (bp)")+
  annotate(geom = "text", x = -800, y = 1000, label = "rho = 0.0058, P = 0.97", size = 2, hjust = 0)+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.25, r = 0.25, unit = "cm"))


#### The direction of transmitted and untransmitted differences (by chromosome arms) ####
######## Define the directional statistical function
diff_direction_chr <- function(data, t_col, u_col) {
  data %>%
    group_by(family_role, chr) %>%
    summarise(
      n_total    = sum(!is.na(.data[[t_col]]) & !is.na(.data[[u_col]])),
      n_T_longer = sum(.data[[t_col]] > .data[[u_col]], na.rm = TRUE),
      n_U_longer = sum(.data[[t_col]] < .data[[u_col]], na.rm = TRUE),
      n_Equal    = sum(.data[[t_col]] == .data[[u_col]], na.rm = TRUE),
      tidy(binom.test(n_T_longer, n_total, p = 0.5, alternative = "greater")),
      .groups = 'drop'
    )
}  


######## Directional statistics of transmitted and untransmitted differences
direction_chr_by_Sib <- diff_direction_chr(tlens_T_U_Sib_diff,"T-Sib","U-Sib")
direction_chr_by_ASD <- diff_direction_chr(tlens_T_U_ASD_diff,"T-ASD","U-ASD")

direction_chr_by_Sib <- direction_chr_by_Sib %>%
  mutate(family_role = factor(family_role, levels = c("Parents","Father","Mother")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%   
  select(family_role,chr,n_total,n_T_longer,n_U_longer,n_Equal,estimate,conf.low,conf.high,statistic,p.value,method,alternative) %>%
  rename(proportion_T_longer = estimate) %>%
  arrange(family_role, chr)

write_xlsx(direction_chr_by_Sib, "direction_chr_by_Sib.xlsx")

direction_chr_by_ASD <- direction_chr_by_ASD %>%
  mutate(family_role = factor(family_role, levels = c("Parents","Father", "Mother")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%   
  select(family_role,chr,n_total,n_T_longer,n_U_longer,n_Equal,estimate,conf.low,conf.high,statistic,p.value,method,alternative) %>%
  rename(proportion_T_longer = estimate) %>%
  arrange(family_role, chr)

write_xlsx(direction_chr_by_ASD, "direction_chr_by_ASD.xlsx")


########## Visualization of the direction of transmitted and untransmitted differences (by chromosome arms) ##########
direction_chr_by_Sib_plot <- direction_chr_by_Sib %>%
  select(family_role, chr, n_T_longer, n_U_longer) %>%
  pivot_longer(cols = c(n_T_longer, n_U_longer),
               names_to = "direction",
               values_to = "counts") %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         direction = factor(direction, levels = c("n_T_longer", "n_U_longer")),
         family_role= factor(family_role, levels = c("Parents","Father","Mother")))

p7 <- ggplot(direction_chr_by_Sib_plot, aes(x = chr_number, y = counts, color = direction))+
  geom_linerange(aes(ymin = 0, ymax = counts), linewidth = 0.35, position = position_dodge(width = 0.5))+
  geom_point(size = 1, position = position_dodge(width = 0.5))+
  facet_wrap2(family_role ~ ., nrow = 3, scales = "free_y", strip.position = "right",
              strip = strip_themed(background_y = elem_list_rect(fill = c("grey80","#FCD28c80","#A1B5A380"))))+
  facetted_pos_scales(y = list(family_role == "Father" ~ scale_y_continuous(expand = c(0,0), limits = c(0,30), breaks = c(10,20,30)),
                               family_role == "Mother"     ~ scale_y_continuous(expand = c(0,0), limits = c(0,30), breaks = c(10,20,30)),
                               family_role == "Parents" ~ scale_y_continuous(expand = c(0,0), limits = c(0,60), breaks = c(20,40,60))))+
  scale_color_manual(values = c("#EB6368","grey80"), labels = c("T-Sib > U-Sib","T-Sib < U-Sib"))+
  scale_x_discrete(limits = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                              "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q","mTL"))+
  labs(x = NULL, y = "Number of samples", color = NULL)+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        legend.position = "top",
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.35),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.25, l = 0.5, r = 0.5, unit = "cm"))


direction_chr_by_ASD_plot <- direction_chr_by_ASD %>%
  select(family_role, chr, n_T_longer, n_U_longer) %>%
  pivot_longer(cols = c(n_T_longer, n_U_longer),
               names_to = "direction",
               values_to = "counts") %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         direction = factor(direction, levels = c("n_T_longer", "n_U_longer")),
         family_role= factor(family_role, levels = c("Parents","Father","Mother")))

p8 <- ggplot(direction_chr_by_ASD_plot, aes(x = chr_number, y = counts, color = direction))+
  geom_linerange(aes(ymin = 0, ymax = counts), linewidth = 0.35, position = position_dodge(width = 0.5))+
  geom_point(size = 1, position = position_dodge(width = 0.5))+
  facet_wrap2(family_role ~ ., nrow = 3, scales = "free_y", strip.position = "right",
              strip = strip_themed(background_y = elem_list_rect(fill = c("grey80","#F9D5B280","#87BEA780"))))+
  facetted_pos_scales(y = list(family_role == "Father" ~ scale_y_continuous(expand = c(0,0), limits = c(0,30), breaks = c(10,20,30)),
                               family_role == "Mother"     ~ scale_y_continuous(expand = c(0,0), limits = c(0,30), breaks = c(10,20,30)),
                               family_role == "Parents" ~ scale_y_continuous(expand = c(0,0), limits = c(0,60), breaks = c(20,40,60))))+
  scale_color_manual(values = c("#EB6368","grey80"), labels = c("T-ASD > U-ASD","T-ASD < U-ASD"))+
  scale_x_discrete(limits = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                              "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q","mTL"))+
  labs(x = NULL, y = "Number of samples", color = NULL)+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        legend.position = "top",
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.35),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))

p7_8 <- plot_grid(p7, p8, labels = c('a','b'), label_size = 8, ncol = 1, rel_widths = c(1,1), align = 'h', axis = 'tb')

ggsave(filename = "Supplementary Figure 11.pdf", plot = p7_8, width = 15, height = 18, units = "cm", dpi = 600)


######################## The direction of transmitted and untransmitted differences (by individual) ########################
direction_sample_by_Sib <- tlens_T_U_Sib_diff %>%
  filter(chr != "mTL") %>%
  group_by(sample) %>%
  summarise(n_total = n(),
            n_T_longer = sum(`T-Sib` > `U-Sib`, na.rm = TRUE),
            n_U_longer = sum(`T-Sib` < `U-Sib`, na.rm = TRUE),
            n_Equal = sum(`T-Sib` == `U-Sib`, na.rm = TRUE),
            tidy(binom.test(n_T_longer, n_total, p = 0.5, alternative = "greater")), .groups = 'drop')

direction_sample_by_Sib <- direction_sample_by_Sib %>% 
  mutate(family_role = str_extract(sample, "[f/m][a/o]"),
         family_role = factor(family_role, levels = c("fa","mo"), labels = c("Father", "Mother"))) %>%
  select(family_role,sample,n_total,n_T_longer,n_U_longer,n_Equal,estimate,conf.low,conf.high,statistic,p.value,method,alternative) %>%
  rename(proportion_T_longer = estimate)

write_xlsx(direction_sample_by_Sib, "direction_sample_by_Sib.xlsx")


direction_sample_by_ASD <- tlens_T_U_ASD_diff %>%
  filter(chr != "mTL") %>%
  group_by(sample) %>%
  summarise(n_total = n(),
            n_T_longer = sum(`T-ASD` > `U-ASD`, na.rm = TRUE),
            n_U_longer = sum(`T-ASD` < `U-ASD`, na.rm = TRUE),
            n_Equal = sum(`T-ASD` == `U-ASD`, na.rm = TRUE),
            tidy(binom.test(n_T_longer, n_total, p = 0.5, alternative = "greater")), .groups = 'drop')

direction_sample_by_ASD <- direction_sample_by_ASD %>% 
  mutate(family_role = str_extract(sample, "[f/m][a/o]"),
         family_role = factor(family_role, levels = c("fa","mo"), labels = c("Father", "Mother"))) %>%
  select(family_role,sample,n_total,n_T_longer,n_U_longer,n_Equal,estimate,conf.low,conf.high,statistic,p.value,method,alternative) %>%
  rename(proportion_T_longer = estimate)

write_xlsx(direction_sample_by_ASD, "direction_sample_by_ASD.xlsx")


##################### Visualization of the direction of transmitted and untransmitted differences (by individual) #####################
direction_sample_by_Sib_plot <- direction_sample_by_Sib %>%
  select(family_role, sample, n_T_longer, n_U_longer) %>%
  pivot_longer(cols = c(n_T_longer, n_U_longer),
               names_to = "direction",
               values_to = "counts") %>%
  mutate(direction = factor(direction, levels = c("n_U_longer","n_T_longer")))

p9 <- ggplot(data = direction_sample_by_Sib_plot, aes(x = sample, y = counts, fill = direction))+
  geom_bar(stat = "identity", position = "fill")+
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.2, color = "black")+
  facet_wrap2(family_role ~ ., nrow = 2, scales = "free_x", strip.position = "right",
              strip = strip_themed(background_y = elem_list_rect(fill = c("#FCD28c80","#A1B5A380"))))+
  scale_fill_manual(values = c("grey80","#EB6368"), labels = c("T-Sib < U-Sib","T-Sib > U-Sib"))+
  scale_y_continuous(expand = c(0,0.03))+
  labs(x = NULL, y = "ChArmTL proportion", fill = NULL)+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 5, hjust = 1, angle = 30),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        legend.position = "top",
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.35),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.25, l = 0.5, r = 0.5, unit = "cm"))


direction_sample_by_ASD_plot <- direction_sample_by_ASD %>%
  select(family_role, sample, n_T_longer, n_U_longer) %>%
  pivot_longer(cols = c(n_T_longer, n_U_longer),
               names_to = "direction",
               values_to = "counts") %>%
  mutate(direction = factor(direction, levels = c("n_U_longer","n_T_longer")))

p10 <- ggplot(data = direction_sample_by_ASD_plot, aes(x = sample, y = counts, fill = direction))+
  geom_bar(stat = "identity", position = "fill")+
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.2, color = "black")+
  facet_wrap2(family_role ~ ., nrow = 2, scales = "free_x", strip.position = "right",
              strip = strip_themed(background_y = elem_list_rect(fill = c("#F9D5B280","#87BEA780"))))+
  scale_fill_manual(values = c("grey80","#EB6368"), labels = c("T-ASD < U-ASD","T-ASD > U-ASD"))+
  scale_y_continuous(expand = c(0,0.03))+
  labs(x = NULL, y = "ChArmTL proportion", fill = NULL)+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 5, hjust = 1, angle = 30),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        legend.position = "top",
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.35),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))

p9_10 <- plot_grid(p9, p10, labels = c('a','b'), label_size = 8, ncol = 1, rel_widths = c(1,1), align = 'h', axis = 'tb')

ggsave(filename = "Supplementary Figure 12.pdf", plot = p9_10, width = 15, height = 15, units = "cm", dpi = 600)


############################## ChArmTL correlation between parents and offspring #############################
######## Build inheritance_pairs
tlens_parent_T <- tlens_hap %>%
  filter(family_role %in% c("Father","Mother"), haplotype %in% c("T-Sib","T-ASD")) %>%
  mutate(parent_role = family_role,
         child_role = ifelse(haplotype == "T-Sib","Sibling","ASD"),
         parent_TL = TL_p75,
         parent_age = age,
         parent_sex = sex,
         parent_haplotype = haplotype) %>%
  select(family, chr, parent_role, child_role, parent_TL, parent_age, parent_sex, parent_haplotype)

tlens_offspring_inherited <- tlens_hap %>%
  filter(family_role %in% c("Sibling","ASD"), haplotype %in% c("Pat","Mat")) %>%
  mutate(child_role = family_role,
         parent_role = ifelse(haplotype == "Pat","Father","Mother"),
         child_TL = TL_p75,
         child_age = age,
         child_sex = sex,
         child_haplotype = haplotype) %>%
  select(family, chr, child_role, parent_role, child_TL, child_age, child_sex, child_haplotype)

#### Merge them into a parent-offspring matching table
inheritance_pairs <- tlens_parent_T %>%
  inner_join(tlens_offspring_inherited, by = c("family","chr","parent_role","child_role")) %>%
  select(family, child_role, parent_role, chr, parent_TL, child_TL, parent_age, child_age,
         child_sex, parent_sex, parent_haplotype, child_haplotype)

write_xlsx(inheritance_pairs, "inheritance_pairs.xlsx")


#### dat_all (including mTL) and dat_arm (only arms)
inheritance_pairs_all <- inheritance_pairs %>%
  mutate(child_age_c = as.numeric(scale(child_age, center = TRUE, scale = FALSE)),
         parent_age_c = as.numeric(scale(parent_age, center = TRUE, scale = FALSE)),
         parent_role = factor(parent_role, levels = c("Father","Mother")),
         child_role = factor(child_role, levels = c("Sibling","ASD")),
         child_sex = factor(child_sex, levels = c("Male","Female")),
         family = factor(family),
         chr = factor(chr))

inheritance_pairs_chr <- inheritance_pairs_all %>%
  filter(chr != "mTL") %>%
  mutate(chr = droplevels(chr),
         parent_TL_c = as.numeric(scale(parent_TL, center = TRUE, scale = FALSE)))


######## main model（arm-level，dat_arm）
m_final <- lmer(child_TL ~ parent_TL_c * parent_role + child_role + child_age_c + parent_age_c + child_sex + (1 | family) + (1 | chr),
                data = inheritance_pairs_chr, REML = TRUE)

summary(m_final)
anova(m_final)

m_final_fixed_eff <- broom.mixed::tidy(m_final, effects = "fixed", conf.int = TRUE, conf.method = "Wald")
write_xlsx(m_final_fixed_eff, "model_final_fixed_effect.xlsx")

#### slope: paternal vs. maternal
m_final_emm_trends <- emtrends(m_final, pairwise ~ parent_role, var = "parent_TL_c", lmer.df  = "satterthwaite", lmerTest.limit = 4000)

summary(m_final_emm_trends$emtrends, infer = c(TRUE, TRUE))    # Check the slope estimates and 95% confidence intervals for father and mother
summary(m_final_emm_trends$contrasts, infer = c(TRUE, TRUE))    # Check the significance test for the slope difference between paternal and maternal (corresponding to the interaction terms)


######## Visualization
#### Calculate the Marginal Effects of the model
pred_data <- ggpredict(m_final, terms = c("parent_TL_c [all]","parent_role")) 
pred_data$x_raw <- pred_data$x + mean(inheritance_pairs_chr$parent_TL, na.rm = TRUE)    # true scale

#### fitting line + 95% confidence interval + background scatter of the original data
p11 <- ggplot() +    
  geom_point(data = inheritance_pairs_chr, aes(x = parent_TL, y = child_TL, color = parent_role),
             alpha = 0.75, size = 0.1, show.legend = FALSE) +
  geom_ribbon(data = pred_data, aes(x = x_raw, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25, key_glyph = "rect") +
  geom_line(data = pred_data, aes(x = x_raw, y = predicted, color = group), linewidth = 0.35, key_glyph = "path") +
  scale_color_manual(name = NULL, values = c("Father" = "#FCD28c", "Mother" = "#A1B5A3"),
                     labels = c("Father: beta = 0.082, P = 0.0004", "Mother: beta = 0.236, P < 0.0001")) +
  scale_fill_manual(name = NULL, values = c("Father" = "#FCD28c", "Mother" = "#A1B5A3"),
                    labels = c("Father: beta = 0.082, P = 0.0004", "Mother: beta = 0.236, P < 0.0001")) +
  scale_x_continuous(expand = c(0,0), breaks = c(2500,5000,7500,10000), labels = c("2,500","5,000","7,500","10,000")) +
  scale_y_continuous(expand = c(0,0), breaks = c(2500,5000,7500,10000), labels = c("2,500","5,000","7,500","10,000")) +
  labs(x = "Parental transmitted ChArmTLs (bp)", y = "Children's inherited ChArmTLs (bp)") +
  annotate("text", x = min(inheritance_pairs_chr$parent_TL, na.rm = TRUE) * 0.85, y = max(inheritance_pairs_chr$child_TL, na.rm = TRUE) * 0.95,
           label = "Interaction (Mother vs. Father):\ndelta beta = 0.154, P < 0.0001", hjust = 0, vjust = 1, size = 2) +
  theme_classic() +
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1), 
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.25, r = 0.5, unit = "cm"))

p3_4 <- plot_grid(p3, p4, labels = c('a','b'), label_size = 8, ncol = 1, rel_heights = c(1,1), align = 'h', axis = 'tb')
p5_6_11 <- plot_grid(p5, p6, p11, labels = c('c','d','e'), label_size = 8, nrow = 1, rel_widths = c(1,1,1), align = 'h', axis = 'tb')

p3_4_5_6_11 <- plot_grid(p3_4, p5_6_11, ncol = 1, rel_heights = c(1,0.75), align = 'v', axis = 'lr')
ggsave(filename = "Figure 4.pdf", plot = p3_4_5_6_11, width = 17, height = 14, units = "cm", dpi = 600)

