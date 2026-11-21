########################## set the work directory ##########################
rm(list=ls())
setwd("Path/to/your_workdir")
getwd()

############################ Load R package and data ###########################
cran_packages <- c(
  "readxl", "ggplot2", "dplyr", "RColorBrewer", "writexl", "lme4", "lmerTest", 
  "tidytext", "ggh4x", "tidyr", "circlize", "ComplexHeatmap", "ggplotify", 
  "Hmisc", "ggcorrplot", "ggrepel", "cowplot"
  )

# Only install packages that are not already installed
new_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Install Bioconductor packages (for ComplexHeatmap)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install ComplexHeatmap if it's not already installed
if (!require("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap")
}

library(readxl)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(writexl)
library(lme4)
library(lmerTest)
library(tidytext)
library(ggh4x)
library(tidyr)
library(circlize)
library(ComplexHeatmap)
library(ggplotify)
library(Hmisc)
library(ggcorrplot)
library(ggrepel)
library(cowplot)

tlens <- read_xlsx("tlens_chr_filtered.xlsx")


############################# Length heterogeneity #############################
######## ChArmTLs
#### plot by chr
tlens_by_chr_plot <- tlens %>%
  mutate(chr_number = substr(chr,4,nchar(as.character(chr))),
         chr_number = factor(chr_number, levels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q",
                                                    "6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                                                    "11p","11q","12p","12q","13q","14q","16p","16q","17p",
                                                    "17q","18p","18q","19p","19q","20p","20q","21q","22q")))

p1 <- ggplot(data = tlens_by_chr_plot, aes(x = chr_number, y = TL_p75))+
  geom_point(aes(color = chr_number), position = position_jitter(0.3), size = 0.1, alpha = 0.5, show.legend = F)+
  geom_boxplot(color = "black", fill = NA, linewidth = 0.3, median.linewidth = 0.3, width = 0.5, outlier.alpha = 0, show.legend = FALSE)+
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Paired"))(length(unique(tlens_by_chr_plot$chr_number)))) +
  scale_y_continuous(expand = c(0,0) ,limits = c(0,13000), breaks = c(0,4000,8000,12000), labels =c("0","4,000","8,000","12,000"))+
  labs(x = NULL, y = "ChArmTLs (bp)", color = NULL)+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.25, l = 0.5, r = 0.5, unit = "cm"))

#### Statistics of ChArmTL
tlens_by_chr_stats <- tlens_by_chr_plot %>%
  group_by(chr) %>%
  summarise(N = sum(!is.na(TL_p75)),
            Mean = mean(TL_p75, na.rm = TRUE),
            SD = sd(TL_p75, na.rm = TRUE),
            Median = median(TL_p75, na.rm = TRUE),
            Q1 = quantile(TL_p75, 0.25, na.rm = TRUE),
            Q3 = quantile(TL_p75, 0.75, na.rm = TRUE),
            IQR = IQR(TL_p75, na.rm = TRUE),
            Min = min(TL_p75, na.rm = TRUE),
            Max = max(TL_p75, na.rm = TRUE),
            CV = sd(TL_p75, na.rm = TRUE) / mean(TL_p75, na.rm = TRUE),
            .groups = "drop") %>%
  bind_rows(tlens %>%
              summarise(N = sum(!is.na(TL_p75)),
                        Mean = mean(TL_p75, na.rm = TRUE),
                        SD = sd(TL_p75, na.rm = TRUE),
                        Median = median(TL_p75, na.rm = TRUE),
                        Q1 = quantile(TL_p75, 0.25, na.rm = TRUE),
                        Q3 = quantile(TL_p75, 0.75, na.rm = TRUE),
                        IQR = IQR(TL_p75, na.rm = TRUE),
                        Min = min(TL_p75, na.rm = TRUE),
                        Max = max(TL_p75, na.rm = TRUE),
                        CV = sd(TL_p75, na.rm = TRUE) / mean(TL_p75, na.rm = TRUE))%>% 
              mutate(chr = "All", .before = 1))

write_xlsx(tlens_by_chr_stats, "tlens_by_chr_stats.xlsx")

#### Inter-arm heterogeneity of ChArmTLs
## Global test of mean difference
tlens_by_chr_plot <- tlens_by_chr_plot %>%
  mutate(family = substr(sample,1,nchar(as.character(sample))-3))    #Add the family ID column

model_full <- lmer(TL_p75 ~ chr + (1 | sample), data = tlens, REML = FALSE)
model_null <- lmer(TL_p75 ~ 1 + (1 | sample), data = tlens, REML = FALSE)

anova(model_null,model_full)    # The likelihood ratio test (LRT) shows that there is a significant inter-arm heterogeneity in ChArmTLs


######## ChArmTL allelic differences
allele_diff <- tlens %>%
  group_by(sample, chr) %>%
  summarise(n_allele = n(),
            len_1 = TL_p75[1],
            len_2 = TL_p75[2],
            allele_diff = abs(TL_p75[1] - TL_p75[2]),
            .groups = "drop") %>%
  filter(n_allele == 2)


#### plot by chr
allele_diff_by_chr_plot <- allele_diff %>%
  mutate(chr_number = substr(chr,4,nchar(as.character(chr))),
         chr_number = factor(chr_number, levels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q",
                                                    "6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                                                    "11p","11q","12p","12q","13q","14q","16p","16q","17p",
                                                    "17q","18p","18q","19p","19q","20p","20q","21q","22q")))

p2 <- ggplot(data = allele_diff_by_chr_plot, aes(x = chr_number, y = allele_diff))+
  geom_point(aes(color = chr_number), position = position_jitter(0.3), size = 0.1, alpha = 0.5, show.legend = F)+
  geom_boxplot(color = "black", fill = NA, linewidth = 0.3, median.linewidth = 0.3, width = 0.5, outlier.alpha = 0, show.legend = FALSE)+
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Paired"))(length(unique(allele_diff_by_chr_plot$chr_number)))) +
  scale_y_continuous(expand = c(0,0) ,limits = c(0,9050), breaks = c(0,3000,6000,9000), labels =c("0","3,000","6,000","9,000"))+
  labs(x = NULL, y = "ChArmTL allelic differences (bp)", color = NULL)+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.25, l = 0.5, r = 0.5, unit = "cm"))

#### Statistics of ChArmTL allelic differences
allele_diff_by_chr_stats <- allele_diff_by_chr_plot %>%
  group_by(chr) %>%
  summarise(N = sum(!is.na(allele_diff)),
            Mean = mean(allele_diff, na.rm = TRUE),
            SD = sd(allele_diff, na.rm = TRUE),
            Median = median(allele_diff, na.rm = TRUE),
            Q1 = quantile(allele_diff, 0.25, na.rm = TRUE),
            Q3 = quantile(allele_diff, 0.75, na.rm = TRUE),
            IQR = IQR(allele_diff, na.rm = TRUE),
            Min = min(allele_diff, na.rm = TRUE),
            Max = max(allele_diff, na.rm = TRUE),
            CV = sd(allele_diff, na.rm = TRUE) / mean(allele_diff, na.rm = TRUE), .groups = "drop") %>%
  bind_rows(allele_diff %>%
              summarise(N = sum(!is.na(allele_diff)),
                        Mean = mean(allele_diff, na.rm = TRUE),
                        SD = sd(allele_diff, na.rm = TRUE),
                        Median = median(allele_diff, na.rm = TRUE),
                        Q1 = quantile(allele_diff, 0.25, na.rm = TRUE),
                        Q3 = quantile(allele_diff, 0.75, na.rm = TRUE),
                        IQR = IQR(allele_diff, na.rm = TRUE),
                        Min = min(allele_diff, na.rm = TRUE),
                        Max = max(allele_diff, na.rm = TRUE),
                        CV = sd(allele_diff, na.rm = TRUE) / mean(allele_diff, na.rm = TRUE))%>% 
              mutate(chr = "All", .before = 1))

write_xlsx(allele_diff_by_chr_stats, "allele_diff_by_chr_stats.xlsx")


################################## rank order ##################################
######## ChArmTLs
#### all sample
tlens_by_sample_stats <- tlens %>%
  group_by(sample) %>%
  summarise(mean_by_sample = mean(TL_p75, na.rm = TRUE),
            sd_by_sample = sd(TL_p75, na.rm = TRUE),
            .groups = "drop")

tlens_z_score <- tlens %>%
  left_join(tlens_by_sample_stats, by = "sample") %>%
  mutate(relative_tlens = TL_p75 - mean_by_sample,
         z_score_tlens = relative_tlens / sd_by_sample) %>%
  select(-mean_by_sample, -sd_by_sample, -relative_tlens)

rank_arm_by_sample <- tlens_z_score %>%
  group_by(family_role, chr) %>%
  summarise(mean_z_score_tlens = mean(z_score_tlens, na.rm = TRUE),
            sd_z_score_tlens = sd(z_score_tlens, na.rm = TRUE),
            .groups = "drop") %>%
  bind_rows(tlens_z_score %>%
              group_by(chr) %>%
              summarise(mean_z_score_tlens = mean(z_score_tlens, na.rm = TRUE),
                        sd_z_score_tlens = sd(z_score_tlens, na.rm = TRUE)) %>%
              mutate(family_role = "All", .before = 1)) %>%
  group_by(family_role) %>%
  mutate(rank_mean_z_score_tlens = rank(mean_z_score_tlens, ties.method = "min"),
         chr_number = substr(chr,4,nchar(as.character(chr))),
         family_role = factor(family_role, levels = c("All","Father","Mother","Sibling","ASD")))

write_xlsx(rank_arm_by_sample, "rank_arm_by_sample.xlsx")

## The overall rank of all samples
tlens_z_score_sample_all <- tlens_z_score %>%
  left_join(rank_arm_by_sample %>%
              ungroup() %>%
              filter(family_role == "All") %>%
              select(-family_role, -chr_number), by = "chr") %>%
  mutate(chr_number = substr(chr,4,nchar(as.character(chr))))

p3 <- ggplot()+
  geom_jitter(data = tlens_z_score_sample_all,
              aes(x = reorder(chr_number, rank_mean_z_score_tlens), y = z_score_tlens, color = chr_number),
              width = 0.25, size = 0.1, alpha = 0.5, show.legend = FALSE)+
  geom_errorbar(data = subset(rank_arm_by_sample, family_role == "All"),
                aes(x = reorder(chr_number, rank_mean_z_score_tlens), ymin = mean_z_score_tlens - sd_z_score_tlens, ymax = mean_z_score_tlens + sd_z_score_tlens),
                width = 0.25, linewidth = 0.35)+
  stat_summary(data = tlens_z_score_sample_all,
               aes(x = reorder(chr_number, rank_mean_z_score_tlens), y = z_score_tlens), 
               geom = "crossbar", fun = "mean", linewidth = 0.15, width = 0.5)+
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "black")+
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Paired"))(length(unique(tlens_z_score_sample_all$chr_number))))+
  labs(x = NULL, y = "ChArmTL Z-scores")+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.5, r = 0.5, unit = "cm"))

## by family_role
tlens_z_score_sample_role <- tlens_z_score %>%
  left_join(rank_arm_by_sample, by = c("family_role","chr")) %>%
  mutate(chr_number = substr(chr,4,nchar(as.character(chr))),
         family_role = factor(family_role, levels = c("All","Father","Mother","Sibling","ASD")))

p4 <- ggplot()+
  geom_jitter(data = tlens_z_score_sample_role,
              aes(x = reorder_within(chr_number, rank_mean_z_score_tlens, family_role), y = z_score_tlens, color = family_role),
              width = 0.25, size = 0.75, show.legend = FALSE)+
  geom_errorbar(data = subset(rank_arm_by_sample, family_role != "All"),
                aes(x = reorder_within(chr_number, rank_mean_z_score_tlens, family_role), ymin = mean_z_score_tlens - sd_z_score_tlens, ymax = mean_z_score_tlens + sd_z_score_tlens),
                width = 0.25, linewidth = 0.35)+
  stat_summary(data = tlens_z_score_sample_role,
               aes(x = reorder_within(chr_number, rank_mean_z_score_tlens, family_role), y = z_score_tlens), 
               geom = "crossbar", fun = "mean", linewidth = 0.15, width = 0.5)+
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "black")+
  facet_wrap2(~ family_role, nrow = 4, scales = "free_x", strip.position = "right",
              strip = strip_themed(background_y = elem_list_rect(fill = c("#FCD28c80","#A1B5A380","#82ADD080","#EB636880"))))+
  scale_color_manual(values = c("#FCD28c","#A1B5A3","#82ADD0","#EB6368"))+
  scale_x_reordered()+
  labs(x = NULL, y = "ChArmTL Z-scores")+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.35),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))

ggsave(filename = "Supplementary Figure 4.pdf", plot = p4, width = 17, height = 15, units = "cm", dpi = 600)

#### by haplotype
tlens_hap <- tlens %>%
  filter(family_role %in% c("Sibling","ASD"), parent_of_origin %in% c("Pat","Mat"))%>%
  mutate(haplotype = case_when(family_role == "Sibling" & parent_of_origin == "Pat" ~ "Sib-Pat",
                               family_role == "Sibling" & parent_of_origin == "Mat" ~ "Sib-Mat",
                               family_role == "ASD" & parent_of_origin == "Pat" ~ "ASD-Pat",
                               family_role == "ASD" & parent_of_origin == "Mat" ~ "ASD-Mat",
                               TRUE ~ NA_character_)) %>%
  bind_rows(tlens %>%
              filter(family_role %in% c("Father","Mother"), trans_status_Sib %in% c("T-Sib","U-Sib"))%>%
              mutate(haplotype = case_when(family_role == "Father" & trans_status_Sib == "T-Sib" ~ "FT-Sib",
                                           family_role == "Father" & trans_status_Sib == "U-Sib" ~ "FU-Sib",
                                           family_role == "Mother" & trans_status_Sib == "T-Sib" ~ "MT-Sib",
                                           family_role == "Mother" & trans_status_Sib == "U-Sib" ~ "MU-Sib",
                                           TRUE ~ NA_character_))) %>%
  bind_rows(tlens %>% 
              filter(family_role %in% c("Father","Mother"), trans_status_ASD %in% c("T-ASD","U-ASD"))%>%
              mutate(haplotype = case_when(family_role == "Father" & trans_status_ASD == "T-ASD" ~ "FT-ASD",
                                           family_role == "Father" & trans_status_ASD == "U-ASD" ~ "FU-ASD",
                                           family_role == "Mother" & trans_status_ASD == "T-ASD" ~ "MT-ASD",
                                           family_role == "Mother" & trans_status_ASD == "U-ASD" ~ "MU-ASD",
                                           TRUE ~ NA_character_))) %>%
  mutate(haplotype = factor(haplotype, levels = c("Sib-Pat","Sib-Mat","ASD-Pat","ASD-Mat",
                                                  "FT-Sib","FU-Sib","MT-Sib","MU-Sib",
                                                  "FT-ASD","FU-ASD","MT-ASD","MU-ASD")))

tlens_by_hap_stats <- tlens_hap %>%
  group_by(sample, haplotype) %>%
  summarise(mean_by_hap = mean(TL_p75, na.rm = TRUE),
            sd_by_hap = sd(TL_p75, na.rm = TRUE),
            .groups = "drop")

tlens_hap_z_score <- tlens_hap %>%
  left_join(tlens_by_hap_stats, by = c("sample","haplotype")) %>%
  mutate(relative_tlens = TL_p75 - mean_by_hap,
         z_score_tlens = relative_tlens / sd_by_hap) %>%
  select(-mean_by_hap, -sd_by_hap, -relative_tlens)

rank_arm_by_hap <- tlens_hap_z_score %>%
  group_by(haplotype, chr) %>%
  summarise(mean_z_score_tlens = mean(z_score_tlens, na.rm = TRUE),
            sd_z_score_tlens = sd(z_score_tlens, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(haplotype) %>%
  mutate(rank_mean_z_score_tlens = rank(mean_z_score_tlens, ties.method = "min"),
         chr_number = substr(chr,4,nchar(as.character(chr))),
         chr_number = factor(chr_number, levels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q",
                                                    "6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                                                    "11p","11q","12p","12q","13q","14q","16p","16q","17p",
                                                    "17q","18p","18q","19p","19q","20p","20q","21q","22q"))) %>%
  arrange(haplotype, chr_number)


write_xlsx(rank_arm_by_hap, "rank_arm_by_hap.xlsx")


######## ChArmTL allelic differences
#### all sample
allele_diff_by_sample_stats <- allele_diff %>%
  group_by(sample) %>%
  summarise(mean_allele_diff = mean(allele_diff, na.rm = TRUE),
            sd_allele_diff = sd(allele_diff, na.rm = TRUE),
            .groups = "drop")

allele_diff_z_score <- allele_diff %>%
  mutate(raw_family_role = tolower(sub(".*_", "", sample)),
         family_role = case_when(raw_family_role == "fa" ~ "Father",
                                 raw_family_role == "mo" ~ "Mother",
                                 raw_family_role == "s1" ~ "Sibling",
                                 raw_family_role == "p1" ~ "ASD",
                                 TRUE ~ raw_family_role)) %>%
  left_join(allele_diff_by_sample_stats, by = "sample") %>%
  mutate(relative_allele_diff = allele_diff - mean_allele_diff,
         z_score_allele_diff = relative_allele_diff / sd_allele_diff) %>%
  select(-raw_family_role, -mean_allele_diff, -sd_allele_diff, -relative_allele_diff)

rank_arm_by_allele_diff <- allele_diff_z_score %>%
  group_by(family_role, chr) %>%
  summarise(mean_z_score_allele_diff = mean(z_score_allele_diff, na.rm = TRUE),
            sd_z_score_allele_diff = sd(z_score_allele_diff, na.rm = TRUE),
            .groups = "drop") %>%
  bind_rows(allele_diff_z_score %>%
              group_by(chr) %>%
              summarise(mean_z_score_allele_diff = mean(z_score_allele_diff, na.rm = TRUE),
                        sd_z_score_allele_diff = sd(z_score_allele_diff, na.rm = TRUE)) %>%
              mutate(family_role = "All", .before = 1)) %>%
  group_by(family_role) %>% 
  mutate(rank_mean_z_score_allele_diff = rank(mean_z_score_allele_diff, ties.method = "min"),
         chr_number = substr(chr,4,nchar(as.character(chr))),
         family_role = factor(family_role, levels = c("All","Father","Mother","Sibling","ASD")))

write_xlsx(rank_arm_by_allele_diff, "rank_arm_by_allele_diff.xlsx")

## The overall rank of all samples
allele_diff_z_score_sample_all <- allele_diff_z_score %>%
  left_join(rank_arm_by_allele_diff %>%
              ungroup() %>%
              filter(family_role == "All") %>%
              select(-family_role, -chr_number), by = "chr") %>%
  mutate(chr_number = substr(chr,4,nchar(as.character(chr))))

p5 <- ggplot()+
  geom_jitter(data = allele_diff_z_score_sample_all,
              aes(x = reorder(chr_number, rank_mean_z_score_allele_diff), y = z_score_allele_diff, color = chr_number),
              width = 0.25, size = 0.1, alpha = 0.5, show.legend = FALSE)+
  geom_errorbar(data = subset(rank_arm_by_allele_diff, family_role == "All"),
                aes(x = reorder(chr_number, rank_mean_z_score_allele_diff), ymin = mean_z_score_allele_diff - sd_z_score_allele_diff, ymax = mean_z_score_allele_diff + sd_z_score_allele_diff),
                width = 0.25, linewidth = 0.35)+
  stat_summary(data = allele_diff_z_score_sample_all,
               aes(x = reorder(chr_number, rank_mean_z_score_allele_diff), y = z_score_allele_diff), 
               geom = "crossbar", fun = "mean", linewidth = 0.15, width = 0.5)+
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "black")+
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Paired"))(length(unique(allele_diff_z_score_sample_all$chr_number))))+
  labs(x = NULL, y = "Z-scores of ChArmTL allelic difference")+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.5, r = 0.5, unit = "cm"))

## by family_role
allele_diff_z_score_sample_role <- allele_diff_z_score %>%
  left_join(rank_arm_by_allele_diff %>%
              select(-chr_number), by = c("family_role","chr")) %>%
  mutate(chr_number = substr(chr,4,nchar(as.character(chr))),
         family_role = factor(family_role, levels = c("All","Father","Mother","Sibling","ASD")))

p6 <- ggplot()+
  geom_jitter(data = allele_diff_z_score_sample_role,
              aes(x = reorder_within(chr_number, rank_mean_z_score_allele_diff, family_role), y = z_score_allele_diff, color = family_role),
              width = 0.25, size = 0.75, show.legend = FALSE)+
  geom_errorbar(data = subset(rank_arm_by_allele_diff, family_role != "All"),
                aes(x = reorder_within(chr_number, rank_mean_z_score_allele_diff, family_role), ymin = mean_z_score_allele_diff - sd_z_score_allele_diff, ymax = mean_z_score_allele_diff + sd_z_score_allele_diff),
                width = 0.25, linewidth = 0.35)+
  stat_summary(data = allele_diff_z_score_sample_role,
               aes(x = reorder_within(chr_number, rank_mean_z_score_allele_diff, family_role), y = z_score_allele_diff), 
               geom = "crossbar", fun = "mean", linewidth = 0.15, width = 0.5)+
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "black")+
  facet_wrap2(~ family_role, nrow = 4, scales = "free_x", strip.position = "right",
              strip = strip_themed(background_y = elem_list_rect(fill = c("#FCD28c80","#A1B5A380","#82ADD080","#EB636880"))))+
  scale_color_manual(values = c("#FCD28c","#A1B5A3","#82ADD0","#EB6368"))+
  scale_x_reordered()+
  labs(x = NULL, y = "Z-scores of ChArmTL allelic difference")+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.35),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))

ggsave(filename = "Supplementary Figure 6.pdf", plot = p6, width = 17, height = 15, units = "cm", dpi = 600)


#################################### heatmap ###################################
######## ChArmTLs
#### Merge data
arm_conservatism_sample <- rank_arm_by_sample %>%
  select(chr_number, family_role, rank_mean_z_score_tlens) %>%
  pivot_wider(names_from = family_role, values_from = rank_mean_z_score_tlens) %>%
  select(chr_number,All,Father,Mother,Sibling,ASD)

arm_conservatism_hap <- rank_arm_by_hap %>%
  select(chr_number, haplotype, rank_mean_z_score_tlens) %>%
  pivot_wider(names_from = haplotype, values_from = rank_mean_z_score_tlens)

arm_conservatism <- arm_conservatism_sample %>%
  full_join(arm_conservatism_hap, by = "chr_number")

## Add two sets of data from Karimian et al.
Karimian_rank_order <- data.frame(chr_number = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q",
                                                 "6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                                                 "11p","11q","12p","12q","13q","14q","16p","16q","17p",
                                                 "17q","18p","18q","19p","19q","20p","20q","21q","22q"),
                                  `Cord_blood` = c(26,NA,32,6,36,13,20,34,31,21,12,29,9,NA,33,10,14,5,28,22,18,24,
                                                   3,35,30,25,19,8,1,15,NA,27,11,23,17,2,7,4),
                                  `Adults_147` = c(31,23,35,8,38,6,25,36,32,12,19,29,10,33,30,18,15,7,22,28,21,20,3,
                                                   37,34,24,16,13,1,4,NA,26,9,27,11,2,14,5))

arm_conservatism <- arm_conservatism %>%
  full_join(Karimian_rank_order, by = "chr_number")

write_xlsx(arm_conservatism, "arm_conservatism.xlsx")

#### heatmap
rank_order_tlens <- as.vector(arm_conservatism[order(arm_conservatism$All),]$chr_number)
#"20q","17p","22q","16p","3q","9p","12p","17q","19p","19q","9q","21q","16q","7p","6p","20p","2q","11p","18p","5q","4p","14q","10p","10q","11q",
#"1p","1q","6q","18q","8p","8q","2p","13q","12q","7q","4q","3p","5p" 

## converted into a matrix
arm_conservatism_matrix <- as.matrix(arm_conservatism[,c(2:20)])    # Convert to matrix format
rownames(arm_conservatism_matrix) <- arm_conservatism$chr_number     # chr_number as rowname
arm_conservatism_matrix_t <- t(arm_conservatism_matrix)
arm_conservatism_matrix_t <- arm_conservatism_matrix_t[,rank_order_tlens]

## Define the colors of the heat map
col_fun_tlens <- colorRamp2(seq(min(arm_conservatism_matrix, na.rm = TRUE), max(arm_conservatism_matrix, na.rm = TRUE), length.out = 7),
                            rev(brewer.pal(7, "RdYlBu")))

## Create row splitting factors (every 4 rows in a group, with the last 2 rows in a group)
row_split_tlens <- factor(rep(c("Group1", "Group2", "Group3", "Group4", "Group5"), times = c(5,4,4,4,2)))

## heat map
p7 <- as.ggplot(Heatmap(arm_conservatism_matrix_t,
                        cluster_rows = FALSE, cluster_columns = FALSE,
                        col = col_fun_tlens, na_col = "grey50",
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        column_names_side = "top",
                        row_names_side = "left",
                        column_names_gp = gpar(fontsize = 5), column_names_rot = 70,
                        row_names_gp = gpar(fontsize = 5),
                        row_split = row_split_tlens,
                        row_title = NULL,
                        heatmap_legend_param = list(
                          title = "Rank",
                          title_position = "topcenter",
                          title_gp = gpar(fontsize = 6),
                          labels_gp = gpar(fontsize = 5),
                          grid_width = unit(2, "mm"))))+
  theme(plot.margin = margin(t = 0.1, b = 0.1, l = 0.1, r = 0.1, unit = "cm"))


######## ChArmTL allelic differences
arm_conservatism_allele_diff <- rank_arm_by_allele_diff %>%
  select(chr_number,family_role, rank_mean_z_score_allele_diff) %>%
  pivot_wider(names_from = family_role, values_from = rank_mean_z_score_allele_diff) %>%
  select(chr_number,All,Father,Mother,Sibling,ASD)

write_xlsx(arm_conservatism_allele_diff, "arm_conservatism_allele_diff.xlsx")

#### heat map
rank_order_allele_diff <- as.vector(arm_conservatism_allele_diff[order(arm_conservatism_allele_diff$All),]$chr_number)
#"9p","20q","19q","11q","17p","7q","16p","2p","12p","22q","4p","6p","5q","3p","3q","14q","8p","17q","10p","21q","18q","18p","2q","13q","16q",
#"7p","5p","19p","1q","1p","12q","8q","10q","4q","6q","9q","11p","20p"

## converted into a matrix
arm_conservatism_allele_diff_matrix <- as.matrix(arm_conservatism_allele_diff[,c(2:6)])    # Convert to matrix format
rownames(arm_conservatism_allele_diff_matrix) <- arm_conservatism_allele_diff$chr_number     # chr_number as rowname
arm_conservatism_allele_diff_matrix_t <- t(arm_conservatism_allele_diff_matrix)
arm_conservatism_allele_diff_matrix_t <- arm_conservatism_allele_diff_matrix_t[,rank_order_allele_diff]

## Define the colors of the heat map
col_fun_allele_diff <- colorRamp2(seq(min(arm_conservatism_allele_diff_matrix_t, na.rm = TRUE),
                                      max(arm_conservatism_allele_diff_matrix_t, na.rm = TRUE),
                                      length.out = 7),
                                  rev(brewer.pal(7, "RdYlBu")))

## heat map
p8 <- as.ggplot(Heatmap(arm_conservatism_allele_diff_matrix_t,
                        cluster_rows = FALSE, cluster_columns = FALSE,
                        col = col_fun_allele_diff, na_col = "grey50",
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        column_names_side = "top",
                        row_names_side = "left",
                        column_names_gp = gpar(fontsize = 5), column_names_rot = 70,
                        row_names_gp = gpar(fontsize = 5),
                        heatmap_legend_param = list(
                          title = "Rank",
                          title_position = "topcenter",
                          title_gp = gpar(fontsize = 6),
                          labels_gp = gpar(fontsize = 5),
                          grid_width = unit(2, "mm"))))+
  theme(plot.margin = margin(t = 0.1, b = 0.1, l = 0.1, r = 0.1, unit = "cm"))


################# The correlation between groups in rank order #################
######## ChArmTLs
arm_conservatism_cor <- rcorr(arm_conservatism_matrix, type = "spearman")    # correlation matrix

write.csv(arm_conservatism_cor$r, "arm_conservatism_cor.csv")
write.csv(arm_conservatism_cor$P, "arm_conservatism_cor_P.csv")

p9 <- ggcorrplot(arm_conservatism_cor$r[19:1,19:1],
                 method = "square",
                 type = "upper",
                 show.diag = TRUE,
                 outline.color = "NA", 
                 lab = TRUE,
                 lab_size = 1,
                 colors = c("#2166ac", "white", "#b2182b"),
                 tl.srt = 45,
                 ggtheme  = theme_minimal())+
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(size = 5, hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        legend.key.width = unit(0.2, "cm"),
        legend.position = c(0.9,0.5),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0.1, b = 0.1, l = 0.1, r = 0.1, unit = "cm"))


#### The correlation between the overall rank order of ChArmTLs and the overall rank order of ChArmTL allelic differences
rank_arm_all <- arm_conservatism %>%
  select(chr_number, All) %>%
  rename(ChArmTL = All) %>%
  full_join(arm_conservatism_allele_diff[,c(1,2)], by = "chr_number") %>%
  rename(allele_diff = All)

## correlation test
cor.test(rank_arm_all$ChArmTL, rank_arm_all$allele_diff, method = "spearman", alternative = "two.sided", exact = FALSE, conf.level = 0.95)

p10 <- ggplot(rank_arm_all, aes(x = allele_diff, y = ChArmTL, label = chr_number))+
  geom_point(size = 0.5, color = "#82ADD0")+
  geom_text_repel(size = 1.5, max.overlaps = Inf, box.padding = 0.05, point.padding = 0.05, segment.color = "grey40")+
  geom_smooth(method = 'lm', formula = y ~ x, se = T, linewidth = 0.35, color = "grey40")+
  scale_x_continuous(expand = c(0,0), limits = c(0,40), breaks = c(10,20,30,40))+
  scale_y_continuous(expand = c(0,0), limits = c(0,40), breaks = c(10,20,30,40))+
  coord_fixed(ratio = 1)+
  labs(x = "Rank by ChArmTL allelic difference", y = "Rank by ChArmTL")+
  annotate(geom = "text", x = 2, y = 38, label = "rho = 0.37, P = 0.023", size = 2, hjust = 0)+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.25, r = 0.5, unit = "cm"))


#################################### cowplot ###################################
######## Figure 1
p1_3 <- plot_grid(p1, p3, ncol = 1, labels = c('a', 'b'), label_size = 8, rel_widths = c(1,1), align = 'v')
p7_9 <- plot_grid(p7, p9, ncol = 2, labels = c('c', 'd'), label_size = 8, rel_widths = c(1.25,1), align = 'h')

Figure_1 <- plot_grid(p1_3, p7_9, nrow = 2, rel_heights = c(1,1.2))

ggsave(filename = "Figure 1.pdf", plot = Figure_1, width = 17, height = 15, units = "cm", dpi = 600)


######## Supplementary Figure 5
p2_5 <- plot_grid(p2, p5, ncol = 1, labels = c('a', 'b'), label_size = 8, rel_widths = c(1,1), align = 'v')
p8_10 <- plot_grid(p8, p10, ncol = 2, labels = c('c', 'd'), label_size = 8, rel_widths = c(2,1), align = 'v')

Figure_S5 <- plot_grid(p2_5, p8_10, nrow = 2, rel_heights = c(1,0.75))

ggsave(filename = "Supplementary Figure 5.pdf", plot = Figure_S5, width = 13, height = 12, units = "cm", dpi = 600)

