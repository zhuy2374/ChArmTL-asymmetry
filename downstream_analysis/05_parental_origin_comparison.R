############################ set the work directory ############################
rm(list=ls())
setwd("Path/to/your_workdir")
getwd()


############################ Load R package and data ###########################
# Define the list of required packages
required_packages <- c("readxl", "dplyr","writexl","tidyverse","ggplot2",
                       "broom","ggrepel","ggh4x","cowplot")

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
library(tidyverse)
library(ggplot2)
library(broom)
library(ggrepel)
library(ggh4x)
library(cowplot)


tlens <- read_xlsx("tlens_chr_filtered.xlsx")
clinical_info <- read_xlsx("31family_info.xlsx", range = "A1:G125")    #clinical phenotype data: family_ID, sample_ID, family_role, sex, age


################## Combine clinical information with ChArmTLs ##################
######## Calculate the haplotype-specific mTL of offspring and extract the ChArmTLs
tlens_offspring_hap <- tlens %>%
  filter(family_role %in% c("ASD","Sibling"), parent_of_origin %in% c("Pat","Mat")) %>%
  select(sample, parent_of_origin, chr, TL_p75) %>%
  bind_rows(tlens %>%
              filter(family_role %in% c("ASD","Sibling"), parent_of_origin %in% c("Pat","Mat")) %>%
              group_by(sample, parent_of_origin) %>%
              summarise(TL_p75 = mean(TL_p75, na.rm = TRUE), .groups = "drop") %>%
              mutate(chr = "mTL") %>%
              select(sample, parent_of_origin, chr, TL_p75)) %>%
  mutate(parent_of_origin = factor(parent_of_origin, levels = c("Pat","Mat")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%
  arrange(sample,parent_of_origin,chr)


########merge data
tlens_offspring_hap <- merge(tlens_offspring_hap, clinical_info, by = "sample") 
tlens_offspring_hap$family_role <- factor(tlens_offspring_hap$family_role, levels = c("Sibling","ASD")) 
write_xlsx(tlens_offspring_hap, "tlens_offspring_hap.xlsx") 


####################### The distribution of Δ(Pat - Mat) #######################
tlens_Pat_Mat_diff <- tlens_offspring_hap %>% 
  pivot_wider(names_from = parent_of_origin, values_from = TL_p75) %>%
  drop_na(Pat, Mat) %>%
  mutate(diff = Pat-Mat)

tlens_Pat_Mat_diff <- tlens_Pat_Mat_diff %>%
  mutate(family_role = "Offspring") %>%
  bind_rows(tlens_Pat_Mat_diff) %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         family_role = factor(family_role, levels = c("Offspring","Sibling","ASD")))

Pat_Mat_diff_stats <- tlens_Pat_Mat_diff %>%
  group_by(family_role, chr) %>%
  summarise(mean_diff = mean(diff, na.rm = TRUE),
            sd_diff = sd(diff, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")))

p1 <- ggplot()+
  geom_errorbar(data = Pat_Mat_diff_stats, aes(x = chr_number, ymin = mean_diff - sd_diff, ymax = mean_diff + sd_diff, color = family_role),
                width = 0.4, linewidth = 0.25, position = position_dodge(width = 0.75))+
  stat_summary(data = tlens_Pat_Mat_diff, aes(x = chr_number, y = diff, color = family_role), 
               geom = "point", fun = "mean", size = 0.35, position = position_dodge(width = 0.75))+
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "black")+
  scale_color_manual(values = c("Offspring" = "grey80", "Sibling" = "#82ADD0", "ASD" = "#EB6368"), name = NULL) +
  scale_y_continuous(breaks = c(-2000,0,2000,4000), labels = c("-2,000","0","2,000","4,000")) +
  scale_x_discrete(limits = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                              "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q", "mTL"))+
  labs(x = NULL, y = "Pat - Mat\nChArmTL (bp)") +
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


######### Analysis of differences between parental and maternal origins ########
######## A custom function compares the telomere lengths of the paternal and maternal origins
pat_mat_difference <- function(data) {
  data %>%
    group_by(family_role, chr) %>%
    summarise(n_pairs = n(),
              mean_Pat = mean(Pat), sd_Pat = sd(Pat),
              mean_Mat = mean(Mat), sd_Mat = sd(Mat),
              tidy(t.test(Pat, Mat, paired = TRUE)), .groups = 'drop')
}


######## Comparison of differences between paternal and maternal origins
diff_pat_mat <- pat_mat_difference(tlens_Pat_Mat_diff)

#### Merge data and add corrected P-values
diff_pat_mat <- diff_pat_mat %>%
  group_by(family_role) %>%
  mutate(p_adj_fdr = case_when(chr == "mTL" ~ p.value,
                               TRUE ~ p.adjust(p.value, method = "fdr")),
         family_role= factor(family_role, levels = c("ASD","Sibling","Offspring")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%
  select(family_role,chr,n_pairs,mean_Pat,sd_Pat,mean_Mat,sd_Mat,estimate,conf.low,conf.high,statistic,p.value,p_adj_fdr,method,alternative) %>%
  rename(mean_diff = estimate) %>%
  arrange(family_role, chr)

write_xlsx(diff_pat_mat, "diff_pat_mat.xlsx")


########### Visualization of parental and maternal origin differences ##########
diff_pat_mat_plot <- diff_pat_mat %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         mean_diff_sign = ifelse(mean_diff > 0, "Pat > Mat","Pat < Mat"),
         mean_diff_sign = factor(mean_diff_sign, levels =c("Pat > Mat","Pat < Mat")),
         significance = case_when(p_adj_fdr < 0.05 ~ "FDR < 0.05",
                                  p.value < 0.05 ~ "Nominal P < 0.05",
                                  TRUE ~ "Not significant"),
         significance = factor(significance, levels = c("FDR < 0.05","Nominal P < 0.05","Not significant")))

p2 <- ggplot(diff_pat_mat_plot, aes(x = chr_number, y = family_role, color = significance, shape = mean_diff_sign, size = abs(mean_diff)))+
  geom_point(stroke = 0.35)+
  scale_color_manual(values = c("FDR < 0.05" = "#EB6368", "Nominal P < 0.05" = "#FCD28c", "Not significant"  = "grey80"), name = NULL)+
  scale_shape_manual(values = c("Pat > Mat" = 16, "Pat < Mat" = 1), name = NULL)+
  scale_size_continuous(range = c(0.5,3), breaks = c(500,1000,1500,2000), name = "|Mean difference|")+
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
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.5, r = 0.25, unit = "cm"))


######## count of significant  arms
Significant_arms <- diff_pat_mat_plot %>%
  filter(!significance == "Not significant", chr != "mTL") %>%
  count(family_role, significance, name = "counts")

p3 <- ggplot(Significant_arms, aes(x = family_role, y = counts, fill = significance))+
  geom_col(width = 0.5, show.legend = FALSE)+
  scale_fill_manual(values = c("FDR < 0.05" = "#EB6368", "Nominal P < 0.05" = "#FCD28c"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,25), breaks = c(0,10,20))+
  coord_flip()+
  labs(x = NULL, y = "Significant arms")+
  geom_text(aes(label = counts, group = significance), position = position_stack(vjust = 0.5), hjust = 2, size  = 1.5) +
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks.x =  element_line(linewidth = 0.2),
        axis.ticks.length.x = unit(0.05,"cm"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 7),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.25, r = 0.5, unit = "cm"))


########## The correlation between Sibling and ASD probands in parental origins differences ##########
diff_corr_sib_asd <- diff_pat_mat %>%
  select(chr,family_role,mean_diff) %>%
  filter(!family_role == "Offspring") %>%
  pivot_wider(names_from = family_role, values_from = mean_diff) %>% 
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]"))) %>%
  filter( chr != "mTL") 

shapiro.test(diff_corr_sib_asd$ASD)
shapiro.test(diff_corr_sib_asd$Sibling)

cor.test(diff_corr_sib_asd$ASD, diff_corr_sib_asd$Sibling, alternative = "two.sided", method = "pearson", exact = TRUE, conf.level = 0.95)

p4 <- ggplot(diff_corr_sib_asd, aes(x = Sibling, y = ASD, label = chr_number))+
  geom_point(size = 0.5, color = "#82ADD0")+
  geom_text_repel(size = 1.5, max.overlaps = Inf, box.padding = 0.05, point.padding = 0.05, segment.color = "grey40")+
  geom_smooth(method = 'lm', formula = y ~ x, se = T, linewidth = 0.35, color = "grey40")+
  geom_rug(outside = F, length = unit(0.02, "npc"), linewidth = 0.2, color = "#82ADD0")+
  scale_x_continuous(expand = c(0,0), limits = c(-200,2700), breaks = c(0,500,1000,1500,2000,2500), labels = c("0","500","1,000","1,500","2,000","2,500"))+
  scale_y_continuous(expand = c(0,0), limits = c(-400,2500), breaks = c(0,500,1000,1500,2000,2500), labels = c("0","500","1,000","1,500","2,000","2,500"))+
  coord_fixed(ratio = 1)+
  labs(x = "Pat - Mat ChArmTLs in sibling (bp)", y = "Pat - Mat ChArmTLs in ASD (bp)")+
  annotate(geom = "text", x = -100, y = 2300, label = "r = 0.42, P = 0.0092", size = 2, hjust = 0)+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.5, r = 0.25, unit = "cm"))


########## Visualization of external queue validation for parental and maternal differences ##########
diff_pat_mat_G4 <- read_xlsx("diff_pat_mat_G4.xlsx")

diff_pat_mat_G4_plot <- diff_pat_mat_G4 %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         mean_diff_sign = ifelse(mean_diff > 0, "Pat > Mat", "Pat < Mat"),
         mean_diff_sign = factor(mean_diff_sign, levels =c("Pat > Mat", "Pat < Mat")),
         significance = case_when(p.value < 0.05 ~ "Nominal P < 0.05",
                                  TRUE ~ "Not significant"),
         significance = factor(significance, levels = c("Nominal P < 0.05", "Not significant")),
         logP = ifelse(p.value < 0.05, -log10(p.value), 1))

p5 <- ggplot(diff_pat_mat_G4_plot, aes(x = chr_number, y = abs(mean_diff), color = mean_diff_sign, shape = significance, size = logP))+
  geom_point(stroke = 0.35)+
  scale_color_manual(values = c("#EB6368","grey80"))+
  scale_shape_manual(values = c("Nominal P < 0.05" = 16, "Not significant" = 1))+
  scale_x_discrete(limits = c("1q","4p","5p","5q","6p","7p","9p","10p","10q","12p",
                              "16p","16q","17q","18p","19p","20p","21q","22q","mTL"))+
  scale_y_continuous(limits = c(0,4500), breaks = c(0,1000,2000,3000,4000), labels = c("0","1,000","2,000","3,000","4,000"))+
  scale_size_continuous(range = c(0.5,3), breaks = c(1, -log10(0.05), 2), labels = c("ns","0.05","0.01"))+
  labs(x = NULL, y = "Mean difference (bp)", size = "P", color = NULL)+
  guides(color = guide_legend(override.aes = list(size = 1.5), position = "top", order = 1), 
         size = guide_legend(override.aes = list(shape = c(1,16,16), stroke = 0.35), position = "top", order = 2),
         shape = "none")+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.25, r = 0.5, unit = "cm"))


#### The direction of paternal and maternal differences (by chromosome arms) ####
######## Define the directional statistical function
diff_direction_chr <- function(data) {
  data %>%
    group_by(family_role, chr) %>%
    summarise(
      n_total = n(),
      n_Pat_longer = sum(Pat > Mat, na.rm = TRUE),
      n_Mat_longer = sum(Pat < Mat, na.rm = TRUE),
      n_Equal = sum(Pat == Mat, na.rm = TRUE),
      tidy(binom.test(n_Pat_longer, n_total, p = 0.5, alternative = "greater")), .groups = 'drop')
}


######## Directional statistics of paternal and maternal differences
direction_chr <- diff_direction_chr(tlens_Pat_Mat_diff)

direction_chr <- direction_chr %>%
  mutate(family_role = factor(family_role, levels = c("Offspring","Sibling","ASD")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%   
  select(family_role,chr,n_total,n_Pat_longer,n_Mat_longer,n_Equal,estimate,conf.low,conf.high,statistic,p.value,method,alternative) %>%
  rename(proportion_Pat_longer = estimate) %>%
  arrange(family_role, chr)

write_xlsx(direction_chr, "direction_chr.xlsx")


########## Visualization of the direction of paternal and maternal differences (by chromosome arms) ##########
direction_chr_plot <- direction_chr %>%
  select(family_role, chr, n_Pat_longer, n_Mat_longer) %>%
  pivot_longer(cols = c(n_Pat_longer, n_Mat_longer),
               names_to = "direction",
               values_to = "counts") %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         direction = factor(direction, levels = c("n_Pat_longer","n_Mat_longer")),
         family_role= factor(family_role, levels = c("Offspring","Sibling","ASD")))

p6 <- ggplot(direction_chr_plot, aes(x = chr_number, y = counts, color = direction))+
  geom_linerange(aes(ymin = 0, ymax = counts), linewidth = 0.35, position = position_dodge(width = 0.5))+
  geom_point(size = 1, position = position_dodge(width = 0.5))+
  facet_wrap2(family_role ~ ., nrow = 3, scales = "free_y", strip.position = "right",
              strip = strip_themed(background_y = elem_list_rect(fill = c("grey80","#82ADD080","#EB636880"))))+
  facetted_pos_scales(y = list(family_role == "Offspring" ~ scale_y_continuous(expand = c(0,0), limits = c(0,60), breaks = c(20,40,60)),
                               family_role == "Sibling" ~ scale_y_continuous(expand = c(0,0), limits = c(0,30), breaks = c(10,20,30)),
                               family_role == "ASD"     ~ scale_y_continuous(expand = c(0,0), limits = c(0,30), breaks = c(10,20,30))))+
  scale_color_manual(values = c("#EB6368","grey80"), labels = c("Pat > Mat","Pat < Mat"))+
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

p1_NULL <- plot_grid(p1, NULL, labels = c('a',''), label_size = 8, nrow = 1, rel_widths = c(1,0.01))
p2_3 <- plot_grid(p2, p3, labels = c('b',''), label_size = 8, nrow = 1, rel_widths = c(1,0.18), align = 'h', axis = 'tb')
p6_NULL <- plot_grid(p6, NULL, labels = c('c',''), label_size = 8, nrow = 1, rel_widths = c(1,0.025), align = 'h', axis = 'tb')
p4_5 <- plot_grid(p4, p5, labels = c('d','e'), label_size = 8, nrow = 1, rel_widths = c(0.75,1), align = 'h')

p1_2_3_4_5_6 <- plot_grid(p1_NULL, p2_3, p6_NULL, p4_5, ncol = 1, align = 'v', axis = 'lr', rel_heights = c(1.2,1,2,1.75))
ggsave(filename = "Figure 3.pdf", plot = p1_2_3_4_5_6, width = 17, height = 20, units = "cm", dpi = 600)


######################### The direction of paternal and maternal differences (by individual) #########################
direction_sample <- tlens_Pat_Mat_diff %>%
  filter(chr != "mTL") %>%
  group_by(sample) %>%
  summarise(n_total = n(),
            n_Pat_longer = sum(Pat > Mat, na.rm = TRUE),
            n_Mat_longer = sum(Pat < Mat, na.rm = TRUE),
            n_Equal = sum(Pat == Mat, na.rm = TRUE),
            tidy(binom.test(n_Pat_longer, n_total, p = 0.5, alternative = "greater")), .groups = 'drop')

direction_sample <- direction_sample %>% 
  mutate(family_role = str_extract(sample, "[p/s]1"),
         family_role = factor(family_role, levels = c("s1","p1"), labels = c("Sibling","ASD"))) %>%
  select(family_role,sample,n_total,n_Pat_longer,n_Mat_longer,n_Equal,estimate,conf.low,conf.high,statistic,p.value,method,alternative) %>%
  rename(proportion_Pat_longer = estimate)

write_xlsx(direction_sample, "direction_sample.xlsx")


########## Visualization of the direction of paternal and maternal differences (by individual) ##########
direction_sample_plot <- direction_sample %>%
  select(family_role, sample, n_Pat_longer, n_Mat_longer) %>%
  pivot_longer(cols = c(n_Pat_longer, n_Mat_longer),
               names_to = "direction",
               values_to = "counts") %>%
  mutate(direction = factor(direction, levels = c("n_Mat_longer","n_Pat_longer")))

p7 <- ggplot(data = direction_sample_plot, aes(x = sample, y = counts, fill = direction))+
  geom_bar(stat = "identity", position = "fill")+
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.2, color = "black")+
  facet_wrap2(family_role ~ ., nrow = 2, scales = "free_x", strip.position = "right",
              strip = strip_themed(background_y = elem_list_rect(fill = c("#82ADD080","#EB636880"))))+
  scale_fill_manual(values = c("grey80","#EB6368"), labels = c("Pat < Mat","Pat > Mat"))+
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
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))

ggsave(filename = "Supplementary Figure 8.pdf", plot = p7, width = 15, height = 8, units = "cm", dpi = 600)

