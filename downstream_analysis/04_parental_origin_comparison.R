########################## set the work directory ##########################
rm(list=ls())
setwd("Path/to/your_workdir")
getwd()


############################ Load R package and data ###########################
# Define the list of required packages
required_packages <- c("readxl", "dplyr","tidyverse","broom",
                       "writexl","ggplot2","cowplot","ggrepel","ggh4x")

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
library(tidyverse)    # pivot_wider function
library(broom)    #For tidy extraction of statistical results
library(writexl)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggh4x)    # facet_wrap2, modify facet strip background colors


tlens <- read_xlsx("tlens_chr_filtered_mean.xlsx")


################### Paternal vs maternal difference analysis ###################
########Data reshaping
tlens_offspring <- tlens %>% 
  filter(clinical_group %in% c("ASD", "Sibling"),
         level %in% c("Pat", "Mat")) 

tlens_offspring_wide <- tlens_offspring %>%
  pivot_wider(names_from = level, values_from = TL_p75) %>%
  drop_na(Pat, Mat)   # Keep only rows with both paternal and maternal data


########Custom function to compare paternal and maternal ChArmTLs
pat_mat_difference <- function(data, label) {
  data %>%
    group_by(chr) %>%
    summarise(n_pairs = n(),  # Count the number of paired samples in this group
              mean_Pat = mean(Pat), sd_Pat = sd(Pat),    #Summarize the Pat group
              mean_Mat = mean(Mat), sd_Mat = sd(Mat),    #Summarize the Mat group
              tidy(t.test(Pat, Mat, paired = TRUE)), .groups = 'drop') %>%
    mutate(analysis_level = label)
}


########Paternal vs maternal difference comparison
diff_offspring <- pat_mat_difference(tlens_offspring_wide, "Offspring")    #Offspring
diff_sib <- pat_mat_difference(subset(tlens_offspring_wide, clinical_group == "Sibling"), "Sibling")
diff_asd <- pat_mat_difference(subset(tlens_offspring_wide, clinical_group == "ASD"), "ASD")

####Merge data and add adjusted P-values
pat_mat_diff <- bind_rows(diff_offspring, diff_sib, diff_asd) %>%
  group_by(analysis_level) %>%
  mutate(p_adj_fdr = case_when(
    chr == "mTL" ~ p.value,
    TRUE ~ p.adjust(p.value, method = "fdr")),    #P-value adjustment
    analysis_level= factor(analysis_level, levels = c("Offspring","Sibling", "ASD")),
       chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                    "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                    "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                    "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%   
  select(analysis_level,chr,n_pairs,mean_Pat,sd_Pat,mean_Mat,sd_Mat,estimate,conf.low,conf.high,statistic,p.value,p_adj_fdr,method,alternative) %>%
  rename(mean_diff = estimate) %>%
  arrange(analysis_level, chr)

write_xlsx(pat_mat_diff, "pat_mat_diff.xlsx")


############################## Visualization of paternal vs maternal differences ##############################
pat_mat_diff_plot <- pat_mat_diff %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         mean_diff_sign = ifelse(mean_diff > 0, "Pat > Mat", "Pat < Mat"),
         significance = ifelse(p.value < 0.05, "Significant", "Not Significant"),
         logP = ifelse(p.value < 0.05, -log10(p.value), 1))      #Effect direction

pat_mat_diff_plot$logP[pat_mat_diff_plot$logP > 10] <- pat_mat_diff_plot$logP[pat_mat_diff_plot$logP > 10] - 6    # Map the position of mTL P-values

pat_mat_diff_plot$analysis_level <- factor(pat_mat_diff_plot$analysis_level, levels = c( "ASD","Sibling","Offspring"))
pat_mat_diff_plot$mean_diff_sign <- factor(pat_mat_diff_plot$mean_diff_sign, levels = c("Pat > Mat", "Pat < Mat"))
pat_mat_diff_plot$significance <- factor(pat_mat_diff_plot$significance, levels = c("Significant", "Not Significant"))

p1 <- ggplot(pat_mat_diff_plot, aes(x = chr_number, y = analysis_level, color = mean_diff_sign, shape = significance, size = logP))+
  geom_point(stroke = 0.35)+
  scale_color_manual(values = c("#EB6368","#82ADD0"))+
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 1))+
  theme_bw()+
  scale_x_discrete(limits = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                              "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q", "mTL"))+
  scale_size_continuous(range = c(0.5,3), breaks = c(2,4,6,8))+    #Control the minimum and maximum bubble sizes
  labs(x = NULL, y = NULL, color = NULL, size = "-log10(P)")+
  guides(color = guide_legend(override.aes = list(size = 1.5), position = "top", order = 1), 
         size = guide_legend(override.aes = list(shape = 16, stroke = 0.35), position = "top", order = 2),
         shape = "none")+    #Adjust point size in the legend and the legend position
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
        legend.box.spacing = unit(0.01, "cm"),    #Distance between the legend box and the plot area
        panel.grid = element_line(linewidth = 0.2),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.1, unit = "cm"))


########Count significant arms
Significant_counts <- pat_mat_diff_plot %>%
  filter(significance == "Significant", chr != "mTL") %>%
  count(analysis_level, name = "counts")

Significant_counts$analysis_level <- factor(Significant_counts$analysis_level, levels = c( "ASD","Sibling","Offspring"))

p2 <- ggplot(Significant_counts, aes(x = analysis_level, y = counts))+
  geom_col(width = 0.5,fill = "#EB6368")+
  theme_bw()+
  coord_flip()+
  scale_y_continuous(expand = c(0,0), limits = c(0,25), breaks = c(0,10,20))+
  labs(x = NULL, y = "Significant counts")+
  theme(axis.line = element_blank(),
        axis.ticks.x =  element_line(linewidth = 0.2),
        axis.ticks.length.x = unit(0.05,"cm"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 7),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.1, r = 0.5, unit = "cm"))


########### Correlation of the mean Pat-Mat difference in ChArmTLs between Sibling and ASD ###########
diff_sib_asd_corr_plot <- pat_mat_diff %>%
  select(chr,analysis_level,mean_diff) %>%
  filter(analysis_level != "Offspring") %>%
  pivot_wider(names_from = analysis_level, values_from = mean_diff) %>% 
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]"))) %>%
  filter( chr != "mTL") 
  
shapiro.test(diff_sib_asd_corr_plot$ASD)    #Normal distribution
shapiro.test(diff_sib_asd_corr_plot$Sibling)    # Non-normal distribution

cor.test(diff_sib_asd_corr_plot$ASD, diff_sib_asd_corr_plot$Sibling, alternative = "two.sided", method = "spearman", exact = TRUE, conf.level = 0.95)

p3 <- ggplot(diff_sib_asd_corr_plot, aes(x = Sibling, y = ASD, label = chr_number))+
  geom_point(size = 0.5, color = "#82ADD0")+
  geom_text_repel(size = 1, max.overlaps = Inf, box.padding = 0.05, point.padding = 0.05, segment.color = "grey40")+
  geom_smooth(method = 'lm', formula = y ~ x, se = T, linewidth = 0.35, color = "#82ADD0")+    #Add a fitted line
  geom_rug(outside = F, length = unit(0.02, "npc"), linewidth = 0.2, color = "#82ADD0")+    #Add a rug plot
  annotate(geom = "text", x = -300, y = 2300, label = "r = 0.52, P = 0.0011", size = 2, hjust = 0)+
  theme_classic()+
  scale_x_continuous(expand = c(0,0), limits = c(-500,3500), breaks = c(0,1000,2000,3000), labels = c("0","1,000","2,000","3,000"))+
  scale_y_continuous(expand = c(0,0), limits = c(-500,2500), breaks = c(0,1000,2000), labels = c("0","1,000","2,000"))+
  coord_fixed(ratio = 1)+
  labs(x = "Mean Difference of Sibling (Pat - Mat)", y = "Mean Difference of ASD (Pat - Mat)")+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))

ggsave(filename = "Supplementary Figure 8.pdf", plot = p3, width = 8, height = 8, units = "cm", dpi = 600)


############## Directionality of paternal vs maternal differences (by chromosome arm) ##############
########Define a directionality statistics function
diff_chr_direction <- function(data, label) {
  data %>%
    group_by(chr) %>%
    summarise(
      n_total = n(),
      n_Pat_longer = sum(Pat > Mat, na.rm = TRUE),    #Count Pat > Mat
      n_Mat_longer = sum(Pat < Mat, na.rm = TRUE),    #Count Pat < Mat
      n_Equal = sum(Pat == Mat, na.rm = TRUE),    #Count ties (continuous variables are rare in theory, but included for completeness
      tidy(binom.test(n_Pat_longer, n_total, p = 0.5, alternative = "greater")), .groups = 'drop') %>%
    mutate(analysis_level = label)
}

########Directionality statistics for paternal vs maternal differences
direction_chr_offspring <- diff_chr_direction(tlens_offspring_wide, "Offspring")    #Offspring
direction_chr_sib <- diff_chr_direction(subset(tlens_offspring_wide, clinical_group == "Sibling"), "Sibling")
direction_chr_asd <- diff_chr_direction(subset(tlens_offspring_wide, clinical_group == "ASD"), "ASD")

########Merge data
direction_chr <- bind_rows(direction_chr_offspring, direction_chr_sib, direction_chr_asd) %>%
  mutate(analysis_level = factor(analysis_level, levels = c("Offspring","Sibling", "ASD")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%   
  select(analysis_level,chr,n_total,n_Pat_longer,n_Mat_longer,n_Equal,estimate,conf.low,conf.high,statistic,p.value,method,alternative) %>%
  rename(proportion_Pat_longer = estimate) %>%
  arrange(analysis_level, chr)

write_xlsx(direction_chr, "direction_chr.xlsx") 


############### Visualization of paternal vs maternal difference directionality (by chromosome arm) ###############
direction_chr_plot <- direction_chr %>%
  select(analysis_level, chr, n_Pat_longer, n_Mat_longer) %>%
  pivot_longer(
    cols = c(n_Pat_longer, n_Mat_longer),
    names_to = "direction",
    values_to = "counts"
  ) %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         direction = factor(direction, levels = c("n_Pat_longer", "n_Mat_longer")))

p4 <- ggplot(direction_chr_plot, aes(x = chr_number, y = counts, color = direction))+
  geom_linerange(aes(ymin = 0, ymax = counts), linewidth = 0.35, position = position_dodge(width = 0.5))+    #Draw line segments
  geom_point(size = 1, position = position_dodge(width = 0.5))+    # Draw points
  facet_wrap2(analysis_level ~ ., nrow = 3, scales = "free_y", strip.position = "right",
              strip = strip_themed(background_y = elem_list_rect(fill = c("white", "#82ADD080","#EB636880"))))+
  scale_color_manual(values = c("#EB6368","#82ADD0"), labels = c("Pat > Mat","Pat < Mat"))+
  theme_bw()+    # Theme styling
  scale_x_discrete(limits = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                              "11p","11q","12p","12q","13q","14q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q","mTL"))+
  facetted_pos_scales(y = list(analysis_level == "Offspring" ~ scale_y_continuous(expand = c(0,0), limits = c(0,60), breaks = c(20,40,60)),
                               analysis_level == "Sibling" ~ scale_y_continuous(expand = c(0,0), limits = c(0,30), breaks = c(10,20,30)),
                               analysis_level == "ASD"     ~ scale_y_continuous(expand = c(0,0), limits = c(0,30), breaks = c(10,20,30))))+
  labs(x = NULL, y = "NUmber of samples", color = NULL)+
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
        legend.box.spacing = unit(0.01, "cm"),    #Distance between the legend box and the plot area
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.35),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))


p1_2 <- plot_grid(p1, p2, labels = c('',''), nrow = 1, rel_widths = c(1,0.18), align = 'h', axis = 'tb')
p4_NULL <- plot_grid(p4, NULL, labels = c('',''), nrow = 1, rel_widths = c(1,0.02), align = 'h', axis = 'tb')

p1_2_4 <- plot_grid(p1_2, p4_NULL, ncol = 1, labels =c('a','b'), label_size = 8, 
                    align = 'v', axis = 'lr', rel_heights = c(1,2))

ggsave(filename = "Figure 3.pdf", plot = p1_2_4, width = 17, height = 15, units = "cm", dpi = 600)


############# Directionality of paternal vs maternal differences (by individual) ############
direction_sample <- tlens_offspring_wide %>%
  filter(chr != "mTL") %>%
  group_by(sample) %>%
  summarise(
    n_total = n(),
    n_Pat_longer = sum(Pat > Mat, na.rm = TRUE),    #Count Pat > Mat
    n_Mat_longer = sum(Pat < Mat, na.rm = TRUE),    #Count Pat < Mat
    n_Equal = sum(Pat == Mat, na.rm = TRUE),    #Count ties (continuous variables are rare in theory, but included for completeness
    tidy(binom.test(n_Pat_longer, n_total, p = 0.5, alternative = "greater")), .groups = 'drop')

direction_sample <- direction_sample %>% 
  mutate(analysis_level = str_extract(sample, "[p/s]1")) %>% 
  mutate(analysis_level = factor(analysis_level, levels = c("s1","p1"), labels = c("Sibling", "ASD"))) %>%   
  select(analysis_level,sample,n_total,n_Pat_longer,n_Mat_longer,n_Equal,estimate,conf.low,conf.high,statistic,p.value,method,alternative) %>%
  rename(proportion_Pat_longer = estimate)

write_xlsx(direction_sample, "direction_sample.xlsx")


############# Visualization of paternal vs maternal difference directionality (by individual) ############# 
direction_sample_plot <- direction_sample %>%
  select(analysis_level, sample, n_Pat_longer, n_Mat_longer) %>%
  pivot_longer(
    cols = c(n_Pat_longer, n_Mat_longer),
    names_to = "direction",
    values_to = "counts"
  ) %>%
  mutate(direction = factor(direction, levels = c("n_Mat_longer","n_Pat_longer")))

p5 <- ggplot(data = direction_sample_plot, aes(x = sample, y = counts, fill = direction))+
  geom_bar(stat = "identity", position = "fill")+
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.2, color = "black")+
  facet_wrap2(analysis_level ~ ., nrow = 2, scales = "free_x", strip.position = "right",
              strip = strip_themed(background_y = elem_list_rect(fill = c("#82ADD080","#EB636880"))))+
  scale_fill_manual(values = c("#82ADD0","#EB6368"), labels = c("Pat < Mat","Pat > Mat"))+
  theme_bw()+    # Theme styling
  scale_y_continuous(expand = c(0,0.03))+
  labs(x = NULL, y = "ChArmTL proportion", fill = NULL)+
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
        legend.box.spacing = unit(0.01, "cm"),    #Distance between the legend box and the plot area
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.35),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))

ggsave(filename = "Supplementary Figure 9.pdf", plot = p5, width = 15, height = 8, units = "cm", dpi = 600)

