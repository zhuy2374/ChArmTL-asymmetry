########################## set the work directory ##########################
rm(list=ls())
setwd("Path/to/your_workdir")
getwd()


############################ Load R package and data ###########################
required_packages <- c("readxl", "dplyr", "ggpubr", "writexl", "tidyverse", "rstatix",
                       "car", "ggplot2", "ggsignif", "broom", "cowplot")

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
library(ggpubr)
library(writexl)
library(rstatix)    #shapiro_test()
library(car)    #Levene test
library(ggplot2)
library(ggsignif)
library(tidyverse)    #pivot_wider()
library(broom)
library(cowplot)


telogator2_0 <- read_delim("merged_tlens_by_allele_CEPH-1463.tsv", delim = "\t", escape_double = FALSE, 
                           col_types = cols(position = col_character(), 
                                            allele_id = col_character(), 
                                            read_lengths = col_character()),
                           trim_ws = TRUE)    # 1605 rows

family_info <- read_excel("family_info_CEPH-1463.xlsx", range = "A1:G24")    #columns:generation. sample, sex, age_at_blood_draw, father, mother, assembly_algorithm

############################### QC and filtering ###############################
######## excluding sample without parental information, 1605 rows ---> 1164 rows
telogator2_0 <- subset(telogator2_0, sample != "200080" & sample != "200100" & sample != "NA12889" & sample != "NA12890" & sample != "NA12891" & sample != "NA12892")    

######## excluding inaccurate estimations
#### unmapped
telogator2_1 <- telogator2_0[!grepl("chrUq|chrUp", telogator2_0$`#chr`), ]    # 1164 rows ---> 1164 rows

#### ambiguous alignments
telogator2_1 <- telogator2_1[!grepl(",|，", telogator2_1$`#chr`), ]    # 1164 rows ---> 1161 rows

#### putative interstitial telomeric sequences
telogator2_1 <- telogator2_1[!grepl("i$", telogator2_1$allele_id), ]    # 1161 rows ---> 1150 rows

#### tvr_len was zero
telogator2_1 <- telogator2_1[telogator2_1$tvr_len != 0, ]    # 1150 rows ---> 1131 rows


######## excluding sex chromosomes
telogator2_2 <- telogator2_1[!grepl("X|Y", telogator2_1$`#chr`), ]    # 1131 rows ---> 1087 rows


######## excluding outliers
summary(telogator2_2$TL_p75)    #Descriptive analysis of telomere length

hist(telogator2_2$TL_p75, prob = T)    #Histogram, non-normal distribution
lines(density(na.omit(telogator2_2$TL_p75)), col = 'blue')   #Density curve


#### Outlier exclusion (By generation)
telogator2_3 <- telogator2_2 %>%
  mutate(QC_status = if_else(TL_p75 <= 100, "outlier", "passed"))   # Initial absolute threshold filtering: values <= 100 are directly marked as outliers

telogator2_3 <- telogator2_3 %>%
  left_join(family_info %>% select(sample, generation), by = "sample") %>%　　# Merge 'generation' info from family_info
  group_by(generation) %>%
  mutate(    # calculate the Q1 and Q3 values for each group dynamically (based solely on the current "passed" data)
    Q1 = quantile(TL_p75[QC_status == "passed"], 0.25, na.rm = TRUE),
    Q3 = quantile(TL_p75[QC_status == "passed"], 0.75, na.rm = TRUE),
    IQR_val = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR_val,    # dynamic calculation of upper and lower bounds
    upper_bound = Q3 + 1.5 * IQR_val,
    QC_status = case_when(    # keep the original outliers, and add new outliers that exceed the dynamic boundaries
      QC_status == "outlier" ~ "outlier",
      TL_p75 < lower_bound | TL_p75 > upper_bound ~ "outlier",
      TRUE ~ "passed")) %>%
  ungroup() %>%
  select(-Q1, -Q3, -IQR_val, -lower_bound, -upper_bound)

table(telogator2_3$QC_status)    # 28 outliers + 1059 passed


########Distribution and normality check after removing outliers
telogator2_4 <- subset(telogator2_3, QC_status == "passed")     # 1087 rows ---> 1059 rows

summary(telogator2_4$TL_p75)

hist(telogator2_4$TL_p75, prob = TRUE)
lines(density(na.omit(telogator2_4$TL_p75)), col = "blue")

write_xlsx(telogator2_4, "tlens_outlier_filtered_CEPH-1463.xlsx")


############################### Stage 2 QC  ###############################
tlens_phased_filtered <- read_xlsx("tlens_phased_filtered.xlsx")    # The result file after excluding inaccurately phased ChArmTLs.  1032 rows

tlens_chr_filtered <- subset(tlens_unphased_filtered, chr != "chr13p" & chr != "chr14p" & chr != "chr15p" & chr != "chr21p" & chr != "chr22p" & chr != "chr15q")    # 1032 rows ---> 923 rows


############ Number of ChArmTLs for different haplotypes per sample ############
alleles_count <- as.data.frame(table(tlens$sample, tlens$haplotype))    # Count the number of ChArmTLs for different haplotypes in each offspring
colnames(alleles_count) <- c("sample", "haplotype","counts")    # Modify default column names

#### Mean and Median
mean(alleles_count$counts)    # 27.14706
median(alleles_count$counts)    # 27

alleles_count %>%
  group_by(haplotype) %>%
  summarise(
    n = n(),
    mean = mean(counts),
    median = median(counts)
  )  

#### Test whether there are differences between haplotypes
alleles_count %>%
  group_by(haplotype) %>%
  shapiro_test(counts)    # Normality test
leveneTest(counts ~ haplotype, data = alleles_count)    # Homogeneity of variance test
t.test(counts ~ haplotype, data = alleles_count)    # Compare differences between groups

# p1 <- ggplot(data = alleles_count, aes(x = haplotype, y = counts, fill = haplotype))+
#   geom_boxplot(linewidth = 0.3, median.linewidth = 0.3, width = 0.5, outlier.alpha = 0, show.legend = FALSE)+
#   scale_fill_manual(values = c("#82ADD0","#EB6368"))+
#   theme_classic()+
#   scale_y_continuous(expand = c(0,0) ,limits = c(15,40), breaks = c(20,25,30,35))+
#   labs(x = NULL, y = 'Number of ChArmTLs')+
#   theme(axis.line = element_line(linewidth = 0.3),
#         axis.ticks =  element_line(linewidth = 0.2),
#         axis.ticks.length = unit(0.05,"cm"),
#         axis.text.x = element_text(size = 6),
#         axis.text.y = element_text(size = 6),
#         axis.title.y = element_text(size = 7),
#         plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))+
#   geom_signif(y_position = 36, xmin = 1, xmax = 2, annotations = c("ns"),
#               tip_length = 0.01, size = 0.2, textsize = 1.5)
# 
# # ggplot_build(p1)$data[[1]]    # Check the end point of the upper whisker (ymax) of the boxplot to determine the y_position for the geom_signif function


############################### Parental origin analysis ###############################
######## Calculate offspring haplotype mTL and extract ChArmTLs
#### Extract ChArmTLs
tlens_hap_arm <- tlens %>%
  select(sample, haplotype, chr, TL_p75)

#### Offspring paternal and maternal mTL
tlens_hap_mTL <- tlens %>%
  group_by(sample, haplotype) %>%
  summarise(TL_p75 = mean(TL_p75, na.rm = TRUE), .groups = "drop") %>%
  mutate(chr = "mTL") %>%
  select(sample, haplotype, chr, TL_p75)

######## Merge data
tlens_hap <- bind_rows(tlens_hap_arm, tlens_hap_mTL) %>%
  mutate(haplotype = factor(haplotype, levels = c("Pat","Mat")),
         chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                      "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                      "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                      "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%
  arrange(sample, haplotype, chr)

write_xlsx(tlens_hap, "tlens_chr_filtered_hap_CEPH-1463.xlsx")    # Save data


######## Data format processing
tlens_hap_wide <- tlens_hap %>%
  pivot_wider(names_from = haplotype, values_from = TL_p75) %>%
  drop_na(Pat, Mat)   # Keep only rows where both paternal and maternal data exist


######## Custom function to compare paternal and maternal ChArmTLs
pat_mat_difference <- function(data) {
  data %>%
    group_by(chr) %>%
    summarise(n_pairs = n(),  # 统计该分组下的配对样本数
              mean_Pat = mean(Pat), sd_Pat = sd(Pat),    #统计Pat组
              mean_Mat = mean(Mat), sd_Mat = sd(Mat),    #统计Mat组
              tidy(t.test(Pat, Mat, paired = TRUE)), .groups = 'drop')
}


######## Paternal and maternal difference comparison
diff_offspring <- pat_mat_difference(tlens_hap_wide)

#### Merge data and add adjusted P-values
pat_mat_diff <- diff_offspring %>%
  mutate(p_adj_fdr = case_when(
    chr == "mTL" ~ p.value,
    TRUE ~ p.adjust(p.value, method = "fdr")),    # P-value adjustment
    chr = factor(chr, levels = c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q",
                                 "chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q",
                                 "chr11p","chr11q","chr12p","chr12q","chr13q","chr14q","chr16p","chr16q","chr17p",
                                 "chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21q","chr22q","mTL"))) %>%   
  select(chr,n_pairs,mean_Pat,sd_Pat,mean_Mat,sd_Mat,estimate,conf.low,conf.high,statistic,p.value,p_adj_fdr,method,alternative) %>%
  rename(mean_diff = estimate) %>%
  arrange(chr)

write_xlsx(pat_mat_diff, "pat_mat_diff_CEPH-1463.xlsx")


############## Visualization of paternal and maternal differences ##############
diff_offspring_plot <- pat_mat_diff %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         mean_diff_sign = ifelse(mean_diff > 0, "Pat > Mat", "Pat < Mat"),
         significance = ifelse(p.value < 0.05, "Significant", "Not Significant"),
         logP = ifelse(p.value < 0.05, -log10(p.value), 1))      # Effect direction

diff_offspring_plot$mean_diff_sign <- factor(diff_offspring_plot$mean_diff_sign, levels = c("Pat > Mat", "Pat < Mat"))
diff_offspring_plot$significance <- factor(diff_offspring_plot$significance, levels = c("Significant", "Not Significant"))

p2 <- ggplot(diff_offspring_plot, aes(x = chr_number, y = abs(mean_diff), color = mean_diff_sign, shape = significance, size = logP))+
  geom_point(stroke = 0.35)+
  scale_color_manual(values = c("#EB6368","#82ADD0"))+
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 1))+
  theme_bw()+
  scale_x_discrete(limits = c("1q","3q","4p","4q","5p","5q","6p","6q","7p","8p","9p","9q","10q",
                              "11p","12p","16p","16q","17q","18p","19p","19q","20p","21q","22q", "mTL"))+    # the 24 paternally biased arms identified in our primary cohort + mTL
  scale_y_continuous(limits = c(0,4000), breaks = c(0,1000,2000,3000,4000), labels = c("0","1,000","2,000","3,000","4,000"))+
  scale_size_continuous(range = c(0.5,3), breaks = c(1, -log10(0.05), 2), labels = c("ns","0.05","0.01"))+    # Control minimum and maximum bubble size
  labs(x = NULL, y = "Mean difference (bp)", size = "P", color = NULL)+
  guides(color = guide_legend(override.aes = list(size = 1.5), position = "top", order = 1), 
         size = guide_legend(override.aes = list(shape = c(1,16,16), stroke = 0.35), position = "top", order = 2),
         shape = "none")+    # Modify scatter size in legend and legend position
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
        legend.box.spacing = unit(0.01, "cm"),    # Distance between legend box and plot area
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))


########Arms with significant differences (5q、7p、21q)
tlens_chr_sign <- tlens_hap %>% 
  filter(chr %in% c("chr5q","chr7p","chr21q")) %>%          #  Select the target arms
  group_by(sample, chr) %>%                # Group by sample and chromosome
  filter(n_distinct(haplotype) == 2) %>%       # Retain the group that contains both Pat and Mat
  ungroup()

p3 <- ggplot(data = tlens_chr_sign, aes(x = haplotype, y = TL_p75, fill = haplotype, color = haplotype))+
  geom_boxplot(linewidth = 0.3, median.linewidth = 0.3, width = 0.5, outlier.alpha = 0, show.legend = FALSE)+
  geom_point(size = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), show.legend = FALSE)+
  facet_wrap(~ chr, ncol = 3, scales = "free_y")+   # 3列，纵轴可自由缩放
  scale_fill_manual(values = c("#82ADD080","#EB636880"))+
  scale_colour_manual(values = c("#82ADD0","#EB6368"))+
  theme_bw()+
  facetted_pos_scales(y = list(chr == "chr5q" ~ scale_y_continuous(expand = c(0,0), limits = c(0,10000), breaks = c(2500,5000,7500,10000), labels =c("2,500","5,000","7,500","10,000")),
                               chr == "chr7p" ~ scale_y_continuous(expand = c(0,0), limits = c(0,12500), breaks = c(2500,5000,7500,10000,12500), labels =c("2,500","5,000","7,500","10,000","12,500")),
                               chr == "chr21q" ~ scale_y_continuous(expand = c(0,0), limits = c(0,10000), breaks = c(2500,5000,7500,10000), labels =c("2,500","5,000","7,500","10,000"))))+
  labs(x = NULL, y = 'ChArmTL (bp)')+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6,),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.35),
        panel.spacing = unit(1, "cm"),    #分面之前的间距
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm")) +
  geom_signif(data = subset(tlens_chr_sign, chr == "chr5q"),  # 指定分面
              comparisons = list(c("Pat", "Mat")), annotations = "*",y_position = (8715 + 200), 
              tip_length = 0, size = 0.3, color = "black", textsize = 2, show.legend = F)+
  geom_signif(data = subset(tlens_chr_sign, chr == "chr7p"),
              comparisons = list(c("Pat", "Mat")), annotations = "*",y_position = (10404 + 200), 
              tip_length = 0, size = 0.3, color = "black", textsize = 2, show.legend = F)+
  geom_signif(data = subset(tlens_chr_sign, chr == "chr21q"),
              comparisons = list(c("Pat", "Mat")), annotations = "**",y_position = (7922 + 200),
              tip_length = 0, size = 0.3, color = "black", textsize = 2, show.legend = F)

# ggplot_build(p3)$data[[1]]    #查看箱体上须终点（ymax），确定geom_signif函数的y_position


p2_3 <- plot_grid(p2, p3, ncol = 1, labels = c('a', 'b'), label_size = 8, rel_heights = c(0.8,1),align = 'h', axis = 'tb')
ggsave(filename = "Figure 4.pdf", plot = p2_3, width = 12, height = 10, units = "cm", dpi = 600)

