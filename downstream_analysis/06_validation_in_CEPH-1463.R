############################ set the work directory ############################
rm(list=ls())
setwd("Path/to/your_workdir")
getwd()


############################ Load R package and data ###########################
required_packages <- c("readxl", "dplyr", "writexl", "rstatix","car", "ggplot2", 
                       "ggsignif", "RColorBrewer", "lmer", "tidyverse", "cowplot")

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
library(rstatix)
library(car)
library(ggplot2)
library(ggsignif)
library(RColorBrewer)
library(lme4)
library(tidyverse)
library(cowplot)

telogator2_0 <- read_delim("merged_tlens_by_allele_G4.tsv", delim = "\t", escape_double = FALSE, 
                           col_types = cols(position = col_character(), 
                                            allele_id = col_character(), 
                                            read_lengths = col_character()),
                           trim_ws = TRUE)    # 1605 rows

telogator2_0 <- telogator2_0 %>% 
  select(Sample, `#chr`, allele_id, TL_p75, tvr_len) %>%
  rename(chr = `#chr`)

colnames(telogator2_0)
# "generation"， "sample", "haplotype"， "chr", "allele_id","TL_p75", "tvr_len"


############################### QC and filtering ###############################
######## excluding sample without parental information
telogator2_1 <- telogator2_0 %>%
  filter(generation != "G1",
         !sample %in% c("200080","200100"))    # 1605 rows ---> 1164 rows


######## excluding inaccurate estimations
telogator2_2 <- telogator2_1 %>%
  filter(!grepl("chrUq|chrUp", chr),     # unmapped, 1164 rows ---> 1164 rows
         !grepl(",|，", chr),            # ambiguous alignments, 1164 rows ---> 1161 rows
         !grepl("i$", allele_id),        # putative interstitial telomeric sequences, 1161 rows ---> 1150 rows
         tvr_len != 0)                   # tvr_len was zero, 1150 rows ---> 1131 rows


######## excluding sex chromosomes
telogator2_3 <- telogator2_2[!grepl("X|Y", telogator2_2$chr), ]    # 1131 rows ---> 1087 rows


######## excluding outliers
summary(telogator2_3$TL_p75)

hist(telogator2_3$TL_p75, prob = T)    #Histogram, non-normal distribution
lines(density(na.omit(telogator2_3$TL_p75)), col = 'blue')


#### Outlier exclusion (By Chromosome arm)
telogator2_4 <- telogator2_3 %>%
  mutate(QC_status = if_else(TL_p75 <= 100, "outlier","passed")) %>%    # Initial absolute threshold filtering: values <= 100 are directly marked as outliers
  group_by(chr) %>%     # group the data by chr for calculation
  mutate(Q1 = quantile(TL_p75[QC_status == "passed"], 0.25, na.rm = TRUE),    # calculate the Q1 and Q3 values for each chr dynamically (based solely on the current "passed" data)
         Q3 = quantile(TL_p75[QC_status == "passed"], 0.75, na.rm = TRUE),
         IQR_val = Q3 - Q1,
         lower_bound = Q1 - 1.5 * IQR_val,    # dynamic calculation of upper and lower bounds
         upper_bound = Q3 + 1.5 * IQR_val,
         QC_status = case_when(QC_status == "outlier" ~ "outlier",    # keep the original outliers, and add new outliers that exceed the dynamic boundaries
                               TL_p75 < lower_bound | TL_p75 > upper_bound ~ "outlier",
                               TRUE ~ "passed")) %>%
  ungroup() %>%
  select(-Q1, -Q3, -IQR_val, -lower_bound, -upper_bound)

table(telogator2_4$QC_status, telogator2_4$generation)    # 34 outliers [8(G2) + 8(G3) + 18(G4)] + 1053 passed


########Distribution and normality check after removing outliers
telogator2_4_passed <- subset(telogator2_4, QC_status == "passed")     # 1087 rows ---> 1053 rows

summary(telogator2_4_passed$TL_p75)

hist(telogator2_4_passed$TL_p75, prob = T)
lines(density(na.omit(telogator2_4_passed$TL_p75)), col = 'blue')

write_xlsx(telogator2_4_passed, "tlens_outlier_filtered_G4.xlsx")


########################### Manual haplotype phasing ###########################


################################## Stage 2 QC ##################################
telogator2_5 <- read_excel("tlens_inaccurate_phased_filtered_G4.xlsx",
                           col_types = c("text", "text", "text", "text", "text","numeric", "numeric"))    # The result file after excluding inaccurately phased ChArmTLs. 1025 rows  

colnames(telogator2_5)    
# "generation", "sample", "haplotype", "chr", "allele_id", "TL_p75", "tvr_len"

tlens_chr_filtered_38 <- telogator2_5 %>%
  filter(!chr %in% c("chr13p","chr14p","chr15p","chr15q","chr21p","chr22p"))    # 1025 rows ---> 913 rows


############ Number of ChArmTLs for different haplotypes per sample ############
ChArmTLs_count_by_haplotype <- tlens_chr_filtered_38 %>%
  count(sample, haplotype, name = "n_ChArmTLs", .drop = TRUE) %>%
  mutate(haplotype = factor(haplotype, levels = c("Pat","Mat")))

write_xlsx(ChArmTLs_count_by_haplotype, "ChArmTLs_count_by_haplotype_G4.xlsx")

#### Mean and Median
ChArmTLs_count_by_haplotype %>%
  summarise(n = n(),
            mean = sprintf("%.2f", mean(n_ChArmTLs, na.rm = TRUE)),
            median = median(n_ChArmTLs))

ChArmTLs_count_by_haplotype %>%
  group_by(haplotype) %>%
  summarise(n = n(),
            mean = sprintf("%.2f", mean(n_ChArmTLs, na.rm = TRUE)),
            median = median(n_ChArmTLs))

#### Test whether there are differences between haplotypes
ChArmTLs_count_by_haplotype %>%
  group_by(haplotype) %>%
  shapiro_test(n_ChArmTLs)
leveneTest(n_ChArmTLs ~ haplotype, data = ChArmTLs_count_by_haplotype)
t.test(n_ChArmTLs ~ haplotype, data = ChArmTLs_count_by_haplotype)

p1 <- ggplot(data = ChArmTLs_count_by_haplotype, aes(x = haplotype, y = n_ChArmTLs, fill = haplotype))+
  geom_boxplot(linewidth = 0.3, median.linewidth = 0.3, width = 0.5, outlier.alpha = 0, show.legend = FALSE)+
  scale_fill_manual(values = c("#82ADD0","#EB6368"))+
  scale_y_continuous(expand = c(0,0) ,limits = c(15,40), breaks = c(20,25,30,35))+
  labs(x = NULL, y = 'Number of ChArmTLs')+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.5, b = 0.25, l = 0.5, r = 0.25, unit = "cm"))+
  geom_signif(y_position = 34, xmin = 1, xmax = 2, annotations = c("ns"),
              tip_length = 0.01, size = 0.2, textsize = 1.5)

# ggplot_build(p1)$data[[1]]    #Check the end point of the upper whisker (ymax) of the boxplot to determine the y_position for the geom_signif function


######### Haplotype coverage of ChArmTLs (Only 18 arms were used for candidate validation) ##########
tlens_chr_filtered_18 <- tlens_chr_filtered_38 %>%
  filter(chr %in% c("chr1q","chr4p","chr5p","chr5q","chr6p","chr7p","chr9p","chr10p","chr10q","chr12p",
                    "chr16p","chr16q","chr17q","chr18p","chr19p","chr20p","chr21q","chr22q"))    # 913 rows ---> 409 rows

write_xlsx(tlens_chr_filtered_18, "tlens_chr_filtered_18_G4.xlsx")

arm_coverage_18 <- tlens_chr_filtered_18 %>%
  count(chr, generation, name = "counts") %>%
  mutate(coverage = counts / 34,
         chr_number = substr(chr, 4, nchar(as.character(chr))),
         generation = factor(generation, levels =c("G2","G3","G4")))

write_xlsx(arm_coverage_18, "arm_coverage_18_G4.xlsx")

p2 <- ggplot(data = arm_coverage_18, aes(x = reorder(chr_number, -coverage), y = coverage))+
  geom_bar(stat = "identity", fill = "#82ADD0")+
  scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0.25,0.5,0.75,1), labels = c(25,50,75,100),
                     sec.axis = sec_axis(transform = ~ .*34, name = "Number of alleles", breaks = c(10,20,30)))+
  labs(x = NULL, y = "Coverage (%)", fill = NULL)+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        axis.title.y.right = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),    # 图例框离 plot 区域的距离
        legend.position = "top",
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.25, l = 0.25, r = 0.5, unit = "cm"))


############################## Length Distribution #############################
######## plot by chr
tlens_chr_filtered_plot <- tlens_chr_filtered_18 %>%
  mutate(chr_number = substr(chr,4,nchar(as.character(chr))),
         chr_number = factor(chr_number, levels = c("1q","4p","5p","5q","6p","7p","9p","10p","10q","12p",
                                                    "16p","16q","17q","18p","19p","20p","21q","22q")))

p3 <- ggplot(data = tlens_chr_filtered_plot, aes(x = chr_number, y = TL_p75))+
  geom_point(aes(color = chr_number), position = position_jitter(0.3), size = 0.1, alpha = 0.5, show.legend = F)+
  geom_boxplot(color = "black", fill = NA, linewidth = 0.25, median.linewidth = 0.25, width = 0.5, outlier.alpha = 0, show.legend = FALSE)+
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Paired"))(length(unique(tlens_chr_filtered_plot$chr_number)))) +
  scale_y_continuous(expand = c(0,0) ,limits = c(0,12000), breaks = c(0,4000,8000,12000), labels =c("0","4,000","8,000","12,000"))+
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
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.5, r = 0.5, unit = "cm"))

p1_2 <- plot_grid(p1, p2, labels = c('a','b'), label_size = 8, nrow = 1, rel_widths = c(0.25,1), align = 'h', axis = 'tb')
p3_NULL <- plot_grid(p3, NULL, labels = c('c',''), label_size = 8, nrow = 1, rel_widths = c(1,0.025), align = 'h', axis = 'tb')

p1_2_3 <- plot_grid(p1_2, p3_NULL, ncol = 1, rel_heights = c(1,1), align = 'v', axis = 'lr')
ggsave(filename = "Supplementary Figure 9.pdf", plot = p1_2_3, width = 15, height = 8, units = "cm", dpi = 600)


######## Statistics of ChArmTL
tlens_by_chr_stats <- tlens_chr_filtered_plot %>%
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
  bind_rows(tlens_chr_filtered_18 %>%
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


######## Inter-arm heterogeneity of ChArmTLs
model_0 <- lmer(TL_p75 ~ 1 + (1 | sample), data = tlens_chr_filtered_18, REML = FALSE)
model_1 <- lmer(TL_p75 ~ chr + (1 | sample), data = tlens_chr_filtered_18, REML = FALSE)

anova(model_0,model_1)    # The likelihood ratio test (LRT) shows that there is a significant inter-arm heterogeneity in ChArmTLs


############################### Parental origin analysis ###############################
######## Calculate the haplotype-specific mTL of offspring and extract the ChArmTLs
tlens_chr_filtered_hap <- tlens_chr_filtered_18 %>% 
  select(sample, haplotype, chr, TL_p75) %>%
  bind_rows(tlens_chr_filtered_38 %>%
              group_by(sample, haplotype) %>%
              summarise(TL_p75 = mean(TL_p75, na.rm = TRUE), .groups = "drop") %>%
              mutate(chr = "mTL") %>%
              select(sample, haplotype, chr, TL_p75)) %>%
  mutate(haplotype = factor(haplotype, levels = c("Pat","Mat")),
         chr = factor(chr, levels = c("chr1q","chr4p","chr5p","chr5q","chr6p","chr7p","chr9p","chr10p","chr10q","chr12p",
                                      "chr16p","chr16q","chr17q","chr18p","chr19p","chr20p","chr21q","chr22q","mTL"))) %>%
  arrange(sample, haplotype, chr)

write_xlsx(tlens_chr_filtered_hap, "tlens_chr_filtered_hap_G4.xlsx")


######################### The statistic of Δ(Pat - Mat) ########################
tlens_Pat_Mat_diff <- tlens_chr_filtered_hap %>%
  pivot_wider(names_from = haplotype, values_from = TL_p75) %>%
  drop_na(Pat, Mat) %>%
  mutate(diff = Pat-Mat,
         chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")))


######### Analysis of differences between parental and maternal origins ########
######## A custom function compares the telomere lengths of the paternal and maternal origins
pat_mat_difference <- function(data) {
  data %>%
    group_by(chr) %>%
    summarise(n_pairs = n(),
              mean_Pat = mean(Pat), sd_Pat = sd(Pat),
              mean_Mat = mean(Mat), sd_Mat = sd(Mat),
              tidy(t.test(Pat, Mat, paired = TRUE, alternative = "greater")), .groups = 'drop')
}


######## Paternal and maternal difference comparison
diff_pat_mat_G4 <- pat_mat_difference(tlens_Pat_Mat_diff)

#### Merge data and add adjusted P-values
diff_pat_mat_G4 <- diff_pat_mat_G4 %>%
  mutate(p_adj_fdr = case_when(chr == "mTL" ~ p.value,
                               TRUE ~ p.adjust(p.value, method = "fdr")),    # P 值矫正
         chr = factor(chr, levels = c("chr1q","chr4p","chr5p","chr5q","chr6p","chr7p","chr9p","chr10p","chr10q","chr12p",
                                      "chr16p","chr16q","chr17q","chr18p","chr19p","chr20p","chr21q","chr22q","mTL"))) %>%
  select(chr,n_pairs,mean_Pat,sd_Pat,mean_Mat,sd_Mat,estimate,conf.low,conf.high,statistic,p.value,p_adj_fdr,method,alternative) %>%
  rename(mean_diff = estimate) %>%
  arrange(chr)

write_xlsx(diff_pat_mat_G4, "diff_pat_mat_G4.xlsx")    # 保存数据


############## Visualization of paternal and maternal differences ##############
diff_pat_mat_G4_plot <- diff_pat_mat_G4 %>%
  mutate(chr_number = if_else(chr == "mTL", "mTL", str_extract(chr, "\\d+[pq]")),
         mean_diff_sign = ifelse(mean_diff > 0, "Pat > Mat", "Pat < Mat"),
         mean_diff_sign = factor(mean_diff_sign, levels =c("Pat > Mat", "Pat < Mat")),
         significance = case_when(p.value < 0.05 ~ "Nominal P < 0.05",
                                  TRUE ~ "Not significant"),
         significance = factor(significance, levels = c("Nominal P < 0.05", "Not significant")),
         logP = ifelse(p.value < 0.05, -log10(p.value), 1))

p4 <- ggplot(diff_pat_mat_plot, aes(x = chr_number, y = abs(mean_diff), color = mean_diff_sign, shape = significance, size = logP))+
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
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))


#### The direction of paternal and maternal differences (by chromosome arms) ####
######## Define the directional statistical function
diff_direction_chr <- function(data) {
  data %>%
    group_by(chr) %>%
    summarise(
      n_total = n(),
      n_Pat_longer = sum(Pat > Mat, na.rm = TRUE),
      n_Mat_longer = sum(Pat < Mat, na.rm = TRUE),
      n_Equal = sum(Pat == Mat, na.rm = TRUE),
      tidy(binom.test(n_Pat_longer, n_total, p = 0.5, alternative = "greater")), .groups = 'drop')
}


######## Directional statistics of paternal and maternal differences
direction_chr_G4 <- diff_direction_chr(tlens_Pat_Mat_diff)

direction_chr_G4 <- direction_chr_G4 %>%
  mutate(chr = factor(chr, levels = c("chr1q","chr4p","chr5p","chr5q","chr6p","chr7p","chr9p","chr10p","chr10q","chr12p",
                                      "chr16p","chr16q","chr17q","chr18p","chr19p","chr20p","chr21q","chr22q","mTL"))) %>%   
  select(chr,n_total,n_Pat_longer,n_Mat_longer,n_Equal,estimate,conf.low,conf.high,statistic,p.value,method,alternative) %>%
  rename(proportion_Pat_longer = estimate) %>%
  arrange(chr)

write_xlsx(direction_chr_G4, "direction_chr_G4.xlsx")

