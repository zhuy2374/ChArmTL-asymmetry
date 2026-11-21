########################## set the work directory ##########################
rm(list=ls())
setwd("Path/to/your_workdir")
getwd()


############################ Load R package and data ###########################
install.packages(c("readxl", "dplyr", "ggplot2", "ggh4x", 
                   "writexl", "cowplot", "rstatix", "car", 
                   "ggsignif"))

library(readxl)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(writexl)
library(cowplot)
library(rstatix)
library(car)
library(ggsignif)


telogator2_4 <- read_xlsx("tlens_outlier_filtered.xlsx")    # 9296 rows

colnames(telogator2_4_passed)    
# "sample", "chr", "allele_id", "TL_p75", "read_TLs", "tvr_len", "number_of_supporting_reads", "family_role", "QC_status"


############### Statistics and visualization of supporting reads ################
########Mean values and quantiles
telogator2_4 %>%
  summarise(n = n(),
            mean = sprintf("%.2f", mean(number_of_supporting_reads, na.rm = TRUE)),
            Q1 = quantile(number_of_supporting_reads, .25),
            median = median(number_of_supporting_reads),
            Q3 = quantile(number_of_supporting_reads, .75))

telogator2_4 %>%
  group_by(family_role) %>%
  summarise(
    n = n(),
    mean = sprintf("%.2f", mean(number_of_supporting_reads)),
    Q1 = quantile(number_of_supporting_reads, .25),
    median = median(number_of_supporting_reads),
    Q3 = quantile(number_of_supporting_reads, .75)
  )

number_of_supporting_reads_plot <- telogator2_4 %>%
  mutate(family_role = factor(family_role, levels = c("Father","Mother","Sibling","ASD")),
         chr_number = substr(chr,4,nchar(as.character(chr))),
         chr_number = factor(chr_number, levels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q",
                                                    "6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
                                                    "11p","11q","12p","12q","13p","13q","14p","14q","15p","15q",
                                                    "16p","16q","17p","17q","18p","18q","19p","19q","20p","20q",
                                                    "21p","21q","22p","22q")))

p1 <- ggplot(number_of_supporting_reads_plot, aes(x = chr_number, y = number_of_supporting_reads, color = family_role))+
  geom_boxplot(linewidth = 0.3, median.linewidth = 0.3, outlier.alpha = 0, show.legend = FALSE)+
  facet_wrap2(family_role ~ ., nrow = 4, strip.position = "right",
              strip = strip_themed(background_y = elem_list_rect(fill = c("#FCD28c80","#A1B5A380","#82ADD080","#EB636880"))))+
  scale_color_manual(values = c("#FCD28c","#A1B5A3","#82ADD0","#EB6368"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,20), breaks = c(5,10,15,20))+
  labs(x = NULL, y = "Number of supporting reads")+
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

ggsave(filename = "Supplementary Figure 1.pdf", plot = p1, width = 17, height = 15, units = "cm", dpi = 600)


###################### Allele count statistics for chr22p ######################
allele_count_22p <- telogator2_4 %>%
  filter(chr == "chr22p") %>%
  count(sample, name = "n_allele") %>%         # 统计每个样本的 22p 的 allele 个数
  mutate(n_allele = as.character(n_allele))    # 转换成字符

write_xlsx(allele_count_22p, "allele_count_22p.xlsx")

p2 <- ggplot(allele_count_22p, aes(x = n_allele))+
  geom_bar(width = 0.7, fill = "#82ADD0")+
  scale_y_continuous(expand = c(0,0), limits = c(0,35), breaks = c(10,20,30))+
  labs(x = "Number of ChArmTL alleles on 22p", y = "Number of samples")+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.5, b = 0.25, l = 0.5, r = 0.5, unit = "cm"))


###### Determine the parental origins and transmission status of ChArmTLs ######


######################## Haplotype coverage of ChArmTLs ########################
telogator2_5 <- read_excel("tlens_unphased_filtered.xlsx", 
                           col_types = c("text","text","text","numeric","numeric","numeric","text","text","text","text"))    # 8206 rows

colnames(telogator2_5)    
# "sample", "chr", "allele_id", "TL_p75", "tvr_len",
# "number_of_supporting_reads", "family_role", "parent_of_origin", "trans_status_ASD", "trans_status_Sib"

arm_coverage <- telogator2_5 %>%
  count(chr, family_role, name = "counts") %>%
  mutate(coverage = counts / 248,
         chr_number = substr(chr, 4, nchar(as.character(chr))),
         family_role = factor(family_role, levels =c("Father","Mother","Sibling","ASD")))

write_xlsx(arm_coverage, "arm_coverage.xlsx")

p3 <- ggplot(data = subset(arm_coverage, chr_number != "22p"), aes(x = reorder(chr_number, -coverage), y = coverage, fill = family_role))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#FCD28c","#A1B5A3","#82ADD0","#EB6368"), name = NULL)+
  scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0.25,0.5,0.75,1), labels = c(25,50,75,100),
                     sec.axis = sec_axis(transform = ~ .*248, name = "Number of ChArmTL alleles", breaks = c(50,100,150,200)))+
  labs(x = NULL, y = "Coverage (%)")+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6, hjust = 1, angle = 70),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        axis.title.y.right = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = margin(t = 0.01, b = 0.01, l = 0.01, r = 0.01, unit = "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        legend.position = "top",
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.35),
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.5, r = 0.5, unit = "cm"))

p2_null <- plot_grid(p2, NULL, nrow = 1, labels = c('a',''), label_size = 8, rel_widths = c(1,1.75))
p2_3 <- plot_grid(p2_null, p3, ncol = 1, labels = c('','b'), label_size = 8, rel_heights = c(1,1.25))

ggsave(filename = "Supplementary Figure 2.pdf", plot = p2_3, width = 17, height = 9, units = "cm", dpi = 600)


########## Count telomere-length observations per sample and haplotype #########
tlens_chr_filtered <- telogator2_5 %>%
  filter(!chr %in% c("chr13p", "chr14p", "chr15p", "chr21p", "chr22p", "chr15q"))    # 7735 rows

write_xlsx(tlens_chr_filtered, "tlens_chr_filtered.xlsx")


########Per sample
ChArmTL_counts_by_sample <- tlens_chr_filtered %>%
  count(sample, name = "n_ChArmTLs") %>%
  mutate(raw_family_role = tolower(sub(".*_", "", sample)), 
         family_role = case_when(raw_family_role == "fa" ~ "Father",
                                 raw_family_role == "mo" ~ "Mother",
                                 raw_family_role == "s1" ~ "Sibling",
                                 raw_family_role == "p1" ~ "ASD",
                                 TRUE ~ raw_family_role),
         family_role = factor(family_role, levels = c("Father","Mother","Sibling","ASD"))) %>% 
  select(-raw_family_role)

write_xlsx(ChArmTL_counts_by_sample, "ChArmTL_counts_by_sample.xlsx")

####Mean values and median
ChArmTL_counts_by_sample %>%
  summarise(n = n(),
            mean = sprintf("%.2f", mean(n_ChArmTLs, na.rm = TRUE)),
            median = median(n_ChArmTLs))

ChArmTL_counts_by_sample %>%
  group_by(family_role) %>%
  summarise(n = n(),
            mean = sprintf("%.2f", mean(n_ChArmTLs, na.rm = TRUE)),
            median = median(n_ChArmTLs))  


####Test whether there are differences between groups
ChArmTL_counts_by_sample %>%
  group_by(family_role) %>%
  shapiro_test(n_ChArmTLs)    # ASD and Father are non-normal distributions
leveneTest(n_ChArmTLs ~ family_role, data = ChArmTL_counts_by_sample)   #Equal variances
kruskal.test(n_ChArmTLs ~ family_role, data = ChArmTL_counts_by_sample)    #Test whether there are differences between groups
dunn_test(ChArmTL_counts_by_sample, n_ChArmTLs ~ family_role, p.adjust.method = "fdr")     # Post hoc test

p4 <- ggplot(data = ChArmTL_counts_by_sample, aes(x = family_role, y = n_ChArmTLs, fill = family_role))+
  geom_boxplot(linewidth = 0.3, median.linewidth = 0.3, width = 0.5, outlier.alpha = 0, show.legend = FALSE)+
  scale_fill_manual(values = c("#FCD28c","#A1B5A3","#82ADD0","#EB6368"))+
  scale_y_continuous(expand = c(0,0) ,limits = c(45,75), breaks = c(55,65,75))+
  labs(x = NULL, y = 'Number of ChArmTLs')+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.5, b = 0.25, l = 0.5, r = 0.25, unit = "cm"))+
  geom_signif(y_position = c(73,69), xmin = c(1,3), xmax = c(2,4), annotations = c("ns","ns"),
              tip_length = 0.01, size = 0.2, textsize = 1.5)

# ggplot_build(p4)$data[[1]] # Check the upper whisker endpoint (ymax) to determine geom_signif y_position


########Offspring haplotypes
ChArmTL_counts_by_offsring <- tlens_chr_filtered %>%
  filter(family_role %in% c("Sibling", "ASD")) %>%
  count(sample, parent_of_origin, name = "n_ChArmTLs", .drop = TRUE) %>%
  mutate(raw_family_role = tolower(sub(".*_", "", sample)),
         family_role = case_when(raw_family_role == "s1" ~ "Sib",
                                 raw_family_role == "p1" ~ "ASD",
                                 TRUE ~ raw_family_role),
         haplotype = paste(family_role, parent_of_origin, sep = "-"),
         haplotype = factor(haplotype, levels = c("Sib-Pat","Sib-Mat","ASD-Pat","ASD-Mat"))) %>%
  select(-raw_family_role)

write_xlsx(ChArmTL_counts_by_offsring, "ChArmTL_counts_by_offsring.xlsx")

####Mean values and median
ChArmTL_counts_by_offsring %>%
  summarise(n = n(),
            mean = sprintf("%.2f", mean(n_ChArmTLs, na.rm = TRUE)),
            median = median(n_ChArmTLs))

ChArmTL_counts_by_offsring %>%
  group_by(haplotype) %>%
  summarise(n = n(),
            mean = sprintf("%.2f", mean(n_ChArmTLs, na.rm = TRUE)),
            median = median(n_ChArmTLs))  

####Test whether there are differences between groups
ChArmTL_counts_by_offsring %>%
  group_by(haplotype) %>%
  shapiro_test(n_ChArmTLs)    # ASD-Pat and ASD-Mat are non-normal distributions
leveneTest(n_ChArmTLs ~ haplotype, data = ChArmTL_counts_by_offsring)     #Equal variances
kruskal.test(n_ChArmTLs ~ haplotype, data = ChArmTL_counts_by_offsring)    #Test whether there are differences between groups
dunn_test(ChArmTL_counts_by_offsring, n_ChArmTLs ~ haplotype, p.adjust.method = "fdr")    # Post hoc test

p5 <- ggplot(data = ChArmTL_counts_by_offsring, aes(x = haplotype, y = n_ChArmTLs, fill = haplotype))+
  geom_boxplot(linewidth = 0.3, median.linewidth = 0.3, width = 0.5, outlier.alpha = 0, show.legend = FALSE)+
  scale_fill_manual(values = c("#82ADD0","#82ADD080","#EB6368","#EB636880"))+
  scale_y_continuous(expand = c(0,0) ,limits = c(20,40), breaks = c(25,30,35,40))+
  labs(x = NULL, y = 'Number of ChArmTLs')+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.5, b = 0.25, l = 0.25, r = 0.5, unit = "cm"))+
  geom_signif(y_position = c(37,37,38,39), xmin = c(1,3,1,2), xmax = c(2,4,3,4), annotations = c("ns","ns","ns","ns"),
              tip_length = 0.01, size = 0.2, textsize = 1.5)

# ggplot_build(p5)$data[[1]]


######## Parental transmission status determined by Sibling
ChArmTL_counts_by_parents_to_Sib <- tlens_chr_filtered %>%
  filter(family_role %in% c("Father","Mother"),
         trans_status_Sib %in% c("T-Sib", "U-Sib")) %>%
  count(sample, trans_status_Sib, name = "n_ChArmTLs", .drop = TRUE) %>%
  mutate(raw_family_role = tolower(sub(".*_", "", sample)),
         family_role = case_when(raw_family_role == "fa" ~ "F",
                                 raw_family_role == "mo" ~ "M",
                                 TRUE ~ raw_family_role),
         haplotype = paste(family_role, trans_status_Sib, sep = ""),
         haplotype = factor(haplotype, levels = c("FT-Sib", "FU-Sib", "MT-Sib", "MU-Sib"))) %>%
  select(-raw_family_role)

write_xlsx(ChArmTL_counts_by_parents_to_Sib, "ChArmTL_counts_by_parents_to_Sib.xlsx")

#### Mean values and median
ChArmTL_counts_by_parents_to_Sib %>%
  summarise(n = n(),
            mean = sprintf("%.2f", mean(n_ChArmTLs, na.rm = TRUE)),
            median = median(n_ChArmTLs))

ChArmTL_counts_by_parents_to_Sib %>%
  group_by(haplotype) %>%
  summarise(n = n(),
            mean = sprintf("%.2f", mean(n_ChArmTLs, na.rm = TRUE)),
            median = median(n_ChArmTLs))

#### Test whether there are differences between groups
ChArmTL_counts_by_parents_to_Sib %>%
  group_by(haplotype) %>%
  shapiro_test(n_ChArmTLs)    # FU-Sib is non-normal distributions

leveneTest(n_ChArmTLs ~ haplotype, data = ChArmTL_counts_by_parents_to_Sib)    # Equal variances

kruskal.test(n_ChArmTLs ~ haplotype, data = ChArmTL_counts_by_parents_to_Sib)    # Test whether there are differences between groups
dunn_test(ChArmTL_counts_by_parents_to_Sib, n_ChArmTLs ~ haplotype, p.adjust.method = "fdr")    # Post hoc test

p6 <- ggplot(data = ChArmTL_counts_by_parents_to_Sib, aes(x = haplotype, y = n_ChArmTLs, fill = haplotype))+
  geom_boxplot(linewidth = 0.3, median.linewidth = 0.3, width = 0.5, outlier.alpha = 0, show.legend = FALSE)+
  scale_fill_manual(values = c("#FCD28c","#FCD28c80","#A1B5A3","#A1B5A380"))+
  scale_y_continuous(expand = c(0,0) ,limits = c(20,40), breaks = c(25,30,35,40))+
  labs(x = NULL, y = 'Number of ChArmTLs')+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.5, r = 0.25, unit = "cm"))+
  geom_signif(y_position = c(36,37,38,39), xmin = c(1,3,1,2), xmax = c(2,4,3,4), annotations = c("ns","ns","ns","ns"),
              tip_length = 0.01, size = 0.2, textsize = 1.5)

# ggplot_build(p6)$data[[1]]


######## Parental transmission status determined by ASD
ChArmTL_counts_by_parents_to_ASD <- tlens_chr_filtered %>%
  filter(family_role %in% c("Father", "Mother"),
         trans_status_ASD %in% c("T-ASD", "U-ASD")) %>%
  count(sample, trans_status_ASD, name = "n_ChArmTLs", .drop = TRUE) %>%
  mutate(raw_family_role = tolower(sub(".*_", "", sample)),
         family_role = case_when(raw_family_role == "fa" ~ "F",
                                 raw_family_role == "mo" ~ "M",
                                 TRUE ~ raw_family_role),
         haplotype = paste(family_role, trans_status_ASD, sep = ""),
         haplotype = factor(haplotype, levels = c("FT-ASD","FU-ASD","MT-ASD","MU-ASD"))) %>%
  select(-raw_family_role)

write_xlsx(ChArmTL_counts_by_parents_to_ASD, "ChArmTL_counts_by_parents_to_ASD.xlsx")

#### Mean values and median
ChArmTL_counts_by_parents_to_ASD %>%
  summarise(n = n(),
            mean = sprintf("%.2f", mean(n_ChArmTLs, na.rm = TRUE)),
            median = median(n_ChArmTLs))

ChArmTL_counts_by_parents_to_ASD %>%
  group_by(haplotype) %>%
  summarise(n = n(),
            mean = sprintf("%.2f", mean(n_ChArmTLs, na.rm = TRUE)),
            median = median(n_ChArmTLs))

#### Test whether there are differences between groups
ChArmTL_counts_by_parents_to_ASD %>%
  group_by(haplotype) %>%
  shapiro_test(n_ChArmTLs)    # FT-ASD, FU-ASD are non-normal distributions

leveneTest(n_ChArmTLs ~ haplotype, data = ChArmTL_counts_by_parents_to_ASD)    #Equal variances

kruskal.test(n_ChArmTLs ~ haplotype, data = ChArmTL_counts_by_parents_to_ASD)    # Test whether there are differences between groups
dunn_test(ChArmTL_counts_by_parents_to_ASD, n_ChArmTLs ~ haplotype, p.adjust.method = "fdr")    # Post hoc test

p7 <- ggplot(data = ChArmTL_counts_by_parents_to_ASD, aes(x = haplotype, y = n_ChArmTLs, fill = haplotype))+
  geom_boxplot(linewidth = 0.3, median.linewidth = 0.3, width = 0.5, outlier.alpha = 0, show.legend = FALSE)+
  scale_fill_manual(values = c("#F9D5B2","#F9D5B280","#87BEA7","#87BEA780"))+
  scale_y_continuous(expand = c(0,0) ,limits = c(20,40), breaks = c(25,30,35,40))+
  labs(x = NULL, y = 'Number of ChArmTLs')+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks =  element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t = 0.25, b = 0.5, l = 0.25, r = 0.5, unit = "cm"))+
  geom_signif(y_position = c(36,37,38,39), xmin = c(1,3,1,2), xmax = c(2,4,3,4), annotations = c("ns","ns","ns","ns"),
              tip_length = 0.01, size = 0.2, textsize = 1.5)

# ggplot_build(p7)$data[[1]]

p4_5_6_7 <- plot_grid(p4, p5, p6, p7, ncol = 2, labels = c('a','b','c','d'), label_size = 8)

ggsave(filename = "Supplementary Figure 3.pdf", plot = p4_5_6_7, width = 12, height = 10, units = "cm", dpi = 600)

