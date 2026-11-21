########################## set the work directory ##########################
rm(list=ls())
setwd("Path/to/your_workdir")    # the folder path where the "merged_tlens_by_allele.tsv" file is stored
getwd()


############################ Load R package and data ###########################
install.packages(c("readr", "dplyr", "stringr","writexl"))

library(readr)
library(dplyr)
library(stringr)
library(writexl)

telogator2_0 <- read_delim("merged_tlens_by_allele.tsv", delim = "\t", escape_double = FALSE, 
                           col_types = cols(position = col_character(), 
                                            allele_id = col_character(), TL_p75 = col_number(), 
                                            read_lengths = col_character(), read_mapq = col_character(), 
                                            tvr_len = col_number()), trim_ws = TRUE)    # 11277 rows


#################### Extract columns for subsequent analysis ####################
colnames(telogator2_0)

telogator2_0 <- telogator2_0 %>% 
  select(Sample, `#chr`, allele_id, TL_p75, read_TLs, tvr_len) %>%
  rename(chr = `#chr`) %>%
  mutate(number_of_supporting_reads = str_count(read_TLs, ",") + 1)


############################### QC and filtering ###############################
######## excluding inaccurate estimations
telogator2_1 <- telogator2_0 %>%
  filter(!grepl("chrUq|chrUp", chr),  # unmapped, 11277 rows ---> 11277 rows
         !grepl(",|，", chr),         # ambiguous alignments, 11277 rows ---> 10255 rows
         !grepl("i$", allele_id),     # putative interstitial telomeric sequences, 10255 rows ---> 10052 rows
         tvr_len != 0)                # tvr_len was zero, 10255 rows ---> 9978 rows


######## excluding sex chromosomes
telogator2_2 <- telogator2_1[!grepl("X|Y", telogator2_1$chr), ]    # 9978 rows ---> 9599 rows


######## excluding outliers
summary(telogator2_2$TL_p75)

hist(telogator2_2$TL_p75, prob = T)    #Histogram, non-normal distribution
lines(density(na.omit(telogator2_2$TL_p75)), col = 'blue')   #Density curve

#### remove outliers (IQR)
telogator2_3 <- telogator2_2 %>%
  mutate(
    raw_family_role = tolower(sub(".*_", "", sample)),              # automatically extract "family_role"
    family_role = case_when(raw_family_role == "fa" ~ "Father",
                            raw_family_role == "mo" ~ "Mother",
                            raw_family_role == "s1" ~ "Sibling",
                            raw_family_role == "p1" ~ "ASD",
                            TRUE ~ raw_family_role),
    QC_status = if_else(TL_p75 <= 100, "outlier", "passed")) %>%    # Initial absolute threshold filtering: values <= 100 are directly marked as outlier
  group_by(chr) %>%
  mutate(
    Q1 = quantile(TL_p75[QC_status == "passed"], 0.25, na.rm = TRUE),
    Q3 = quantile(TL_p75[QC_status == "passed"], 0.75, na.rm = TRUE),
    IQR_val = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR_val,
    upper_bound = Q3 + 1.5 * IQR_val,
    QC_status = case_when(
      QC_status == "outlier" ~ "outlier",
      TL_p75 < lower_bound | TL_p75 > upper_bound ~ "outlier",
      TRUE ~ "passed")) %>%
  ungroup() %>%
  select(-raw_family_role, -Q1, -Q3, -IQR_val, -lower_bound, -upper_bound)

table(telogator2_3$QC_status, telogator2_3$family_role)    #303 outliers [72(ASD) + 81(Father) + 63(Mother) + 87(Sibling)] + 9296 passed

telogator2_4 <- subset(telogator2_3, QC_status == "passed")    # 9599 rows ---> 9296 rows
write_xlsx(telogator2_4, "tlens_outlier_filtered.xlsx")

####Distribution and normality check after removing outliers
summary(telogator2_4$TL_p75)

hist(telogator2_4$TL_p75, prob = T)  #Histogram, approximately normal distribution
lines(density(na.omit(telogator2_4$TL_p75)), col = 'blue')   #Density curve

