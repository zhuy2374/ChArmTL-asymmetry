########################## set the work directory ##########################
rm(list=ls())
setwd("Path/to/your_workdir")    # the folder path where the "merged_tlens_by_allele.tsv" file is stored
getwd()


############################ Load R package and data ###########################
install.packages(c("readr", "dplyr", "writexl"))

library(readr)
library(dplyr)
library(writexl)

telogator2_0 <- read_delim("merged_tlens_by_allele.tsv", delim = "\t", escape_double = FALSE, 
                           col_types = cols(position = col_character(), 
                                            allele_id = col_character(), 
                                            read_lengths = col_character()),
                           trim_ws = TRUE)    # 11277 rows


#################### Extract columns for subsequent analysis ####################
colnames(telogator2_0)

telogator2_0 <- telogator2_0 %>% 
  select(Sample, `#chr`, allele_id, TL_p75, read_TLs, tvr_len)    # select "Sample", "#chr", "allele_id", "TL_p75", "read_TLs", "tvr_len"


############################### QC and filtering ###############################
######## excluding inaccurate estimations
#### unmapped
telogator2_1 <- telogator2_0[!grepl("chrUq|chrUp", telogator2_0$`#chr`), ]    # 11277 rows ---> 11277 rows

#### ambiguous alignments
telogator2_1 <- telogator2_1[!grepl(",|，", telogator2_1$`#chr`), ]    # 11277 rows ---> 10255 rows

#### putative interstitial telomeric sequences
telogator2_1 <- telogator2_1[!grepl("i$", telogator2_1$allele_id), ]    # 10255 rows ---> 10052 rows

#### tvr_len was zero
telogator2_1 <- telogator2_1[telogator2_1$tvr_len != 0, ]    # 10255 rows ---> 9978 rows


######## excluding sex chromosomes
telogator2_2 <- telogator2_1[!grepl("X|Y", telogator2_1$`#chr`), ]    # 9978 rows ---> 9599 rows


######## excluding outliers
summary(telogator2_2$TL_p75)    #Descriptive analysis of telomere length

hist(telogator2_2$TL_p75, prob = T)    #Histogram, non-normal distribution
lines(density(na.omit(telogator2_2$TL_p75)), col = 'blue')   #Density curve

#### remove outliers (IQR), separating offspring and parents
telogator2_3 <- telogator2_2 %>%
  mutate(
    clinical_group = sub(".*_", "", sample),    # automatically extract "family_role"
    QC_status = if_else(TL_p75 <= 100, "outlier", "passed"),    # Initial absolute threshold filtering: values <= 100 are directly marked as outliers
    role_group = case_when(    # classify the family roles into groups to facilitate subsequent batch calculations
      clinical_group %in% c("p1", "s1") ~ "offspring",
      clinical_group %in% c("fa", "mo") ~ "parent",
      TRUE ~ "other")) %>%
  group_by(role_group) %>%     # group the data by role_group for calculation
  mutate(    # calculate the Q1 and Q3 values for each group dynamically (based solely on the current "passed" data
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
  select(-role_group, -Q1, -Q3, -IQR_val, -lower_bound, -upper_bound)    # remove the auxiliary columns generated during the intermediate process

table(telogator2_3$QC_status)    #242 outliers + 9357 passed

telogator2_4 <- subset(telogator2_3, QC_status == "passed")    # 9599 rows ---> 9357 rows
write_xlsx(telogator2_4, "tlens_outlier_filtered.xlsx")

####Distribution and normality check after removing outliers
summary(telogator2_4$TL_p75)

hist(telogator2_4$TL_p75, prob = T)  #Histogram, approximately normal distribution
lines(density(na.omit(telogator2_4$TL_p75)), col = 'blue')   #Density curve

