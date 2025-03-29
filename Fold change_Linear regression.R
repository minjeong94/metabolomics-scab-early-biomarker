library(glmnet)
library(tidyverse)
library(e1071)
# Set working directory
setwd("C:/Users/Min Jeong Kang/Desktop")

#List of Biomarkers identified by LR-L2
combined_1_2dpi <- read_csv("Biomarkers_1_2dpi_1229.csv")
C30_list_1_2dpi <- Combined_1_2dpi %>%
  filter(dataset == "C30")
hilic_list_1_2dpi<- Combined_1_2dpi %>%
  filter(dataset == "hilic")


Biomarkers_1_2dpi <-  union(C30_list_1_2dpi$metabolite, hilic_list_1_2dpi$metabolite)

#-----------------------------------------------------------------------------------------
#-------Check the significance with linear regression-------------------------------------
#-----------------------------------------------------------------------------------------


library(tidyverse)
library(mlogit)
library(fixest)

# Set working directory
setwd("C:/Users/Min Jeong Kang/Desktop")

# Read the data
data <- read_csv("scab_obj2_final.csv")

data_sig <- data %>%
  filter(hpi !="0day") %>%
  mutate(stage = case_when(
    hpi %in% c("1day", "2day") ~ "1st",
    hpi %in% c("3day", "4day") ~ "2nd",
    hpi %in% c("5day", "7day") ~ "3rd"
  ))

combined_list_1_2dpi <- Biomarkers_1_2dpi


#--Perform the linear regression with the biomarkers selected by LR-L2--

# Change the stage and list of biomarkers for 2nd (3-4DPI) and 3rd(5-7 DPI) inoculation times----
linear_result <- lapply(
  combined_list_1_2dpi,
  function(x) {
    fml <- formula(paste0(x, " ~ treatment"))
    res <- lm(fml, filter(data_sig, treatment != "control", stage == "1st"))
    final_tibble <- tibble(
      metabolite = x,
      beta = summary(res)$coefficients[2,1],
      p_value = summary(res)$coefficients[2,4],
      significance = ifelse(p_value > 0.05, FALSE, TRUE)
    )
  }
) %>% bind_rows(.)


print(linear_result, n =1000)
sig_biomarkers <- linear_result %>%
  filter(significance == "TRUE")


#--------------------------------------------------------------------------
#---------------------Fold change------------------------------------------
#--------------------------------------------------------------------------
# Change the stage and the list of biomarkers for 2nd (3-4DPI) and 3rd(5-7 DPI) inoculation times-------
data_sig_1st <- data_sig %>% filter(stage == "1st")
filtered_data_1st <- data_sig_1st[, c("sample", "treatment", "hpi","bio_rep", combined_list_1_2dpi), drop = FALSE]
glimpse(filtered_data_1st)

# Load dplyr library
library(dplyr)


#ADP:Valine <---- This part should be changed depending on the list of biomarkers
filtered_Mdata_1st <- filtered_data_1st %>%
  group_by(sample, bio_rep, treatment, hpi) %>%
  summarise(across(ADP:Valine , mean, na.rm = TRUE), .groups = "drop")



glimpse(filtered_Mdata_1st)


fold_change <- function(data) {
  
  res <- lapply(
    combined_list_1_2dpi,
    function(x) {
      data <- data %>% 
        select(sample, bio_rep, treatment, hpi, x) %>% 
        group_by(bio_rep, hpi)  %>%
        summarise(
          fold_change_res = log2(get(x)[treatment == "resistant"]/get(x)[treatment == "control"]),
          fold_change_sus = log2(get(x)[treatment == "susceptible"]/get(x)[treatment == "control"])
        ) %>% 
        bind_cols(tibble(metabolite = rep(x, 6)))
    }
  ) %>% bind_rows(.)
  
}

foldchange_1st <- fold_change(filtered_Mdata_1st)

foldchange_1st_av <- foldchange_1st %>%
  group_by(metabolite) %>%
  summarize(
    avg_fold_change_res = mean(fold_change_res, na.rm = TRUE),
    avg_fold_change_sus = mean(fold_change_sus, na.rm = TRUE)
  )
sig_biomarkers_1_2dpi

Combined_1_2dpi

Biomarkers_1_2dpi_fc_coef <- Combined_1_2dpi %>%
  left_join(foldchange_1st_av, by = "metabolite")
print(Biomarkers_1_2dpi_fc_coef, n=39)


# Perform a left join to add p_value and significance columns
Biomarkers_1_2dpi_fc_coef_pvalue <- Biomarkers_1_2dpi_fc_coef %>%
  left_join(sig_biomarkers_1_2dpi %>% select(metabolite, p_value, significance), by = "metabolite") %>%
  # Replace NA in significance with "NS"
  mutate(significance = ifelse(is.na(significance), "NS", as.character(significance)))

# View the result
print(Biomarkers_1_2dpi_fc_coef_pvalue, n=39)


# Export the data into xlsx, Save the file accordingly
library(writexl)
library(openxlsx)

write_xlsx(Biomarkers_1_2dpi_fc_coef_pvalue, path = ".xlsx")
write.xlsx(Biomarkers_1_2dpi_fc_coef_pvalue, file = ".xlsx")

