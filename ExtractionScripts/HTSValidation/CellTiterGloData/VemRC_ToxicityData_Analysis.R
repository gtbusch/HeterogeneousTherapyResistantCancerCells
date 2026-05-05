library(tidyverse)
library(dplyr)
library(here)
library(readxl)

##

## Task 0 - Load the Data

current_path <- here()

parent_path <- dirname(dirname(dirname(current_path)))

new_path <- file.path(parent_path, "RawData", "HTSValidation", "CellTiterGloData")

data_path <- file.path(parent_path, "ExtractedData", "HTSValidation", "CellTiterGloData")

sara_data = read_excel(file.path(new_path,"Saracatinib_ToxicityData.xlsx"))
jq1_data = read_excel(file.path(new_path,"JQ1_ToxicityData.xlsx"))
otx_data = read_excel(file.path(new_path,"OTX015_ToxicityData.xlsx"))
das_data = read_excel(file.path(new_path,"Dasatinib_ToxicityData.xlsx"))
clo_data = read_excel(file.path(new_path,"Clofarabine_ToxicityData.xlsx"))
camp_data = read_excel(file.path(new_path,"Camptothecin_ToxicityData.xlsx"))
epi_data = read_excel(file.path(new_path,"EpirubicinHCl_ToxicityData.xlsx"))
s3i_data = read_excel(file.path(new_path,"S3I201_ToxicityData.xlsx"))
jq1_coculture_data = read_excel(file.path(new_path, "JQ1_CoCulture_ToxicityData.xlsx"))

all_data = rbind(jq1_data, clo_data, epi_data, s3i_data, sara_data, otx_data, das_data, camp_data)

saveRDS(all_data, file.path(data_path, "validationdata_20240403.rds")) 

revision_round1data = read_excel(file.path(new_path, "AnnexinRound1_ToxicityData.xlsx"))

revision_round1data = rbind(revision_round1data)

saveRDS(revision_round1data, file.path(data_path, "revision_round1data.rds"))

revision_round2data = read_excel(file.path(new_path, "AnnexinRound2_ToxicityData.xlsx"))

revision_round2data = rbind(revision_round2data)

saveRDS(revision_round2data, file.path(data_path, "revision_round2data.rds"))

jq1_coculture_data = rbind(jq1_coculture_data)

saveRDS(jq1_coculture_data, file.path(data_path, "jq1_coculture_data.rds"))

revision_round3data = read_excel(file.path(new_path, "AnnexinRound3_ToxicityData.xlsx"))

revision_round3data = rbind(revision_round3data)

saveRDS(revision_round3data, file.path(data_path, "revision_round3data.rds"))
