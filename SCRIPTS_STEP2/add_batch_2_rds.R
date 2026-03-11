# Load necessary packages
library(Seurat)
library(dplyr)

# Read your RDS file
#RDS_FILE <- "./obj_lineage_Non_lymphoid_cleaned.rds"
#RDS_FILE    <- "../MERGED_PIPELINE_OUTPUT_V5/object_final_annotated_panGI.rds"
RDS_FILE    <- "./obj_lineage_Lymphoid_cleaned.rds"
obj <- readRDS(RDS_FILE)

obj@meta.data <- obj@meta.data %>%
  mutate(
    Lane = as.character(Lane),  # ensure character type
    Day = case_when(
      Lane %in% c("1", "2") ~ "Day1",
      Lane %in% c("3", "4") ~ "Day2",
      Lane %in% c("5", "6") ~ "Day3",
      TRUE ~ NA_character_
    )
  )


#saveRDS(obj, file = "./obj_lineage_Non_lymphoid_withDay.rds")
#saveRDS(obj, file = "./object_final_annotated_panGI_withDay.rds")
saveRDS(obj, file = "./obj_lineage_Lymphoid_withDay.rds")